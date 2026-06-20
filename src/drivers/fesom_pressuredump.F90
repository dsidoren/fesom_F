program fesom_pressuredump
    ! M2.1 pressure/EOS + M2.2 PGF + M2.3 vel_rhs operator byte-gate driver. Loads pi
    ! 1-rank with the SAME rotation as the FESOM2 pi run (alpha/beta/gamma=50/15/-90,
    ! cyclic 360), computes geometry (incl. M2.3 coriolis), builds the ALE depths
    ! Z_3d_n/zbar_3d_n + node thickness hnode (linfs full cells; M1.2/M1.3 proved
    ! them), PRESCRIBES analytic T/S + UV/eta_n/UV_rhsAB fields (identical formulas to
    ! the FESOM2 oracle src/fesom_pressure_dump.F90, from the byte-identical coords),
    ! then runs oce_pressure_bv (raw + smoothed N^2) -> oce_pgf -> oce_dyn_velrhs and
    ! dumps density_m_rho0 / hpressure / bvfreq / pgf_x/pgf_y / coriolis / uv_rhs* in
    ! the oracle's FADVHDMP format. tools/pressure_diff.py compares for max|delta|=0.
    !
    !   FESOM3_MESH_DIR      mesh dir   (default: pi)
    !   FESOM3_PRESSURE_OUT  out path   (default: pressure_f3.bin)
    use mpi
    use mod_precision,    only: WP, MP
    use mod_constants,    only: density_0
    use mod_param_phys,   only: N2smth_h
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex
    use mod_mesh_read,    only: read_mesh
    use mod_mesh_areas,   only: compute_geometry
    use mod_dyn,          only: t_dyn
    use oce_pressure_bv,  only: pressure_bv
    use oce_pgf,          only: pressure_force_4_linfs_fullcell
    use oce_dyn_velrhs,   only: compute_vel_rhs
    use oce_dyn_visc,     only: viscosity_filter
    use oce_dyn_ivertvisc, only: impl_vert_visc_ale
    use oce_ssh_rhs,      only: init_stiff_mat_ale, compute_ssh_rhs_ale
    use oce_ssh_solve,    only: solve_ssh_ale
    use oce_ale,          only: update_vel, compute_hbar_ale, update_eta_n, vert_vel_ale
    use mod_param_phys,   only: alpha, theta
    use mod_advhor_dump,  only: advhor_dump_open, advhor_dump_close, wr_r1, wr_r2, wr_r3
    implicit none

    ! Pinned shared constant (NOT pi's namelist dt=86400/36=2400 s): the vel_rhs gate
    ! only needs BOTH codes to use the same dt in the dt*(...)/elem_area scaling.
    real(kind=WP), parameter :: dt_velrhs = 1800.0_WP
    ! M2.6 SSH stiffness dt: the FESOM2 oracle built ssh_stiff at ocean_setup (line 140)
    ! with the pi NAMELIST dt = 86400/step_per_day, step_per_day=36 -> 2400 s, BEFORE the
    ! shim overrides dt to dt_velrhs. So the stiffness matrix uses this value; compute it
    ! the SAME way FESOM2 does (gen_model_setup.F90:92) — NOT a 2400.0 literal — so the
    ! -no-prec-div bits agree (the only dt-dependent part of M2.6).
    real(kind=WP), parameter :: dt_ssh = 86400.0_WP / real(36, WP)

    character(len=512) :: mesh_dir, out_path
    type(t_partit)      :: partit
    type(t_mesh)        :: mesh
    integer :: nsw, n, nz, nl, u, nzmin, nzmax, e
    real(kind=WP) :: lon, lat
    real(kind=MP) :: zbar_srf, zbar_bot
    real(kind=WP), allocatable :: temp(:,:), salt(:,:), density_ref(:,:)
    real(kind=WP), allocatable :: density(:,:), hpressure(:,:), bvfreq(:,:)
    real(kind=WP), allocatable :: bvfreq_raw(:,:), pgf_x(:,:), pgf_y(:,:)
    ! --- M2.3 vel_rhs (Coriolis + AB2 + PGF + SSH gradient) + M2.4 momentum advection ---
    type(t_dyn) :: dyn
    real(kind=WP), allocatable :: uv_in(:,:,:), uv_rhsAB_prev(:,:,:)
    real(kind=WP), allocatable :: uv_rhs_eul(:,:,:), uv_rhs_ab2(:,:,:), uv_rhsAB_cor(:,:,:)
    real(kind=WP), allocatable :: uvnode_rhs_dump(:,:,:)   ! M2.4 momadv intermediate
    ! --- M2.4 biharmonic viscosity (opt_visc=7, visc_filt_bidiff) ---
    real(kind=WP), allocatable :: uv_rhs_visc(:,:,:), u_c_dump(:,:), v_c_dump(:,:)
    integer       :: ed, el(2), ndu, ng0, ng1, ng2
    real(kind=WP) :: du, dumax
    ! --- M2.5 implicit vertical viscosity (TDMA / Thomas solve, impl_vert_visc_ale) ---
    real(kind=WP), allocatable :: Av(:,:), stress_surf(:,:), uv_rhs_ivv(:,:,:)
    integer       :: nwp, nwm
    ! --- M2.6 SSH (stiffness matrix + ssh_rhs + CG solve) ---
    real(kind=WP), allocatable :: ssh_diag(:), ssh_Aeta(:)
    integer       :: n_ssh_iter, row, ni, ne2
    real(kind=WP) :: sum_ssh_rhs, resid_inf
    ! --- M2.7 ALE (linfs) velocity / SSH / thickness-W update ---
    ! Saved copies of the prescribed inputs the M2.7 kernels overwrite in place
    ! (eta_n by update_eta_n; w_e/w_i by compute_Wvel_split) so the M2.3/M2.4/M2.5
    ! INPUT dump records still echo the prescription, not the M2.7 output.
    real(kind=WP), allocatable :: hbar_in(:), eta_n_in(:), w_e_in(:,:), w_i_in(:,:)
    real(kind=WP), allocatable :: uv_upd(:,:,:)
    integer       :: n_cflsplit

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call get_environment_variable('FESOM3_PRESSURE_OUT', out_path)
    if (len_trim(out_path) == 0) out_path = 'pressure_f3.bin'

    call par_init(partit)
    if (partit%npes /= 1) then
        if (partit%mype == 0) write(*,'(a)') 'fesom_pressuredump: requires 1 rank'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype); error stop 1
    end if

    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    nl = mesh%nl
    write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') 'fesom_pressuredump: nod2D=', mesh%nod2D, &
        ' elem2D=', mesh%elem2D, ' edge2D=', mesh%edge2D, ' nl=', nl, ' CW swaps=', nsw

    ! hnode (per-node layer thickness, linfs full cells): hnode(nz,n)=zbar(nz)-zbar(nz+1)
    ! for nz in [ulevels_nod2D(n), nlevels_nod2D(n)-1]. Proven by the M1.3 gate.
    allocate(mesh%hnode(nl-1, mesh%nod2D)); mesh%hnode = 0.0_MP
    do n = 1, mesh%nod2D
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n)-1
            mesh%hnode(nz, n) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
    end do

    ! helem (per-element layer thickness, linfs full cells): helem(nz,e)=zbar(nz)-zbar(nz+1)
    ! for nz in [ulevels(e), nlevels(e)-1]; the element analog of hnode (M1 advection gate).
    ! zbar_e_bot (partial-cell bottom depth) = zbar(nlevels(e)) for full cells (FESOM2
    ! init_bottom_elem_thickness, use_partial_cell=.false.). Both consumed by M2.5
    ! impl_vert_visc_ale to rebuild the per-column zbar_n/Z_n.
    allocate(mesh%helem(nl-1, mesh%elem2D)); mesh%helem = 0.0_MP
    allocate(mesh%zbar_e_bot(mesh%elem2D)); mesh%zbar_e_bot = 0.0_MP
    do e = 1, mesh%elem2D
        do nz = mesh%ulevels(e), mesh%nlevels(e)-1
            mesh%helem(nz, e) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
        mesh%zbar_e_bot(e) = mesh%zbar(mesh%nlevels(e))
    end do

    ! zbar_3d_n / Z_3d_n: per-node ALE interface/mid depths, built exactly as FESOM2
    ! init_ale (oce_ale.F90:531-566) at the initial state (no cavity, full cells).
    ! Proven field-by-field by the M1.2 gate.
    allocate(mesh%zbar_3d_n(nl, mesh%nod2D), mesh%Z_3d_n(nl-1, mesh%nod2D))
    mesh%zbar_3d_n = 0.0_MP; mesh%Z_3d_n = 0.0_MP
    do n = 1, mesh%nod2D
        nzmin = mesh%ulevels_nod2D(n)
        nzmax = mesh%nlevels_nod2D(n)
        zbar_srf = mesh%zbar(nzmin)
        zbar_bot = mesh%zbar(nzmax)
        mesh%zbar_3d_n(1:nzmin-1, n)       = mesh%zbar(1:nzmin-1)
        mesh%zbar_3d_n(nzmin, n)           = zbar_srf
        mesh%zbar_3d_n(nzmin+1:nzmax-1, n) = mesh%zbar(nzmin+1:nzmax-1)
        mesh%zbar_3d_n(nzmax, n)           = zbar_bot
        mesh%Z_3d_n(1:nzmin-1, n)          = mesh%Z(1:nzmin-1)
        mesh%Z_3d_n(nzmin, n)              = mesh%zbar_3d_n(nzmin,n)   + (mesh%zbar_3d_n(nzmin+1,n)-zbar_srf)/2
        mesh%Z_3d_n(nzmin+1:nzmax-2, n)    = mesh%Z(nzmin+1:nzmax-2)
        mesh%Z_3d_n(nzmax-1, n)            = mesh%zbar_3d_n(nzmax-1,n) + (zbar_bot-mesh%zbar_3d_n(nzmax-1,n))/2
    end do

    ! Prescribe analytic temperature + salinity (at nodes). MUST match the FESOM2
    ! oracle src/fesom_pressure_dump.F90 exactly:
    !   T(nz,n) = 12.0 + 8.0*cos(lat)*cos(lon) - 0.20*nz
    !   S(nz,n) = 34.5 + 0.5*sin(2*lon)*cos(lat) + 0.03*nz   (always > 0 for sqrt(s))
    ! Strong horizontal (cos lat/lon) + vertical (nz) structure exercises the EOS,
    ! the top-down hpressure integration, the N^2 difference, and the horizontal
    ! N^2 smoother (which only changes a horizontally-varying field).
    allocate(temp(nl-1, mesh%nod2D), salt(nl-1, mesh%nod2D), density_ref(nl-1, mesh%nod2D))
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        do nz = 1, nl-1
            temp(nz, n) = 12.0_WP + 8.0_WP*cos(lat)*cos(lon) - 0.20_WP*real(nz, WP)
            salt(nz, n) = 34.5_WP + 0.5_WP*sin(2.0_WP*lon)*cos(lat) + 0.03_WP*real(nz, WP)
        end do
    end do
    ! density_ref = density_0 everywhere (use_density_ref=.false. on pi; the
    ! density anomaly subtracts this ARRAY, not the scalar).
    density_ref = density_0

    allocate(density(nl-1, mesh%nod2D), hpressure(nl, mesh%nod2D), bvfreq(nl, mesh%nod2D))
    allocate(bvfreq_raw(nl, mesh%nod2D))

    ! Pass 1: horizontal N^2 smoothing OFF -> raw bvfreq (+ density_m_rho0, hpressure).
    ! The caller pre-zeros the outputs (pressure_bv leaves below-bottom entries as-is).
    N2smth_h = .false.
    density = 0.0_WP; hpressure = 0.0_WP; bvfreq = 0.0_WP
    call pressure_bv(temp, salt, density_ref, mesh, density, hpressure, bvfreq)
    bvfreq_raw = bvfreq

    ! Pass 2: horizontal N^2 smoothing ON -> smoothed bvfreq.
    N2smth_h = .true.
    bvfreq = 0.0_WP
    call pressure_bv(temp, salt, density_ref, mesh, density, hpressure, bvfreq)

    ! M2.2: hydrostatic PGF from the M2.1 hpressure (the smoothing pass leaves
    ! hpressure unchanged - it only touches bvfreq, L13). Caller pre-zeros pgf_x/pgf_y
    ! (oce_pgf writes only levels ule..nle; FESOM2 leaves below-bottom at alloc-zero).
    allocate(pgf_x(nl-1, mesh%elem2D), pgf_y(nl-1, mesh%elem2D))
    pgf_x = 0.0_WP; pgf_y = 0.0_WP
    call pressure_force_4_linfs_fullcell(hpressure, mesh, pgf_x, pgf_y)

    ! ============== M2.3 vel_rhs assembly + M2.4 momentum advection ==============
    ! Coriolis + AB2 + PGF + SSH-gradient + momentum advection (momadv_opt=2). Build a
    ! t_dyn, prescribe analytic UV (elements) / eta_n (nodes) / previous-step UV_rhsAB +
    ! w_e (nodes, the explicit vertical velocity momadv reads; all identical to the
    ! FESOM2 oracle, from the byte-identical coords), copy in the M2.2 pgf_x/pgf_y, and
    ! run oce_dyn_velrhs TWICE: lfirst=.true. (first Euler step, ff=1.0) then
    ! lfirst=.false. (AB2 step, ff=ab2=1.6). The oracle drives the REAL FESOM2
    ! compute_vel_rhs the same way (now with momadv_opt=2 -> momentum_adv_scalar).
    allocate(dyn%uv(2, nl-1, mesh%elem2D), dyn%uv_rhs(2, nl-1, mesh%elem2D))
    allocate(dyn%uv_rhsAB(1, 2, nl-1, mesh%elem2D))     ! (AB_order-1, 2, nl-1, elem2D)
    allocate(dyn%eta_n(mesh%nod2D))
    allocate(dyn%w_e(nl, mesh%nod2D))                   ! explicit vertical velocity (nodes)
    allocate(dyn%w_i(nl, mesh%nod2D))                   ! M2.5 implicit vertical velocity (nodes)
    allocate(dyn%work%uvnode_rhs(2, nl-1, mesh%nod2D))  ! momadv nodal scratch
    allocate(dyn%work%pgf_x(nl-1, mesh%elem2D), dyn%work%pgf_y(nl-1, mesh%elem2D))
    allocate(dyn%work%u_c(nl-1, mesh%elem2D), dyn%work%v_c(nl-1, mesh%elem2D))  ! visc scratch
    allocate(Av(nl, mesh%elem2D), stress_surf(2, mesh%elem2D))  ! M2.5 prescribed inputs
    dyn%AB_order          = 2
    dyn%momadv_opt        = 2     ! M2.4: enable momentum_adv_scalar (pi production value)
    ! M2.4 biharmonic viscosity (opt_visc=7) — pin the pi/reduced-M2 namelist values
    ! (visc_gamma0=0.003 OVERRIDES the t_dyn default 0.03; gamma_h=0 -> pure biharmonic).
    ! MUST match the FESOM2 oracle src/fesom_pressure_dump.F90.
    dyn%opt_visc          = 7
    dyn%visc_gamma0       = 0.003_WP
    dyn%visc_gamma1       = 0.1_WP
    dyn%visc_gamma2       = 0.285_WP
    dyn%visc_gamma0_h     = 0.0_WP
    dyn%visc_gamma1_h     = 0.0_WP
    dyn%w_e               = 0.0_WP
    dyn%w_i               = 0.0_WP
    dyn%work%uvnode_rhs   = 0.0_WP
    dyn%work%u_c          = 0.0_WP
    dyn%work%v_c          = 0.0_WP
    dyn%work%pgf_x = pgf_x        ! reuse the proven M2.2 PGF (exact copy)
    dyn%work%pgf_y = pgf_y

    ! analytic SSH + explicit vertical velocity w_e at nodes (w_e: sign varies in
    ! space AND depth via cos(0.3*nz) -> non-trivial w*du/dz; MUST match the oracle).
    ! M2.5 also prescribes the IMPLICIT vertical velocity w_i (nodes): a DISTINCT analytic
    ! formula (sin(lon)*cos(2*lat)*cos(0.25*nz)) so its sign varies in space AND depth ->
    ! exercises both the min(0,wu/wd) and max(0,wu/wd) upwind branches of the vertical-
    ! advection update in impl_vert_visc_ale. MUST match the FESOM2 oracle.
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        dyn%eta_n(n) = 0.5_WP*cos(lat)*sin(lon) + 0.3_WP*sin(2.0_WP*lat)
        do nz = 1, nl
            dyn%w_e(nz,n) = 1.0e-4_WP*sin(2.0_WP*lon)*cos(lat)*cos(0.3_WP*real(nz,WP))
            dyn%w_i(nz,n) = 2.0e-4_WP*sin(lon)*cos(2.0_WP*lat)*cos(0.25_WP*real(nz,WP))
        end do
    end do
    ! analytic element velocity + previous-step AB array (from element node-1 coords).
    ! Horizontal amplitude 2.00/1.50 m/s (was 0.10/0.08 at M2.3/M2.4): a strong-current
    ! stress test so the across-edge velocity jump |du| spans ALL THREE biharmonic-
    ! viscosity flow-aware branches of max(gamma0, max(gamma1*|du|, gamma2*|du|^2)). With
    ! the pi gammas the winners cross at |du|=gamma0/gamma1=0.03 (gamma0->gamma1) and
    ! |du|=gamma1/gamma2=0.35 (gamma1->gamma2), so max|du|~0.71 (~2x the upper crossover)
    ! exercises gamma0, gamma1 AND gamma2 with a clear margin (L11/L17). The depth terms
    ! (-0.005*nz / +0.004*nz) cancel in du=UV(el1)-UV(el2). M2.3/M2.4 re-gate with the new
    ! values (identical formula on both sides). MUST match the FESOM2 oracle.
    do e = 1, mesh%elem2D
        lon = mesh%coord_nod2D(1, mesh%elem2D_nodes(1,e))
        lat = mesh%coord_nod2D(2, mesh%elem2D_nodes(1,e))
        do nz = 1, nl-1
            dyn%uv(1,nz,e) =  2.00_WP*cos(lat)*sin(lon)        - 0.005_WP*real(nz,WP)
            dyn%uv(2,nz,e) = -1.50_WP*sin(lat)*cos(2.0_WP*lon) + 0.004_WP*real(nz,WP)
            dyn%uv_rhsAB(1,1,nz,e) =  1.0e6_WP*sin(lon)*cos(lat)        + 1.0e4_WP*real(nz,WP)
            dyn%uv_rhsAB(1,2,nz,e) = -1.0e6_WP*cos(lon)*sin(2.0_WP*lat) - 1.0e4_WP*real(nz,WP)
        end do
        ! M2.5 prescribed inputs (elements): Av strictly POSITIVE vertical viscosity
        ! ([1e-3,1.4e-2] m^2/s, realistic) over all nz=1..nl (the TDMA reads Av(nzmax)),
        ! stress_surf sign-varying wind stress (~0.1 N/m^2 -> exercises the surface BC for
        ! both signs). MUST match the FESOM2 oracle.
        do nz = 1, nl
            Av(nz,e) = 5.0e-3_WP + 4.0e-3_WP*cos(lat)*cos(lon) + 1.0e-4_WP*real(nz,WP)
        end do
        stress_surf(1,e) =  0.10_WP*cos(lat)*sin(lon)
        stress_surf(2,e) = -0.08_WP*sin(lat)*cos(2.0_WP*lon)
    end do

    allocate(uv_in(2,nl-1,mesh%elem2D), uv_rhsAB_prev(2,nl-1,mesh%elem2D))
    allocate(uv_rhs_eul(2,nl-1,mesh%elem2D), uv_rhs_ab2(2,nl-1,mesh%elem2D))
    allocate(uv_rhsAB_cor(2,nl-1,mesh%elem2D), uvnode_rhs_dump(2,nl-1,mesh%nod2D))
    allocate(uv_rhs_visc(2,nl-1,mesh%elem2D))
    allocate(u_c_dump(nl-1,mesh%elem2D), v_c_dump(nl-1,mesh%elem2D))
    uv_in         = dyn%uv
    uv_rhsAB_prev = dyn%uv_rhsAB(1,:,:,:)

    dyn%uv_rhs = 0.0_WP
    call compute_vel_rhs(dyn, mesh, dt_velrhs, lfirst=.true.)    ! Euler start (ff=1.0)
    uv_rhs_eul      = dyn%uv_rhs
    uv_rhsAB_cor    = dyn%uv_rhsAB(1,:,:,:)         ! Coriolis + momentum advection (M2.4)
    uvnode_rhs_dump = dyn%work%uvnode_rhs          ! momadv nodal scratch (post-normalize)
    call compute_vel_rhs(dyn, mesh, dt_velrhs, lfirst=.false.)   ! AB2 step (ff=ab2)
    uv_rhs_ab2   = dyn%uv_rhs

    ! M2.4 biharmonic viscosity (opt_visc=7): a SEPARATE operator run AFTER
    ! compute_vel_rhs (FESOM2 oce_ale.F90:3822). It ADDS the biharmonic increment into
    ! dyn%uv_rhs (reads dyn%uv only; the incoming uv_rhs=uv_rhs_ab2 is just accumulated
    ! into). The oracle drives the REAL FESOM2 visc_filt_bidiff the same way.
    call viscosity_filter(7, dyn, mesh, dt_velrhs)
    uv_rhs_visc = dyn%uv_rhs                        ! post-viscosity UV_rhs (gate target)
    u_c_dump    = dyn%work%u_c                       ! first-stage Laplacian (intermediate)
    v_c_dump    = dyn%work%v_c

    ! gate-strength diagnostic (NOT dumped): over interior edges x levels, which term
    ! WINS the viscosity coefficient max(gamma0, max(gamma1*|du|, gamma2*|du|^2))? A
    ! non-zero share for each of gamma0/gamma1/gamma2 confirms all three flow-aware
    ! branches of visc_filt_bidiff are genuinely selected (L11/L17), not just computed.
    dumax = 0.0_WP; ndu = 0; ng0 = 0; ng1 = 0; ng2 = 0
    do ed = 1, mesh%edge2D
        if (ed > mesh%edge2D_in) cycle
        el    = mesh%edge_tri(:, ed)
        nzmin = maxval(mesh%ulevels(el))
        nzmax = minval(mesh%nlevels(el))
        do nz = nzmin, nzmax-1
            du = sqrt((dyn%uv(1,nz,el(1))-dyn%uv(1,nz,el(2)))**2 &
                    + (dyn%uv(2,nz,el(1))-dyn%uv(2,nz,el(2)))**2)
            dumax = max(dumax, du); ndu = ndu + 1
            if (0.003_WP >= max(0.1_WP*du, 0.285_WP*du*du)) then
                ng0 = ng0 + 1
            else if (0.1_WP*du >= 0.285_WP*du*du) then
                ng1 = ng1 + 1
            else
                ng2 = ng2 + 1
            end if
        end do
    end do
    write(*,'(a,es10.3,a,f5.1,a,f5.1,a,f5.1,a)') &
        'fesom_pressuredump: visc strength: max|du|=', dumax, &
        ' m/s ; gamma0/gamma1/gamma2 selected on ', &
        100.0_WP*real(ng0,WP)/real(max(ndu,1),WP), '/', &
        100.0_WP*real(ng1,WP)/real(max(ndu,1),WP), '/', &
        100.0_WP*real(ng2,WP)/real(max(ndu,1),WP), '% of edge-levels'

    ! ============== M2.5 implicit vertical viscosity (TDMA / Thomas solve) ==============
    ! A SEPARATE operator run AFTER viscosity_filter in the timestep (FESOM2
    ! oce_ale.F90:3874-3876). It CONSUMES the post-viscosity uv_rhs_visc as the explicit
    ! rhs, solves the per-column tridiagonal (implicit vertical viscosity Av + vertical
    ! advection w_i + wind-stress top BC + bottom drag), and OVERWRITES dyn%uv_rhs with the
    ! solution. The oracle drives the REAL FESOM2 impl_vert_visc_ale on the same inputs.
    allocate(uv_rhs_ivv(2,nl-1,mesh%elem2D))
    call impl_vert_visc_ale(dyn, mesh, dt_velrhs, Av, stress_surf)
    uv_rhs_ivv = dyn%uv_rhs                           ! post-TDMA UV_rhs (gate target)

    ! gate-strength diagnostic (NOT dumped): the TDMA must non-trivially change the rhs,
    ! and w_i must take BOTH signs so both upwind branches of the advection update fire.
    nwp = count(dyn%w_i >  0.0_WP)
    nwm = count(dyn%w_i <  0.0_WP)
    write(*,'(a,es10.3,a,es10.3,a,i0,a,i0)') &
        'fesom_pressuredump: ivertvisc: max|uv_rhs_ivv|=', maxval(abs(uv_rhs_ivv)), &
        ' ; max|d(uv_rhs)|=', maxval(abs(uv_rhs_ivv - uv_rhs_visc)), &
        ' ; w_i>0/<0 = ', nwp, '/', nwm

    ! ============== M2.6 SSH: stiffness matrix + ssh_rhs + CG solve ===================
    ! The free-surface implicit solve, run AFTER impl_vert_visc_ale and BEFORE the
    ! velocity/SSH update in the timestep (FESOM2 oce_ale.F90:3920-3930). For which_ale=
    ! 'linfs' the stiffness matrix is built ONCE (init_stiff_mat_ale; update_stiff_mat_ale
    ! is skipped, oce_ale.F90:3921), then compute_ssh_rhs_ale assembles the rhs from the
    ! prescribed UV + post-TDMA UV_rhs (dyn%uv_rhs = uv_rhs_ivv) and solve_ssh_ale CG-solves
    ! for d_eta. The FESOM2 oracle drives the REAL routines on the same inputs.
    alpha = 1.0_WP; theta = 1.0_WP      ! pi/reduced-M2: full implicitness (FESOM2 default)
    allocate(dyn%ssh_rhs(mesh%nod2D), dyn%ssh_rhs_old(mesh%nod2D), dyn%d_eta(mesh%nod2D))
    dyn%ssh_rhs_old = 0.0_WP            ! (1-alpha)*ssh_rhs_old = 0 since alpha=1
    dyn%d_eta       = 0.0_WP            ! CG initial guess x0 = 0

    call init_stiff_mat_ale(mesh, dt_ssh)              ! build ssh_stiff (sparsity + values)
    call compute_ssh_rhs_ale(dyn, mesh)               ! assemble dyn%ssh_rhs
    call solve_ssh_ale(dyn, mesh, n_iter=n_ssh_iter)  ! CG -> dyn%d_eta

    ! gate-strength diagnostics (NOT dumped):
    !  (1) the stiffness diagonal + a full matvec A*eta_n (DUMPED, gated) localise the
    !      matrix assembly; computed here only for the prints below,
    !  (2) sum(ssh_rhs) telescopes to ~1e-13 (a proper edge-divergence; plan M2.6),
    !  (3) the CG converged in n_ssh_iter iters and ||A*d_eta - ssh_rhs||_inf shows d_eta
    !      actually solves the system (non-vacuity, the L11 weak-gate guard).
    allocate(ssh_diag(mesh%nod2D), ssh_Aeta(mesh%nod2D))
    do row = 1, mesh%nod2D
        ni  = mesh%ssh_stiff%rowptr_loc(row)
        ne2 = mesh%ssh_stiff%rowptr_loc(row+1) - 1
        ssh_diag(row) = mesh%ssh_stiff%values(ni)
        ssh_Aeta(row) = sum(mesh%ssh_stiff%values(ni:ne2) * dyn%eta_n(mesh%ssh_stiff%colind_loc(ni:ne2)))
    end do
    sum_ssh_rhs = sum(dyn%ssh_rhs(1:mesh%nod2D))
    resid_inf   = 0.0_WP
    do row = 1, mesh%nod2D
        ni  = mesh%ssh_stiff%rowptr_loc(row)
        ne2 = mesh%ssh_stiff%rowptr_loc(row+1) - 1
        resid_inf = max(resid_inf, abs( &
            sum(mesh%ssh_stiff%values(ni:ne2) * dyn%d_eta(mesh%ssh_stiff%colind_loc(ni:ne2))) - dyn%ssh_rhs(row)))
    end do
    write(*,'(a,i0,a,i0,a,es10.3,a,es10.3,a,es10.3)') &
        'fesom_pressuredump: ssh: nza=', mesh%ssh_stiff%nza, ' ; CG iters=', n_ssh_iter, &
        ' ; max|d_eta|=', maxval(abs(dyn%d_eta)), ' ; sum(ssh_rhs)=', sum_ssh_rhs, &
        ' ; ||A d_eta - rhs||inf=', resid_inf

    ! ============== M2.7 ALE (linfs) velocity / SSH / thickness-W update ==============
    ! The post-CG tail of the timestep (FESOM2 oce_ale.F90:3946-4081), run AFTER the SSH
    ! solve. Build the ALE thickness state the update reads/writes: mesh%hbar (prescribed,
    ! the previous-step elevation), hbar_old/dhe (outputs), hnode_new=hnode (linfs never
    ! evolves it), and dyn%w/cfl_z. Force the pi wsplit config (use_wsplit=.true.,
    ! wsplit_maxcfl=1.0). alpha=theta=1.0 already set above. The chain is, in order:
    !   update_vel (UV += UV_rhs + SSH-grad) -> compute_hbar_ale (hbar/dhe/ssh_rhs_old)
    !   -> update_eta_n (eta_n=hbar) -> vert_vel_ale (w + cfl_z + Wvel split). dyn%uv_rhs
    ! is the post-TDMA uv_rhs_ivv; dyn%uv is the prescribed UV (update_vel overwrites it).
    allocate(mesh%hbar(mesh%nod2D), mesh%hbar_old(mesh%nod2D), mesh%dhe(mesh%elem2D))
    allocate(mesh%hnode_new(nl-1, mesh%nod2D))
    allocate(dyn%w(nl, mesh%nod2D), dyn%cfl_z(nl, mesh%nod2D))
    mesh%hbar_old  = 0.0_MP
    mesh%dhe       = 0.0_MP
    mesh%hnode_new = mesh%hnode            ! linfs: hnode_new == hnode (never evolves)
    dyn%w          = 0.0_WP
    dyn%cfl_z      = 0.0_WP
    dyn%use_wsplit    = .true.             ! pi production value (namelist.dyn)
    dyn%wsplit_maxcfl = 1.0_WP

    ! analytic previous-step elevation hbar (sign-varying, ~0.6 m). MUST match the oracle.
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        mesh%hbar(n) = 0.4_WP*sin(lon)*cos(lat) - 0.2_WP*cos(2.0_WP*lat)
    end do

    ! save the prescribed inputs the M2.7 kernels overwrite in place
    allocate(hbar_in(mesh%nod2D), eta_n_in(mesh%nod2D))
    allocate(w_e_in(nl,mesh%nod2D), w_i_in(nl,mesh%nod2D), uv_upd(2,nl-1,mesh%elem2D))
    hbar_in  = mesh%hbar
    eta_n_in = dyn%eta_n
    w_e_in   = dyn%w_e
    w_i_in   = dyn%w_i

    call update_vel(dyn, mesh, dt_velrhs)           ! UV += UV_rhs + [-g*theta*dt*grad(d_eta)]
    uv_upd = dyn%uv
    call compute_hbar_ale(dyn, mesh, dt_velrhs)     ! ssh_rhs_old, hbar_old, hbar, dhe
    call update_eta_n(dyn, mesh)                    ! eta_n = alpha*hbar + (1-alpha)*hbar_old
    call vert_vel_ale(dyn, mesh, dt_velrhs)         ! w + cfl_z + w_e/w_i split

    ! gate-strength diagnostic (NOT dumped): the Wvel split fires only where CFL_z >
    ! wsplit_maxcfl. A non-zero share confirms the split formula (dd, Wvel_e/Wvel_i) is
    ! genuinely exercised (not just the trivial Wvel_e=Wvel branch), the L11/L17 guard.
    n_cflsplit = 0
    do n = 1, mesh%nod2D
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n)
            if (dyn%cfl_z(nz,n) > dyn%wsplit_maxcfl) n_cflsplit = n_cflsplit + 1
        end do
    end do
    write(*,'(a,es10.3,a,es10.3,a,es10.3,a,i0)') &
        'fesom_pressuredump: ale: max|uv_upd|=', maxval(abs(uv_upd)), &
        ' ; max|hbar|=', maxval(abs(mesh%hbar)), ' ; max|w|=', maxval(abs(dyn%w)), &
        ' ; CFL_z>maxcfl on ', n_cflsplit

    ! --- dump (same FADVHDMP format / names as the FESOM2 oracle) ---
    call advhor_dump_open(u, trim(out_path), mesh%nod2D, mesh%elem2D, mesh%edge2D, nl)
    call wr_r2(u, 'temp',           real(temp,        MP))
    call wr_r2(u, 'salt',           real(salt,        MP))
    call wr_r2(u, 'density_ref',    real(density_ref, MP))
    call wr_r2(u, 'zbar_3d_n',      mesh%zbar_3d_n(1:nl,   1:mesh%nod2D))
    call wr_r2(u, 'Z_3d_n',         mesh%Z_3d_n(1:nl-1,    1:mesh%nod2D))
    call wr_r2(u, 'hnode',          mesh%hnode(1:nl-1,     1:mesh%nod2D))
    call wr_r2(u, 'density_m_rho0', real(density,            MP))
    call wr_r2(u, 'hpressure',      real(hpressure(1:nl-1,:), MP))
    call wr_r2(u, 'bvfreq_raw',     real(bvfreq_raw,         MP))
    call wr_r2(u, 'bvfreq',         real(bvfreq,             MP))
    call wr_r2(u, 'pgf_x',          real(pgf_x(1:nl-1, :),   MP))
    call wr_r2(u, 'pgf_y',          real(pgf_y(1:nl-1, :),   MP))
    ! M2.3 vel_rhs + M2.4 momadv: gated coriolis + prescribed inputs (incl. w_e) +
    ! the momadv nodal intermediate (uvnode_rhs) + the full-assembly outputs.
    call wr_r1(u, 'coriolis',       real(mesh%coriolis(1:mesh%elem2D), MP))
    call wr_r1(u, 'eta_n',          real(eta_n_in,         MP))   ! prescribed (pre-M2.7)
    call wr_r3(u, 'uv_in',          real(uv_in,            MP))
    call wr_r3(u, 'uv_rhsAB_prev',  real(uv_rhsAB_prev,    MP))
    call wr_r2(u, 'w_e',            real(w_e_in(1:nl, :),  MP))   ! prescribed (pre-split)
    call wr_r3(u, 'uvnode_rhs',     real(uvnode_rhs_dump,  MP))
    call wr_r3(u, 'uv_rhsAB_cor',   real(uv_rhsAB_cor,     MP))
    call wr_r3(u, 'uv_rhs_eul',     real(uv_rhs_eul,       MP))
    call wr_r3(u, 'uv_rhs_ab2',     real(uv_rhs_ab2,       MP))
    ! M2.4 biharmonic viscosity (opt_visc=7): the first-stage Laplacian intermediate
    ! u_c/v_c + the post-viscosity UV_rhs (gate target).
    call wr_r2(u, 'visc_u_c',       real(u_c_dump,         MP))
    call wr_r2(u, 'visc_v_c',       real(v_c_dump,         MP))
    call wr_r3(u, 'uv_rhs_visc',    real(uv_rhs_visc,      MP))
    ! M2.5 implicit vertical viscosity (TDMA): prescribed inputs (Av/stress_surf/w_i) +
    ! the post-solve UV_rhs (gate target).
    call wr_r2(u, 'Av',             real(Av(1:nl, :),      MP))
    call wr_r2(u, 'stress_surf',    real(stress_surf,      MP))
    call wr_r2(u, 'w_i',            real(w_i_in(1:nl, :),  MP))   ! prescribed (pre-split)
    call wr_r3(u, 'uv_rhs_ivv',     real(uv_rhs_ivv,       MP))
    ! M2.6 SSH: the stiffness diagonal + matvec A*eta_n (localise the matrix assembly),
    ! the assembled ssh_rhs, and the CG solution d_eta (the gate target).
    call wr_r1(u, 'ssh_stiff_diag', real(ssh_diag,    MP))
    call wr_r1(u, 'ssh_Aeta',       real(ssh_Aeta,    MP))
    call wr_r1(u, 'ssh_rhs',        real(dyn%ssh_rhs, MP))
    call wr_r1(u, 'd_eta',          real(dyn%d_eta,   MP))
    ! M2.7 ALE velocity/SSH/thickness-W update: the prescribed hbar input + the updated
    ! UV (update_vel) + the divergence ssh_rhs_old + new hbar/dhe (compute_hbar_ale) + the
    ! blended eta_n + the vertical velocity w (vert_vel_ale) + hnode_new (=hnode, linfs) +
    ! cfl_z (compute_CFLz) + the explicit/implicit Wvel split (compute_Wvel_split).
    call wr_r1(u, 'hbar_in',        real(hbar_in,          MP))
    call wr_r3(u, 'uv_upd',         real(uv_upd,           MP))
    call wr_r1(u, 'ssh_rhs_old',    real(dyn%ssh_rhs_old,  MP))
    call wr_r1(u, 'hbar',           real(mesh%hbar,        MP))
    call wr_r1(u, 'dhe',            real(mesh%dhe,         MP))
    call wr_r1(u, 'eta_n_upd',      real(dyn%eta_n,        MP))
    call wr_r2(u, 'w',              real(dyn%w(1:nl, :),   MP))
    call wr_r2(u, 'hnode_new',      real(mesh%hnode_new(1:nl-1, :), MP))
    call wr_r2(u, 'cfl_z',          real(dyn%cfl_z(1:nl, :), MP))
    call wr_r2(u, 'w_split_e',      real(dyn%w_e(1:nl, :), MP))
    call wr_r2(u, 'w_split_i',      real(dyn%w_i(1:nl, :), MP))
    call advhor_dump_close(u)
    write(*,'(a)') 'fesom_pressuredump: wrote '//trim(out_path)

    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_pressuredump
