program fesom_pressuredump
    ! M2.1 pressure/EOS operator byte-gate driver. Loads pi 1-rank with the SAME
    ! rotation as the FESOM2 pi run (alpha/beta/gamma=50/15/-90, cyclic 360),
    ! computes geometry, builds the ALE depths Z_3d_n/zbar_3d_n + node thickness
    ! hnode (linfs full cells; M1.2/M1.3 proved them), PRESCRIBES an analytic T/S
    ! field (identical formula to the FESOM2 oracle src/fesom_pressure_dump.F90,
    ! computed from the byte-identical coordinates) + density_ref=density_0, runs
    ! oce_pressure_bv (once with horizontal N^2 smoothing OFF, once ON), and dumps
    ! density_m_rho0 / hpressure / bvfreq (raw + smoothed) in the oracle's FADVHDMP
    ! format. tools/pressure_diff.py then compares for max|delta|=0.
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
    use oce_pressure_bv,  only: pressure_bv
    use oce_pgf,          only: pressure_force_4_linfs_fullcell
    use mod_advhor_dump,  only: advhor_dump_open, advhor_dump_close, wr_r2
    implicit none

    character(len=512) :: mesh_dir, out_path
    type(t_partit)      :: partit
    type(t_mesh)        :: mesh
    integer :: nsw, n, nz, nl, u, nzmin, nzmax
    real(kind=WP) :: lon, lat
    real(kind=MP) :: zbar_srf, zbar_bot
    real(kind=WP), allocatable :: temp(:,:), salt(:,:), density_ref(:,:)
    real(kind=WP), allocatable :: density(:,:), hpressure(:,:), bvfreq(:,:)
    real(kind=WP), allocatable :: bvfreq_raw(:,:), pgf_x(:,:), pgf_y(:,:)

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
    call advhor_dump_close(u)
    write(*,'(a)') 'fesom_pressuredump: wrote '//trim(out_path)

    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_pressuredump
