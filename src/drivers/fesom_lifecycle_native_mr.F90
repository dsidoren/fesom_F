program fesom_lifecycle_native_mr
    ! M3f-4: the MULTI-RANK native-flux forced lifecycle driver — the multi-rank analog of
    ! the 1-rank fesom_lifecycle_native (and the native-flux analog of fesom_lifecycle_mr).
    ! Runs on dist_<NP> (NP>=1) the REAL free-running coupled sea-ice + ocean lifecycle with
    ! the WHOLE air-sea forcing computed NATIVELY each step (NO prescribed atmosphere, NO
    ! prescribed fluxes):
    !
    !   native atmosphere  8 NCAR fields (mod_forcing_read) + NCAR bulk Ch/Ce + stress_atmoce
    !                      + wind-on-ice stress (mod_forcing_bulk) + runoff/Ssurf monthly
    !                      climatology (mod_forcing_other), ALL over owned+halo (nNodL)
    !   ocean2ice          ocean surface -> srfoce_u/v/temp/salt/ssh                  (M3b)
    !   ice_timestep       EVP -> FCT -> cut_off -> thermodynamics                     (M3f-1)
    !   oce_fluxes_mom     ice+atm momentum stress -> stress_surf                      (M3e)
    !   oce_fluxes         heat/freshwater/salt budget -> the tracer/SSH BCs           (M3e)
    !   step_oce           the whole ocean step on the NATIVE fluxes (through partit)  (M2.12)
    !
    ! Each kernel runs the M2.12 optional-partit path: at npes==1 the proven 1-rank code
    ! VERBATIM, at npes>1 owned/halo loop bounds + the FESOM2 halo exchanges. The whole atm
    ! forcing + ice + budget chain is partition-independent per-node (computed over owned+halo
    ! from halo-valid inputs, no extra exchange); ocean2ice exchanges u_w/v_w, the EVP exchanges
    ! uice/vice, the FCT exchanges its solve intermediates, and oce_fluxes' integrate_nod /
    ! mesh%ocean_area do the cross-rank allreduce_sum. step_oce is the byte-proven M2.12c MR step.
    !
    ! The gate (tools/run_lifecycle_fullynative_gate_multirank.sh) compares the per-substep
    ! gid-keyed NODE dumps (mod_dump) PER RANK vs the REAL FESOM2 forced lifecycle on the SAME
    ! dist_<NP> partition (the L8 same-partition rule). NOTE: the oracle flux dump
    ! (fesom_flux_dump) is skipped at npes/=1, so there is no per-step flux self-check at
    ! multi-rank — the node-substep gate is the validation (the 1-rank fully-native gate's
    ! self-check already proved the native fluxes byte-exact).
    !
    !   FESOM3_MESH_DIR     mesh dir (must contain dist_<NP>/)  default: CORE2
    !   FESOM3_IC_FILE      IC netcdf                           default: pool phc3.0_winter.nc
    !   FESOM3_FORCING_DIR  native NCAR forcing dir             (REQUIRED — fully native)
    !   FESOM3_RUNOFF_FILE  runoff climatology                  default: pool CORE2_runoff.nc
    !   FESOM3_SSS_FILE     SSS restoring climatology           default: pool PHC2_salx.nc
    !   FESOM3_WHICHEVP     0=EVP / 1=mEVP                       (default 0)
    !   FESOM_DUMP_FILE     per-rank node dump prefix (mod_dump; <prefix>.<mype5>)
    !   FESOM_DUMP_MAXSTEPS dump step cap
    !   FESOM3_NSTEPS       number of steps                     (default: 3)
    use mpi
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use mod_precision,      only: WP, MP
    use mod_constants,      only: density_0
    use mod_param_phys,     only: N2smth_h, alpha, theta
    use mod_param_phys,     only: mix_coeff_PP, A_ver, K_ver, Kv0_const, mix_scheme_nmb
    use mod_param_phys,     only: use_instabmix, instabmix_kv, use_momix, use_windmix
    ! M4f: GM bolus (FESOM3_FER_GM) + Redi isopycnal diffusion (FESOM3_REDI) at MULTI-RANK —
    ! the work_core GM+Redi config (the M4 routines already thread the optional partit + exchanges).
    use mod_param_phys,     only: Fer_GM, Redi, K_GM_max, K_GM_min, K_GM_bvref, &
                                  K_GM_rampmax, K_GM_rampmin, K_GM_resscalorder, K_GM_cm, &
                                  K_GM_cmin, K_GM_Ktaper, scaling_Ferreira, scaling_Rossby, &
                                  scaling_resolution, scaling_FESOM14, scaling_GMzexp, &
                                  scaling_GINsea, GMzexp_zref, GMzexp_smin
    use mod_param_phys,     only: Redi_Kmax, Redi_Kmin, Redi_Ktaper, K_hor, &
                                  scaling_ODM95, ODM95_Scr, ODM95_Sd, scaling_LDD97
    use mod_mesh,           only: t_mesh
    use mod_partit,         only: t_partit
    use mod_partitioning,   only: par_init, par_ex, set_partition
    use mod_mesh_read,      only: read_mesh
    use mod_mesh_areas,     only: compute_geometry
    use mod_dyn,            only: t_dyn
    use mod_tracer,         only: t_tracer
    use mod_ice,            only: t_ice
    use mod_ice_setup,      only: ice_setup
    use mod_ice_dyn,        only: ocean2ice
    use mod_ice_step,       only: ice_timestep
    use mod_ice_thermo,     only: t_atmflux
    use mod_ice_oce_coupling, only: oce_fluxes_mom, oce_fluxes
    use mod_forcing_read    ! native CORE2 atmospheric forcing read
    use mod_forcing_bulk,   only: forcing_bulk_ncar, forcing_wind_stress, forcing_ice_stress
    use mod_forcing_other,  only: read_other_NetCDF   ! runoff + SSS climatology read
    use oce_initial_state,  only: t_ic3d_config, do_ic3d
    use oce_muscl_adv,      only: muscl_adv_init
    use oce_ssh_rhs,        only: init_stiff_mat_ale
    use mod_step_oce,       only: step_oce
    use mod_dump,           only: dump_init, dump_finalize
    implicit none

    real(kind=WP), parameter :: dt = 86400.0_WP / real(48, WP)   ! CORE2 dt = 1800 s

    character(len=512) :: mesh_dir, ic_file, env, whichevp_str
    type(t_partit)     :: partit
    type(t_mesh)       :: mesh
    type(t_dyn)        :: dyn
    type(t_tracer)     :: tracers
    type(t_ice)        :: ice
    type(t_atmflux)    :: atm
    type(t_ic3d_config):: ic
    integer :: n, nz, nl, nzmin, nzmax, e, tr_num, nsteps, ios, env_len, nsw, whichevp
    integer :: nNodO, nNodL, nEdgeO, nElemO, nElemF
    real(kind=MP) :: zbar_srf, zbar_bot
    real(kind=WP), allocatable :: Ki(:,:), real_salt_flux(:), stress_surf(:,:)
    real(kind=WP) :: is_nonlinfs
    logical :: use_fer_gm, use_redi   ! M4f: GM bolus / Redi isopycnal diffusion toggles
    ! native CORE2 forcing read (the whole atmosphere, over owned+halo).
    character(len=512)  :: forcing_dir, runoff_file, sss_file
    type(t_atm_forcing) :: frc
    integer             :: fld, fyear
    real(kind=WP)       :: rdate_cold
    real(kind=WP), allocatable :: nuw(:), nvw(:), nta(:), nsh(:), nswr(:), nlw(:), npr(:), nps(:)
    real(kind=WP), allocatable :: ncd(:), nch(:), nce(:), nsx(:), nsy(:), nix(:), niy(:)
    real(kind=WP), allocatable :: nro(:), nss(:)
    integer, allocatable :: idlist(:)

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2'
    call get_environment_variable('FESOM3_IC_FILE', ic_file)
    if (len_trim(ic_file) == 0) &
        ic_file = '/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc'
    nsteps = 3
    call get_environment_variable('FESOM3_NSTEPS', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) read(env, *, iostat=ios) nsteps
    whichevp = 0
    call get_environment_variable('FESOM3_WHICHEVP', whichevp_str)
    if (len_trim(whichevp_str) > 0) read(whichevp_str, *, iostat=ios) whichevp
    call get_environment_variable('FESOM3_FER_GM', env, length=env_len, status=ios)
    use_fer_gm = (ios == 0 .and. env_len > 0)
    call get_environment_variable('FESOM3_REDI', env, length=env_len, status=ios)
    use_redi = (ios == 0 .and. env_len > 0)

    !===========================================================================
    ! model_init: MR mesh remap + geometry (set_partition -> read_mesh dispatches to
    ! read_mesh_local at npes>1; npes==1 reads the global mesh).
    call par_init(partit)
    call set_partition(partit, trim(mesh_dir))
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    nl = mesh%nl
    if (partit%npes == 1) then
        nNodO = mesh%nod2D;   nNodL  = mesh%nod2D
        nEdgeO = mesh%edge2D; nElemO = mesh%elem2D; nElemF = mesh%elem2D
    else
        nNodO  = partit%myDim_nod2D
        nNodL  = partit%myDim_nod2D + partit%eDim_nod2D
        nEdgeO = partit%myDim_edge2D
        nElemO = partit%myDim_elem2D
        nElemF = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
    end if
    if (partit%mype == 0) &
        write(*,'(a,i0,a,i0,a,i0,a,i0)') 'fesom_lifecycle_native_mr: nod2D=', mesh%nod2D, &
            ' elem2D=', mesh%elem2D, ' nl=', nl, ' CW swaps=', nsw

    !===========================================================================
    ! ALE depth/thickness state (linfs full cells; local sizes) — as fesom_lifecycle_mr.
    allocate(mesh%hnode(nl-1, nNodL)); mesh%hnode = 0.0_MP
    do n = 1, nNodL
        if (mesh%nlevels_nod2D(n) <= 0) cycle
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n)-1
            mesh%hnode(nz, n) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
    end do
    allocate(mesh%helem(nl-1, nElemF)); mesh%helem = 0.0_MP
    allocate(mesh%zbar_e_bot(nElemF)); mesh%zbar_e_bot = 0.0_MP
    do e = 1, nElemF
        if (mesh%nlevels(e) <= 0) cycle
        do nz = mesh%ulevels(e), mesh%nlevels(e)-1
            mesh%helem(nz, e) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
        mesh%zbar_e_bot(e) = mesh%zbar(mesh%nlevels(e))
    end do
    allocate(mesh%zbar_3d_n(nl, nNodL), mesh%Z_3d_n(nl-1, nNodL))
    mesh%zbar_3d_n = 0.0_MP; mesh%Z_3d_n = 0.0_MP
    do n = 1, nNodL
        if (mesh%nlevels_nod2D(n) <= 0) cycle
        nzmin = mesh%ulevels_nod2D(n); nzmax = mesh%nlevels_nod2D(n)
        zbar_srf = mesh%zbar(nzmin); zbar_bot = mesh%zbar(nzmax)
        mesh%zbar_3d_n(1:nzmin-1, n)       = mesh%zbar(1:nzmin-1)
        mesh%zbar_3d_n(nzmin, n)           = zbar_srf
        mesh%zbar_3d_n(nzmin+1:nzmax-1, n) = mesh%zbar(nzmin+1:nzmax-1)
        mesh%zbar_3d_n(nzmax, n)           = zbar_bot
        mesh%Z_3d_n(1:nzmin-1, n)          = mesh%Z(1:nzmin-1)
        mesh%Z_3d_n(nzmin, n)              = mesh%zbar_3d_n(nzmin,n)   + (mesh%zbar_3d_n(nzmin+1,n)-zbar_srf)/2
        mesh%Z_3d_n(nzmin+1:nzmax-2, n)    = mesh%Z(nzmin+1:nzmax-2)
        mesh%Z_3d_n(nzmax-1, n)            = mesh%zbar_3d_n(nzmax-1,n) + (zbar_bot-mesh%zbar_3d_n(nzmax-1,n))/2
    end do
    allocate(mesh%hbar(nNodL), mesh%hbar_old(nNodL), mesh%dhe(nElemF))
    allocate(mesh%hnode_new(nl-1, nNodL))
    mesh%hbar = 0.0_MP; mesh%hbar_old = 0.0_MP; mesh%dhe = 0.0_MP
    mesh%hnode_new = mesh%hnode

    !===========================================================================
    ! dynamics state — COLD START (local sizes).
    allocate(dyn%uv(2, nl-1, nElemF), dyn%uv_rhs(2, nl-1, nElemF))
    allocate(dyn%uv_rhsAB(1, 2, nl-1, nElemF))
    allocate(dyn%uvnode(2, nl-1, nNodL))
    allocate(dyn%eta_n(nNodL), dyn%d_eta(nNodL))
    allocate(dyn%ssh_rhs(nNodL), dyn%ssh_rhs_old(nNodL))
    allocate(dyn%w(nl, nNodL), dyn%w_e(nl, nNodL), dyn%w_i(nl, nNodL))
    allocate(dyn%cfl_z(nl, nNodL))
    allocate(dyn%work%density_ref(nl-1, nNodL), dyn%work%density_m_rho0(nl-1, nNodL))
    allocate(dyn%work%hpressure(nl, nNodL), dyn%work%bvfreq(nl, nNodL))
    allocate(dyn%work%pgf_x(nl-1, nElemF), dyn%work%pgf_y(nl-1, nElemF))
    allocate(dyn%work%u_c(nl-1, nElemF), dyn%work%v_c(nl-1, nElemF))
    allocate(dyn%work%uvnode_rhs(2, nl-1, nNodL))
    allocate(dyn%work%Kv(nl, nNodL), dyn%work%Av(nl, nElemF))
    dyn%uv = 0.0_WP; dyn%uv_rhs = 0.0_WP; dyn%uv_rhsAB = 0.0_WP; dyn%uvnode = 0.0_WP
    dyn%eta_n = 0.0_WP; dyn%d_eta = 0.0_WP; dyn%ssh_rhs = 0.0_WP; dyn%ssh_rhs_old = 0.0_WP
    dyn%w = 0.0_WP; dyn%w_e = 0.0_WP; dyn%w_i = 0.0_WP; dyn%cfl_z = 0.0_WP
    dyn%work%density_m_rho0 = 0.0_WP; dyn%work%hpressure = 0.0_WP; dyn%work%bvfreq = 0.0_WP
    dyn%work%pgf_x = 0.0_WP; dyn%work%pgf_y = 0.0_WP
    dyn%work%u_c = 0.0_WP; dyn%work%v_c = 0.0_WP; dyn%work%uvnode_rhs = 0.0_WP
    dyn%work%Kv = 0.0_WP; dyn%work%Av = 0.0_WP
    dyn%work%density_ref = density_0
    dyn%AB_order      = 2
    dyn%momadv_opt    = 2
    dyn%opt_visc      = 7
    dyn%visc_gamma0   = 0.003_WP
    dyn%visc_gamma1   = 0.1_WP
    dyn%visc_gamma2   = 0.285_WP
    dyn%visc_gamma0_h = 0.0_WP
    dyn%visc_gamma1_h = 0.0_WP
    dyn%use_wsplit    = .false.
    dyn%wsplit_maxcfl = 1.0_WP

    !===========================================================================
    ! 2-tracer ocean state (data(1)=T, data(2)=S; local sizes) + do_ic3d phc3.0 IC.
    tracers%num_tracers = 2
    allocate(tracers%data(2))
    allocate(tracers%data(1)%values(nl-1, nNodL), tracers%data(2)%values(nl-1, nNodL))
    tracers%data(1)%values = 0.0_WP;  tracers%data(1)%ID = 1
    tracers%data(2)%values = 0.0_WP;  tracers%data(2)%ID = 2
    ic%n_ic3d      = 2
    ic%idlist(1:2) = [2, 1]
    ic%filelist(1) = trim(ic_file)
    ic%filelist(2) = trim(ic_file)
    ic%varlist(1)  = 'salt'
    ic%varlist(2)  = 'temp'
    ic%t_insitu    = .true.
    ic%ic_cyclic   = .true.
    ic%dummy       = 1.e10_WP
    call do_ic3d(tracers, ic, mesh, partit)
    if (partit%mype == 0) &
        write(*,'(a,2es12.4,a,2es12.4)') 'fesom_lifecycle_native_mr: IC(rank0 owned) T=', &
            minval(tracers%data(1)%values(:,1:nNodO)), maxval(tracers%data(1)%values(:,1:nNodO)), &
            '  S=', minval(tracers%data(2)%values(:,1:nNodO)), maxval(tracers%data(2)%values(:,1:nNodO))

    !===========================================================================
    ! tracer advection machinery (cold start = values; local sizes; adv_flux_hor owned edges).
    do tr_num = 1, 2
        allocate(tracers%data(tr_num)%valuesAB(nl-1, nNodL))
        allocate(tracers%data(tr_num)%valuesold(2, nl-1, nNodL))
        tracers%data(tr_num)%AB_order   = 2
        tracers%data(tr_num)%tra_adv_hor = 'MFCT'
        tracers%data(tr_num)%tra_adv_ver = 'QR4C'
        tracers%data(tr_num)%tra_adv_lim = 'FCT'
        tracers%data(tr_num)%tra_adv_ph  = 0.0_WP
        tracers%data(tr_num)%tra_adv_pv  = 1.0_WP
        tracers%data(tr_num)%i_vert_diff = .true.
        tracers%data(tr_num)%valuesAB        = 0.0_WP
        tracers%data(tr_num)%valuesold(1,:,:) = tracers%data(tr_num)%values
        tracers%data(tr_num)%valuesold(2,:,:) = tracers%data(tr_num)%values
    end do
    allocate(tracers%work%fct_LO          (nl-1, nNodL))
    allocate(tracers%work%adv_flux_hor    (nl-1, nEdgeO))
    allocate(tracers%work%adv_flux_ver    (nl,   nNodL))
    allocate(tracers%work%fct_ttf_max     (nl-1, nNodL))
    allocate(tracers%work%fct_ttf_min     (nl-1, nNodL))
    allocate(tracers%work%fct_plus        (nl-1, nNodL))
    allocate(tracers%work%fct_minus       (nl-1, nNodL))
    allocate(tracers%work%del_ttf         (nl-1, nNodL))
    allocate(tracers%work%del_ttf_advhoriz(nl-1, nNodL))
    allocate(tracers%work%del_ttf_advvert (nl-1, nNodL))
    call muscl_adv_init(tracers%work, mesh, partit)

    !===========================================================================
    ! reduced-M2 module config (= the FESOM2 oracle / CORE2 namelist).
    alpha = 1.0_WP; theta = 1.0_WP
    N2smth_h     = .true.
    mix_scheme_nmb = 2           ! reduced-M2 mixing = PP (KPP is M5c/M5d)
    mix_coeff_PP = 0.01_WP
    A_ver        = 1.0e-4_WP
    K_ver        = 1.0e-5_WP
    Kv0_const    = .true.
    use_instabmix = .true.
    instabmix_kv  = 0.1_WP
    use_momix     = .false.
    use_windmix   = .false.

    !===========================================================================
    ! M4f GM/Redi (FESOM3_FER_GM / FESOM3_REDI): enable Gent-McWilliams bolus advection
    ! (+ Redi isopycnal diffusion) at MULTI-RANK. Same work_core config + inits as the 1-rank
    ! fesom_lifecycle_native (M4e), but the GM/Redi work arrays are LOCAL-sized (nNodL/nElemF)
    ! like the rest of the MR state. step_oce threads the optional partit to every M4 routine
    ! (producers + init_Redi_GM/fer_solve_Gamma/fer_gamma2vel + the Redi diff terms), each of
    ! which already does owned-loop bounds + the FESOM2 halo exchanges (sigma_xy/neutral_slope/
    ! slope_tapered/fer_c/fer_K/Ki/fer_gamma exchange_nod, fer_uv exchange_elem, fer_w in
    ! vert_vel_ale). sw_alpha/sw_beta/tr_z are computed over owned+halo (no exchange needed).
    if (use_fer_gm .or. use_redi) then
        Fer_GM = .true.        ! Redi rides the GM coupling (Ki = max(fer_scal*Redi_Kmax, K_GM_min))
        K_GM_max = 1000.0_WP; K_GM_min = 2.0_WP; K_GM_bvref = 1
        K_GM_rampmax = -1.0_WP; K_GM_rampmin = -1.0_WP; K_GM_resscalorder = 2.0_WP
        K_GM_cm = 3.0_WP; K_GM_cmin = 0.1_WP; K_GM_Ktaper = .false.
        scaling_Ferreira = .false.; scaling_Rossby = .false.; scaling_resolution = .true.
        scaling_FESOM14 = .false.; scaling_GMzexp = .true.; scaling_GINsea = .false.
        GMzexp_zref = 500.0_WP; GMzexp_smin = 0.6_WP
        allocate(dyn%fer_uv(2, nl-1, nElemF), dyn%fer_w(nl, nNodL))
        allocate(dyn%work%sw_alpha(nl-1, nNodL), dyn%work%sw_beta(nl-1, nNodL))
        allocate(dyn%work%sigma_xy(2, nl-1, nNodL))
        allocate(dyn%work%fer_K(nl, nNodL), dyn%work%fer_c(nNodL), dyn%work%fer_scal(nNodL))
        allocate(dyn%work%fer_gamma(2, nl, nNodL))
        dyn%fer_uv = 0.0_WP; dyn%fer_w = 0.0_WP
        dyn%work%sw_alpha = 0.0_WP; dyn%work%sw_beta = 0.0_WP; dyn%work%sigma_xy = 0.0_WP
        dyn%work%fer_K = 500.0_WP; dyn%work%fer_c = 1.0_WP; dyn%work%fer_scal = 0.0_WP
        dyn%work%fer_gamma = 0.0_WP
        if (use_redi) then
            Redi = .true.; Redi_Ktaper = .true.; Redi_Kmax = 0.0_WP; Redi_Kmin = 100.0_WP
            scaling_ODM95 = .true.; ODM95_Scr = 0.2e-2_WP; ODM95_Sd = 1.0e-3_WP; scaling_LDD97 = .false.
            K_hor = 0.0_WP
            allocate(dyn%work%Ki(nl-1, nNodL), dyn%work%fer_tapfac(nl-1, nNodL))
            allocate(dyn%work%neutral_slope(3, nl-1, nNodL), dyn%work%slope_tapered(3, nl-1, nNodL))
            allocate(tracers%work%tr_z(nl, nNodL))
            dyn%work%Ki = 0.0_WP; dyn%work%fer_tapfac = 1.0_WP
            dyn%work%neutral_slope = 0.0_WP; dyn%work%slope_tapered = 0.0_WP
            tracers%work%tr_z = 0.0_WP
            if (partit%mype == 0) write(*,'(a)') &
                'fesom_lifecycle_native_mr: Fer_GM + Redi ENABLED (work_core GM+Redi config)'
        else
            Redi = .false.
            if (partit%mype == 0) write(*,'(a)') &
                'fesom_lifecycle_native_mr: Fer_GM ENABLED (work_core GM config; Redi off)'
        end if
    end if

    ! SSH stiffness (built ONCE; dt = CORE2 namelist timestep).
    call init_stiff_mat_ale(mesh, dt, partit)

    !===========================================================================
    ! ice setup (allocate local-sized + ice_mass_matrix_fill + cold-start ice_initial_state)
    ! + the CORE2 &ice_dyn / &ice_therm config overrides (same as fesom_lifecycle_native).
    call ice_setup(ice, tracers, mesh, dt, partit)
    ice%whichEVP      = whichevp
    ice%cd_oce_ice    = 0.0055_WP
    ice%delta_min     = 1.0e-11_WP
    ice%ice_diff      = 0.0_WP
    ice%ice_gamma_fct = 0.5_WP
    ice%thermo%Sice              = 4.0_WP
    ice%thermo%iclasses          = 7
    ice%thermo%new_iclasses      = .false.
    ice%thermo%h_cutoff          = 3.0_WP
    ice%thermo%h0                = 0.5_WP
    ice%thermo%h0_s              = 0.5_WP
    ice%thermo%hmin              = 0.01_WP
    ice%thermo%Armin             = 0.01_WP
    ice%thermo%emiss_ice         = 0.97_WP
    ice%thermo%emiss_wat         = 0.97_WP
    ice%thermo%albsn             = 0.81_WP
    ice%thermo%albsnm            = 0.77_WP
    ice%thermo%albi              = 0.7_WP
    ice%thermo%albim             = 0.68_WP
    ice%thermo%albw              = 0.1_WP
    ice%thermo%open_water_albedo = 0
    ice%thermo%con               = 2.1656_WP
    ice%thermo%consn             = 0.31_WP
    ice%thermo%snowdist          = .true.
    ice%thermo%c_melt            = 0.5_WP

    !===========================================================================
    ! atmflux: allocate (local nNodL) + scalar config (the CORE2 namelist scalars).
    call alloc_atm(atm, nNodL)
    atm%Ch_atm_ice    = 0.00175_WP
    atm%Ce_atm_ice    = 0.00175_WP
    atm%ref_sss       = 34.0_WP
    atm%ref_sss_local = .true.
    atm%use_virt_salt = .true.
    atm%l_snow        = .true.
    atm%surf_relax_S  = 1.929e-06_WP

    !===========================================================================
    ! native CORE2 forcing read setup (REQUIRED — fully native). frc%nnod = nNodL so the
    ! 8 NCAR fields read + bilinear + g2r-rotate fill the OWNED+HALO nodes (each node
    ! interpolates from the full global raw grid — partition-independent, no exchange).
    call get_environment_variable('FESOM3_FORCING_DIR', forcing_dir)
    if (len_trim(forcing_dir) == 0) then
        if (partit%mype == 0) write(*,'(a)') &
            'fesom_lifecycle_native_mr: FULLY NATIVE requires FESOM3_FORCING_DIR'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype); error stop 1
    end if
    fyear = 1948
    frc%nfld = 8;  frc%nnod = nNodL
    frc%iyear = 1948; frc%imm = 1; frc%idd = 1; frc%freq = 1; frc%tmid = 1
    frc%ic_cyclic = .true.; frc%rotated_grid = .true.
    frc%i_xwind = 1; frc%i_ywind = 2
    frc%f(1)%file_base = trim(forcing_dir)//'/u_10.';        frc%f(1)%varname = 'U_10_MOD'
    frc%f(2)%file_base = trim(forcing_dir)//'/v_10.';        frc%f(2)%varname = 'V_10_MOD'
    frc%f(3)%file_base = trim(forcing_dir)//'/q_10.';        frc%f(3)%varname = 'Q_10_MOD'
    frc%f(4)%file_base = trim(forcing_dir)//'/ncar_rad.';    frc%f(4)%varname = 'SWDN_MOD'
    frc%f(5)%file_base = trim(forcing_dir)//'/ncar_rad.';    frc%f(5)%varname = 'LWDN_MOD'
    frc%f(6)%file_base = trim(forcing_dir)//'/t_10.';        frc%f(6)%varname = 'T_10_MOD'
    frc%f(7)%file_base = trim(forcing_dir)//'/ncar_precip.'; frc%f(7)%varname = 'RAIN'
    frc%f(8)%file_base = trim(forcing_dir)//'/ncar_precip.'; frc%f(8)%varname = 'SNOW'
    rdate_cold = real(forcing_julday(fyear,1,1,'noleap'),WP)
    call forcing_alloc(frc)
    do fld = 1, frc%nfld
        call forcing_read_grid(frc, fld, fyear)
        call forcing_build_bilin(frc, fld, mesh, partit)
    end do
    do fld = 1, frc%nfld
        call forcing_getcoeffld(frc, fld, fyear, rdate_cold, mesh, partit)
    end do
    call forcing_rotate_wind(frc, mesh, partit)
    allocate(nuw(nNodL), nvw(nNodL), nta(nNodL), nsh(nNodL), &
             nswr(nNodL), nlw(nNodL), npr(nNodL), nps(nNodL), &
             ncd(nNodL), nch(nNodL), nce(nNodL), &
             nsx(nNodL), nsy(nNodL), nix(nNodL), niy(nNodL))
    ! native runoff + SSS climatology (read ONCE — both constant for a Jan run). PASS partit
    ! so read_other_NetCDF interpolates over owned+halo (myDim+eDim; read_dist_partition sets
    ! myDim — the L42 trap is the 1-rank case where myDim=0, here npes>1 sets it).
    call get_environment_variable('FESOM3_RUNOFF_FILE', runoff_file)
    if (len_trim(runoff_file) == 0) &
        runoff_file = '/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/CORE2_runoff.nc'
    call get_environment_variable('FESOM3_SSS_FILE', sss_file)
    if (len_trim(sss_file) == 0) &
        sss_file = '/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/PHC2_salx.nc'
    allocate(nro(nNodL), nss(nNodL))
    if (partit%npes == 1) then
        ! 1-rank: omit partit (read_other_NetCDF then counts mesh%nod2D — the L42 trap fix).
        call read_other_NetCDF(trim(runoff_file), 'Foxx_o_roff', 1, nro, .false., .true., mesh)
        call read_other_NetCDF(trim(sss_file), 'SALT', 1, nss, .true., .true., mesh)
    else
        call read_other_NetCDF(trim(runoff_file), 'Foxx_o_roff', 1, nro, .false., .true., mesh, partit)
        call read_other_NetCDF(trim(sss_file), 'SALT', 1, nss, .true., .true., mesh, partit)
    end if
    nro = nro / 1000.0_WP
    if (partit%mype == 0) &
        write(*,'(a,2es12.4,a,2es12.4,a)') 'fesom_lifecycle_native_mr: native runoff[', &
            minval(nro(1:nNodO)), maxval(nro(1:nNodO)), '] Ssurf[', &
            minval(nss(1:nNodO)), maxval(nss(1:nNodO)), '] (rank0 owned)'

    !===========================================================================
    ! ocean-step surface BC arrays. heat_flux/water_flux/virtual_salt/relax_salt/stress_surf
    ! are produced NATIVELY each step into atm%* / stress_surf; Ki/real_salt_flux/is_nonlinfs
    ! stay 0 (M4 / linfs). Local sizes.
    is_nonlinfs = 0.0_WP
    allocate(Ki(nl-1, nNodL), real_salt_flux(nNodL), stress_surf(2, nElemF))
    Ki = 0.0_WP; real_salt_flux = 0.0_WP; stress_surf = 0.0_WP

    !===========================================================================
    ! per-rank dump (mod_dump: gid-keyed probes, the rank owning a probe writes it).
    allocate(idlist(max(nNodO, nElemO)))
    do n = 1, size(idlist); idlist(n) = n; end do
    if (partit%npes == 1) then
        call dump_init(partit%mype, nNodO, idlist(1:nNodO), nElemO, idlist(1:nElemO))
    else
        call dump_init(partit%mype, nNodO, partit%myList_nod2D, nElemO, partit%myList_elem2D)
    end if

    !===========================================================================
    ! runloop: ocean2ice -> [native atm] -> ice_timestep -> oce_fluxes_mom -> oce_fluxes ->
    ! step_oce (the FESOM2 runloop order; update_atm_forcing replaced by apply_native_forcing).
    do n = 1, nsteps
        call ocean2ice(ice, dyn, tracers, mesh, partit)
        ! native atmosphere (after ocean2ice -> srfoce live; before EVP -> ice%uice/vice are
        ! the previous step's, as FESOM2 update_atm_forcing uses for the wind-on-ice stress).
        call apply_native_forcing(n)
        call ice_timestep(ice, mesh, atm, partit)
        call oce_fluxes_mom(ice, atm, stress_surf, mesh, partit)
        call oce_fluxes(ice, tracers, atm, mesh, partit)
        call step_oce(n, dt, (n == 1), dyn, tracers, mesh, Ki, &
                      atm%heat_flux, atm%water_flux, atm%virtual_salt, atm%relax_salt, &
                      real_salt_flux, is_nonlinfs, stress_surf, partit)
        if (partit%mype == 0) &
            write(*,'(a,i0,a,es12.4,a,es12.4,a,es12.4)') 'fesom_lifecycle_native_mr: step ', n, &
                '  max|eta_n(owned)|=', maxval(abs(dyn%eta_n(1:nNodO))), &
                '  max|uv(owned)|=', maxval(abs(dyn%uv(:,:,1:nElemO))), &
                '  max|a_ice(owned)|=', maxval(abs(ice%data(1)%values(1:nNodO)))
    end do

    call dump_finalize()
    if (partit%mype == 0) &
        write(*,'(a,i0,a)') 'fesom_lifecycle_native_mr: done (', nsteps, ' steps, NATIVE fluxes).'
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)

contains

    ! compute the native CORE2 atmosphere over owned+halo (nNodL): the 8 NCAR fields
    ! (timeinterp at the per-step rdate) + the NCAR bulk (Ch/Ce, live srfoce) + stress_atmoce
    ! + the wind-on-ice stress (previous-step uice/vice). runoff (nro) + Ssurf (nss) were read
    ! once at setup. Every routine threads partit -> the OWNED+HALO nodes are filled.
    subroutine compute_native_forcing(n)
        integer, intent(in) :: n
        real(kind=WP) :: tnew, rcur
        tnew = real(n,WP)*dt                          ! daynew=1 for n<48 (within forcing day 1)
        rcur = real(forcing_julday(fyear,1,1,'noleap'),WP) + real(1-1,WP) &
             + tnew/86400._WP - dt/86400._WP/2._WP
        call forcing_timeinterp(frc, rcur, partit)
        nuw = frc%atmdata(1,1:nNodL);  nvw = frc%atmdata(2,1:nNodL);  nsh = frc%atmdata(3,1:nNodL)
        nswr = frc%atmdata(4,1:nNodL);  nlw = frc%atmdata(5,1:nNodL)
        nta = frc%atmdata(6,1:nNodL) - 273.15_WP
        npr = frc%atmdata(7,1:nNodL) / 1000._WP
        nps = frc%atmdata(8,1:nNodL) / 1000._WP
        ncd = 0.0_WP; nch = 0.0_WP; nce = 0.0_WP
        call forcing_bulk_ncar(10.0_WP, 10.0_WP, 10.0_WP, nta, nsh, nuw, nvw, &
                               ice%srfoce_temp, ice%srfoce_u, ice%srfoce_v, ncd, nch, nce, mesh, partit)
        call forcing_wind_stress(0.0_WP, nuw, nvw, ice%srfoce_u, ice%srfoce_v, ncd, nsx, nsy, mesh, partit)
        call forcing_ice_stress(0.0012_WP, nuw, nvw, ice%uice, ice%vice, nix, niy, mesh, partit)
    end subroutine compute_native_forcing

    ! compute natively + WRITE into atm%* / ice%stress_atmice for the step.
    subroutine apply_native_forcing(n)
        integer, intent(in) :: n
        call compute_native_forcing(n)
        atm%shortwave       = nswr
        atm%longwave        = nlw
        atm%Tair            = nta
        atm%shum            = nsh
        atm%prec_rain       = npr
        atm%prec_snow       = nps
        atm%u_wind          = nuw
        atm%v_wind          = nvw
        atm%Ch_atm_oce_arr  = nch
        atm%Ce_atm_oce_arr  = nce
        atm%stress_atmoce_x = nsx
        atm%stress_atmoce_y = nsy
        ice%stress_atmice_x(1:nNodL) = nix
        ice%stress_atmice_y(1:nNodL) = niy
        atm%runoff          = nro
        atm%Ssurf           = nss
    end subroutine apply_native_forcing

    subroutine alloc_atm(a, nn)
        type(t_atmflux), intent(inout) :: a
        integer,         intent(in)    :: nn
        allocate(a%shortwave(nn), a%longwave(nn), a%Tair(nn), a%shum(nn))
        allocate(a%prec_rain(nn), a%prec_snow(nn), a%runoff(nn))
        allocate(a%u_wind(nn), a%v_wind(nn), a%Ch_atm_oce_arr(nn), a%Ce_atm_oce_arr(nn))
        allocate(a%evaporation(nn), a%ice_sublimation(nn), a%flice(nn), a%real_salt_flux(nn))
        allocate(a%fw_ice(nn), a%fw_snw(nn))
        allocate(a%hf_Qlat(nn), a%hf_Qsen(nn), a%hf_Qradtot(nn))
        allocate(a%hf_Qswr(nn), a%hf_Qlwr(nn), a%hf_Qlwrout(nn))
        allocate(a%heat_flux(nn), a%water_flux(nn), a%heat_flux_in(nn))
        allocate(a%virtual_salt(nn), a%relax_salt(nn))
        allocate(a%stress_node_surf(2, nn))
        allocate(a%stress_atmoce_x(nn), a%stress_atmoce_y(nn), a%Ssurf(nn))
        a%shortwave = 0.0_WP; a%longwave = 0.0_WP; a%Tair = 0.0_WP; a%shum = 0.0_WP
        a%prec_rain = 0.0_WP; a%prec_snow = 0.0_WP; a%runoff = 0.0_WP
        a%u_wind = 0.0_WP; a%v_wind = 0.0_WP; a%Ch_atm_oce_arr = 0.0_WP; a%Ce_atm_oce_arr = 0.0_WP
        a%evaporation = 0.0_WP; a%ice_sublimation = 0.0_WP; a%flice = 0.0_WP; a%real_salt_flux = 0.0_WP
        a%fw_ice = 0.0_WP; a%fw_snw = 0.0_WP
        a%hf_Qlat = 0.0_WP; a%hf_Qsen = 0.0_WP; a%hf_Qradtot = 0.0_WP
        a%hf_Qswr = 0.0_WP; a%hf_Qlwr = 0.0_WP; a%hf_Qlwrout = 0.0_WP
        a%heat_flux = 0.0_WP; a%water_flux = 0.0_WP; a%heat_flux_in = 0.0_WP
        a%virtual_salt = 0.0_WP; a%relax_salt = 0.0_WP
        a%stress_node_surf = 0.0_WP
        a%stress_atmoce_x = 0.0_WP; a%stress_atmoce_y = 0.0_WP; a%Ssurf = 0.0_WP
    end subroutine alloc_atm

end program fesom_lifecycle_native_mr
