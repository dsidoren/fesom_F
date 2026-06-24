program fesom_lifecycle
    ! M2.11c lifecycle byte-gate driver — the FIRST real multi-step time-stepping run
    ! in the port (not a prescribe-and-stop shim). Mirrors FESOM2's fesom_init ->
    ! fesom_runloop(N) -> fesom_finalize on the CORE2 mesh, 1-rank, reduced-M2 (linfs /
    ! PP / no-GM / no-Redi / opt_visc=7 / use_wsplit=.false.):
    !
    !   model_init : par_init -> read CORE2 mesh (rotation 50/15/-90) -> compute_geometry
    !                -> ALE depths (linfs full cells) -> do_ic3d phc3.0 IC (the M2.11b read)
    !                -> cold-start state (UV/eta/w/uv_rhsAB = 0; valuesold = values)
    !                -> work/FCT arrays + muscl_adv_init -> SSH stiffness -> reduced-M2 config
    !   runloop    : for n=1..N call step_oce(n, dt, lfirst=(n==1), ...) with the SAME
    !                persistent dyn/tracers/mesh state (multi-step AB2 velocity + tracer
    !                rotation evolve in place). step_oce emits the per-substep node dumps
    !                (mod_dump) the byte-gate compares against the REAL FESOM2 oce_timestep_ale.
    !   finalize   : dump_finalize + par_ex.
    !
    ! M2.11c-1 (this gate) is UNFORCED: the oracle runs use_ice=.false. so all surface
    ! fluxes stay 0; FESOM3 prescribes Ki/heat_flux/water_flux/virtual_salt/relax_salt/
    ! stress_surf = 0. M2.11c-2 prescribes the oracle's per-step fluxes from a dump file
    ! (FESOM3_FLUX_FILE) — that hook is read_step_fluxes below.
    !
    !   FESOM3_MESH_DIR   mesh dir         (default: CORE2)
    !   FESOM3_IC_FILE    IC netcdf        (default: pool phc3.0_winter.nc)
    !   FESOM_DUMP_FILE   node dump prefix (mod_dump; set by run script)
    !   FESOM_DUMP_MAXSTEPS  dump step cap (mod_dump)
    !   FESOM3_NSTEPS     number of steps  (default: 3)
    use mpi
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use mod_precision,      only: WP, MP
    use mod_constants,      only: density_0
    use mod_param_phys,     only: N2smth_h, alpha, theta
    use mod_param_phys,     only: mix_coeff_PP, A_ver, K_ver, Kv0_const
    use mod_param_phys,     only: use_instabmix, instabmix_kv, use_momix, use_windmix
    use mod_param_phys,     only: Fer_GM, Redi, K_GM_max, K_GM_min, K_GM_bvref, &
                                  K_GM_rampmax, K_GM_rampmin, K_GM_resscalorder, K_GM_cm, &
                                  K_GM_cmin, K_GM_Ktaper, scaling_Ferreira, scaling_Rossby, &
                                  scaling_resolution, scaling_FESOM14, scaling_GMzexp, &
                                  scaling_GINsea, GMzexp_zref, GMzexp_smin
    use mod_param_phys,     only: Redi_Kmax, Redi_Kmin, Redi_Ktaper, K_hor, &
                                  scaling_ODM95, ODM95_Scr, ODM95_Sd, scaling_LDD97
    use mod_mesh,           only: t_mesh
    use mod_partit,         only: t_partit
    use mod_partitioning,   only: par_init, par_ex
    use mod_mesh_read,      only: read_mesh
    use mod_mesh_areas,     only: compute_geometry
    use mod_dyn,            only: t_dyn
    use mod_tracer,         only: t_tracer
    use oce_initial_state,  only: t_ic3d_config, do_ic3d
    use oce_muscl_adv,      only: muscl_adv_init
    use oce_ssh_rhs,        only: init_stiff_mat_ale
    use mod_step_oce,       only: step_oce
    use mod_dump,           only: dump_init, dump_finalize
    implicit none

    ! CORE2 namelist timestep: dt = 86400/step_per_day, step_per_day=48 -> 1800 s.
    ! Computed the SAME way FESOM2 does (gen_model_setup.F90:92), NOT an 1800.0 literal.
    real(kind=WP), parameter :: dt = 86400.0_WP / real(48, WP)

    character(len=512) :: mesh_dir, ic_file, env
    type(t_partit)     :: partit
    type(t_mesh)       :: mesh
    type(t_dyn)        :: dyn
    type(t_tracer)     :: tracers
    type(t_ic3d_config):: ic
    integer :: nsw, n, nz, nl, nzmin, nzmax, e, tr_num, nsteps, ios, env_len
    real(kind=MP) :: zbar_srf, zbar_bot
    real(kind=WP), allocatable :: Ki(:,:), heat_flux(:), water_flux(:), virtual_salt(:)
    real(kind=WP), allocatable :: relax_salt(:), real_salt_flux(:), stress_surf(:,:)
    real(kind=WP) :: is_nonlinfs
    integer, allocatable :: idlist(:)
    ! M2.11c-2 forced: read the oracle's per-step surface fluxes (the M3 air-sea gap) from
    ! FESOM3_FLUX_FILE (written by the oracle fesom_flux_dump shim). If unset -> unforced (0).
    character(len=512) :: flux_file
    logical :: forced
    integer :: flux_unit
    integer(int32) :: fstep, fnn, fne
    ! M4c/M4d: enable GM bolus (FESOM3_FER_GM) / Redi isopycnal diffusion (FESOM3_REDI).
    logical :: use_fer_gm, use_redi

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2'
    call get_environment_variable('FESOM3_IC_FILE', ic_file)
    if (len_trim(ic_file) == 0) &
        ic_file = '/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc'
    nsteps = 3
    call get_environment_variable('FESOM3_NSTEPS', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) read(env, *, iostat=ios) nsteps
    call get_environment_variable('FESOM3_FER_GM', env, length=env_len, status=ios)
    use_fer_gm = (ios == 0 .and. env_len > 0)
    call get_environment_variable('FESOM3_REDI', env, length=env_len, status=ios)
    use_redi = (ios == 0 .and. env_len > 0)

    call par_init(partit)
    if (partit%npes /= 1) then
        if (partit%mype == 0) write(*,'(a)') 'fesom_lifecycle: requires 1 rank'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype); error stop 1
    end if

    !===========================================================================
    ! model_init: mesh + geometry
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    nl = mesh%nl
    write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') 'fesom_lifecycle: nod2D=', mesh%nod2D, &
        ' elem2D=', mesh%elem2D, ' edge2D=', mesh%edge2D, ' nl=', nl, ' CW swaps=', nsw

    !===========================================================================
    ! ALE depth/thickness state (linfs full cells) — identical to fesom_stepdump.
    allocate(mesh%hnode(nl-1, mesh%nod2D)); mesh%hnode = 0.0_MP
    do n = 1, mesh%nod2D
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n)-1
            mesh%hnode(nz, n) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
    end do
    allocate(mesh%helem(nl-1, mesh%elem2D)); mesh%helem = 0.0_MP
    allocate(mesh%zbar_e_bot(mesh%elem2D)); mesh%zbar_e_bot = 0.0_MP
    do e = 1, mesh%elem2D
        do nz = mesh%ulevels(e), mesh%nlevels(e)-1
            mesh%helem(nz, e) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
        mesh%zbar_e_bot(e) = mesh%zbar(mesh%nlevels(e))
    end do
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
    ! ALE elevation state — COLD START (eta=hbar=0).
    allocate(mesh%hbar(mesh%nod2D), mesh%hbar_old(mesh%nod2D), mesh%dhe(mesh%elem2D))
    allocate(mesh%hnode_new(nl-1, mesh%nod2D))
    mesh%hbar = 0.0_MP; mesh%hbar_old = 0.0_MP; mesh%dhe = 0.0_MP
    mesh%hnode_new = mesh%hnode            ! linfs: hnode_new == hnode

    !===========================================================================
    ! dynamics state — COLD START (UV/eta_n/w/uv_rhsAB/d_eta/ssh_rhs_old = 0).
    allocate(dyn%uv(2, nl-1, mesh%elem2D), dyn%uv_rhs(2, nl-1, mesh%elem2D))
    allocate(dyn%uv_rhsAB(1, 2, nl-1, mesh%elem2D))     ! (AB_order-1, 2, nl-1, elem2D)
    allocate(dyn%uvnode(2, nl-1, mesh%nod2D))
    allocate(dyn%eta_n(mesh%nod2D), dyn%d_eta(mesh%nod2D))
    allocate(dyn%ssh_rhs(mesh%nod2D), dyn%ssh_rhs_old(mesh%nod2D))
    allocate(dyn%w(nl, mesh%nod2D), dyn%w_e(nl, mesh%nod2D), dyn%w_i(nl, mesh%nod2D))
    allocate(dyn%cfl_z(nl, mesh%nod2D))
    allocate(dyn%work%density_ref(nl-1, mesh%nod2D), dyn%work%density_m_rho0(nl-1, mesh%nod2D))
    allocate(dyn%work%hpressure(nl, mesh%nod2D), dyn%work%bvfreq(nl, mesh%nod2D))
    allocate(dyn%work%pgf_x(nl-1, mesh%elem2D), dyn%work%pgf_y(nl-1, mesh%elem2D))
    allocate(dyn%work%u_c(nl-1, mesh%elem2D), dyn%work%v_c(nl-1, mesh%elem2D))
    allocate(dyn%work%uvnode_rhs(2, nl-1, mesh%nod2D))
    allocate(dyn%work%Kv(nl, mesh%nod2D), dyn%work%Av(nl, mesh%elem2D))
    dyn%uv = 0.0_WP; dyn%uv_rhs = 0.0_WP; dyn%uv_rhsAB = 0.0_WP; dyn%uvnode = 0.0_WP
    dyn%eta_n = 0.0_WP; dyn%d_eta = 0.0_WP
    dyn%ssh_rhs = 0.0_WP; dyn%ssh_rhs_old = 0.0_WP
    dyn%w = 0.0_WP; dyn%w_e = 0.0_WP; dyn%w_i = 0.0_WP; dyn%cfl_z = 0.0_WP
    dyn%work%density_m_rho0 = 0.0_WP; dyn%work%hpressure = 0.0_WP; dyn%work%bvfreq = 0.0_WP
    dyn%work%pgf_x = 0.0_WP; dyn%work%pgf_y = 0.0_WP
    dyn%work%u_c = 0.0_WP; dyn%work%v_c = 0.0_WP; dyn%work%uvnode_rhs = 0.0_WP
    dyn%work%Kv = 0.0_WP; dyn%work%Av = 0.0_WP
    dyn%work%density_ref = density_0              ! use_density_ref=.false.
    dyn%AB_order      = 2
    dyn%momadv_opt    = 2
    dyn%opt_visc      = 7
    dyn%visc_gamma0   = 0.003_WP
    dyn%visc_gamma1   = 0.1_WP
    dyn%visc_gamma2   = 0.285_WP
    dyn%visc_gamma0_h = 0.0_WP
    dyn%visc_gamma1_h = 0.0_WP
    dyn%use_wsplit    = .false.          ! CORE2 production (and M1.4/M2.9b precedent)
    dyn%wsplit_maxcfl = 1.0_WP

    !===========================================================================
    ! 2-tracer state (data(1)=T ID 1, data(2)=S ID 2) + do_ic3d phc3.0 IC.
    tracers%num_tracers = 2
    allocate(tracers%data(2))
    allocate(tracers%data(1)%values(nl-1, mesh%nod2D), tracers%data(2)%values(nl-1, mesh%nod2D))
    tracers%data(1)%values = 0.0_WP;  tracers%data(1)%ID = 1
    tracers%data(2)%values = 0.0_WP;  tracers%data(2)%ID = 2
    ! IC config = work_core/namelist.tra &tracer_init3d (idlist=2,1 -> salt first, temp second)
    ic%n_ic3d      = 2
    ic%idlist(1:2) = [2, 1]
    ic%filelist(1) = trim(ic_file)
    ic%filelist(2) = trim(ic_file)
    ic%varlist(1)  = 'salt'
    ic%varlist(2)  = 'temp'
    ic%t_insitu    = .true.
    ic%ic_cyclic   = .true.
    ic%dummy       = 1.e10_WP
    call do_ic3d(tracers, ic, mesh)
    write(*,'(a,2es12.4,a,2es12.4)') 'fesom_lifecycle: IC T=', &
        minval(tracers%data(1)%values), maxval(tracers%data(1)%values), &
        '  S=', minval(tracers%data(2)%values), maxval(tracers%data(2)%values)

    !===========================================================================
    ! tracer advection machinery: valuesAB / valuesold (cold start = values) +
    ! FCT work arrays + muscl_adv_init, pi/CORE2 FCT config (MFCT/QR4C/FCT, opth=0/optv=1).
    do tr_num = 1, 2
        allocate(tracers%data(tr_num)%valuesAB(nl-1, mesh%nod2D))
        allocate(tracers%data(tr_num)%valuesold(2, nl-1, mesh%nod2D))
        tracers%data(tr_num)%AB_order   = 2
        tracers%data(tr_num)%tra_adv_hor = 'MFCT'
        tracers%data(tr_num)%tra_adv_ver = 'QR4C'
        tracers%data(tr_num)%tra_adv_lim = 'FCT'
        tracers%data(tr_num)%tra_adv_ph  = 0.0_WP    ! MFCT
        tracers%data(tr_num)%tra_adv_pv  = 1.0_WP    ! QR4C
        tracers%data(tr_num)%i_vert_diff = .true.
        tracers%data(tr_num)%valuesAB        = 0.0_WP
        tracers%data(tr_num)%valuesold(1,:,:) = tracers%data(tr_num)%values   ! cold start
        tracers%data(tr_num)%valuesold(2,:,:) = tracers%data(tr_num)%values
    end do
    allocate(tracers%work%fct_LO          (nl-1, mesh%nod2D))
    allocate(tracers%work%adv_flux_hor    (nl-1, mesh%edge2D))
    allocate(tracers%work%adv_flux_ver    (nl,   mesh%nod2D))
    allocate(tracers%work%fct_ttf_max     (nl-1, mesh%nod2D))
    allocate(tracers%work%fct_ttf_min     (nl-1, mesh%nod2D))
    allocate(tracers%work%fct_plus        (nl-1, mesh%nod2D))
    allocate(tracers%work%fct_minus       (nl-1, mesh%nod2D))
    allocate(tracers%work%del_ttf         (nl-1, mesh%nod2D))
    allocate(tracers%work%del_ttf_advhoriz(nl-1, mesh%nod2D))
    allocate(tracers%work%del_ttf_advvert (nl-1, mesh%nod2D))
    call muscl_adv_init(tracers%work, mesh)

    !===========================================================================
    ! surface forcing — M2.11c-1 UNFORCED: all zero (oracle runs use_ice=.false.).
    is_nonlinfs = 0.0_WP
    allocate(Ki(nl-1, mesh%nod2D), heat_flux(mesh%nod2D), water_flux(mesh%nod2D))
    allocate(virtual_salt(mesh%nod2D), relax_salt(mesh%nod2D), real_salt_flux(mesh%nod2D))
    allocate(stress_surf(2, mesh%elem2D))
    Ki = 0.0_WP; heat_flux = 0.0_WP; water_flux = 0.0_WP; virtual_salt = 0.0_WP
    relax_salt = 0.0_WP; real_salt_flux = 0.0_WP; stress_surf = 0.0_WP

    !===========================================================================
    ! reduced-M2 module config (= the FESOM2 oracle / CORE2 namelist).
    alpha = 1.0_WP; theta = 1.0_WP
    N2smth_h     = .true.
    mix_coeff_PP = 0.01_WP
    A_ver        = 1.0e-4_WP
    K_ver        = 1.0e-5_WP
    Kv0_const    = .true.
    use_instabmix = .true.
    instabmix_kv  = 0.1_WP
    use_momix     = .false.
    use_windmix   = .false.

    !===========================================================================
    ! M4c GM (FESOM3_FER_GM set): enable Gent-McWilliams bolus advection. Force the
    ! work_core GM config (Fer_GM=T, Redi=F isolating gate; scaling_resolution + GMzexp,
    ! K_GM_max=1000/cm=3/cmin=0.1/resscalorder=2/ramp off/Ktaper off). Allocate + init the GM
    ! work arrays to the FESOM2 oce_setup_step.F90:962-967 / :670/:688 values (fer_K=500/
    ! fer_c=1/fer_scal=0/fer_gamma=0, fer_uv/fer_w=0) so below-bottom byte-matches. The GM
    ! chain + bolus run inside step_oce / solve_tracers_ale (guarded by Fer_GM).
    if (use_fer_gm .or. use_redi) then
        Fer_GM = .true.        ! Redi rides the GM coupling (Ki = max(fer_scal*Redi_Kmax, K_GM_min))
        K_GM_max = 1000.0_WP; K_GM_min = 2.0_WP; K_GM_bvref = 1
        K_GM_rampmax = -1.0_WP; K_GM_rampmin = -1.0_WP; K_GM_resscalorder = 2.0_WP
        K_GM_cm = 3.0_WP; K_GM_cmin = 0.1_WP; K_GM_Ktaper = .false.
        scaling_Ferreira = .false.; scaling_Rossby = .false.; scaling_resolution = .true.
        scaling_FESOM14 = .false.; scaling_GMzexp = .true.; scaling_GINsea = .false.
        GMzexp_zref = 500.0_WP; GMzexp_smin = 0.6_WP
        allocate(dyn%fer_uv(2, nl-1, mesh%elem2D), dyn%fer_w(nl, mesh%nod2D))
        allocate(dyn%work%sw_alpha(nl-1, mesh%nod2D), dyn%work%sw_beta(nl-1, mesh%nod2D))
        allocate(dyn%work%sigma_xy(2, nl-1, mesh%nod2D))
        allocate(dyn%work%fer_K(nl, mesh%nod2D), dyn%work%fer_c(mesh%nod2D), dyn%work%fer_scal(mesh%nod2D))
        allocate(dyn%work%fer_gamma(2, nl, mesh%nod2D))
        dyn%fer_uv = 0.0_WP; dyn%fer_w = 0.0_WP
        dyn%work%sw_alpha = 0.0_WP; dyn%work%sw_beta = 0.0_WP; dyn%work%sigma_xy = 0.0_WP
        dyn%work%fer_K = 500.0_WP; dyn%work%fer_c = 1.0_WP; dyn%work%fer_scal = 0.0_WP
        dyn%work%fer_gamma = 0.0_WP
        if (use_redi) then
            ! M4d Redi (work_core): Redi=T; Redi_Kmax<=0 -> synced to K_GM_max; Redi_Kmin=100;
            ! Redi_Ktaper sqrt-split; ODM95 slope taper (work_core ODM95_Scr=0.2e-2); K_hor=0
            ! (namelist.tra). Inits = FESOM2 oce_setup_step.F90 (Ki=K_hor*..=0; fer_tapfac=1).
            Redi = .true.; Redi_Ktaper = .true.; Redi_Kmax = 0.0_WP; Redi_Kmin = 100.0_WP
            scaling_ODM95 = .true.; ODM95_Scr = 0.2e-2_WP; ODM95_Sd = 1.0e-3_WP; scaling_LDD97 = .false.
            K_hor = 0.0_WP
            allocate(dyn%work%Ki(nl-1, mesh%nod2D), dyn%work%fer_tapfac(nl-1, mesh%nod2D))
            allocate(dyn%work%neutral_slope(3, nl-1, mesh%nod2D), dyn%work%slope_tapered(3, nl-1, mesh%nod2D))
            allocate(tracers%work%tr_z(nl, mesh%nod2D))
            dyn%work%Ki = 0.0_WP; dyn%work%fer_tapfac = 1.0_WP
            dyn%work%neutral_slope = 0.0_WP; dyn%work%slope_tapered = 0.0_WP
            tracers%work%tr_z = 0.0_WP
            write(*,'(a)') 'fesom_lifecycle: Fer_GM + Redi ENABLED (work_core GM+Redi config)'
        else
            Redi = .false.
            write(*,'(a)') 'fesom_lifecycle: Fer_GM ENABLED (work_core GM config; Redi off)'
        end if
    end if

    ! SSH stiffness (built ONCE; dt = CORE2 namelist timestep).
    call init_stiff_mat_ale(mesh, dt)

    !===========================================================================
    ! per-substep node dump (mod_dump; FESOM_DUMP_FILE set by the run script).
    allocate(idlist(max(mesh%nod2D, mesh%elem2D)))
    do n = 1, size(idlist); idlist(n) = n; end do   ! 1-rank: identity global ids
    call dump_init(partit%mype, mesh%nod2D, idlist(1:mesh%nod2D), &
                   mesh%elem2D, idlist(1:mesh%elem2D))

    !===========================================================================
    ! M2.11c-2 forced: open the oracle's per-step flux dump (FESOM3_FLUX_FILE). Unset
    ! -> unforced (the M2.11c-1 zero-flux gate). The flux record order MUST mirror the
    ! oracle fesom_flux_dump shim: int32 step/nn/ne, then heat_flux/water_flux/virtual_salt/
    ! relax_salt (nod2D) + stress_surf (2,elem2D), all real64 (== WP at the anchor).
    call get_environment_variable('FESOM3_FLUX_FILE', flux_file)
    forced = (len_trim(flux_file) > 0)
    if (forced) then
        open(newunit=flux_unit, file=trim(flux_file), status='old', form='unformatted', &
             access='stream', action='read')
        write(*,'(a,a)') 'fesom_lifecycle: FORCED — prescribing per-step fluxes from ', trim(flux_file)
    end if

    !===========================================================================
    ! runloop: N steps, state persists (multi-step AB2 evolution). lfirst=(n==1).
    do n = 1, nsteps
        if (forced) then
            read(flux_unit) fstep, fnn, fne
            read(flux_unit) heat_flux
            read(flux_unit) water_flux
            read(flux_unit) virtual_salt
            read(flux_unit) relax_salt
            read(flux_unit) stress_surf
        end if
        call step_oce(n, dt, (n == 1), dyn, tracers, mesh, Ki, &
                      heat_flux, water_flux, virtual_salt, relax_salt, &
                      real_salt_flux, is_nonlinfs, stress_surf)
        write(*,'(a,i0,a,es12.4,a,es12.4)') 'fesom_lifecycle: step ', n, &
            '  max|eta_n|=', maxval(abs(dyn%eta_n)), '  max|uv|=', maxval(abs(dyn%uv))
    end do
    if (forced) close(flux_unit)

    call dump_finalize()
    write(*,'(a,i0,a)') 'fesom_lifecycle: done (', nsteps, ' steps + per-substep node dump).'
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_lifecycle
