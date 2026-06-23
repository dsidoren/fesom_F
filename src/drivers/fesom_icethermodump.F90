program fesom_icethermodump
    ! M3d sea-ice thermodynamics byte-gate driver. Extends fesom_icefctdump (M3c): same init
    ! (mesh + do_ic3d IC + ssh_stiff + ice_setup + cold-start ice), then ocean2ice +
    ! EVPdynamics_solve + ice_TG_rhs + ice_fct_solve (the M3b/M3c-proven chain) to get the
    ! advected ice state, then PRESCRIBES the analytic atmospheric forcing (shortwave/longwave/
    ! Tair/shum/prec/wind/Ch-Ce_atm_oce), runs the M3d kernels cut_off + thermodynamics, and
    ! dumps a_ice/m_ice/m_snow + flx_h/flx_fw + t_skin. tools/run_icethermo_gate_core2.sh
    ! compares vs the FESOM2 oracle (src/fesom_ice_dump.F90::ice_thermo_dump_write) for
    ! max|delta|=0, CORE2 1-rank, BOTH whichEVP=0 and 1 (the advection upstream).
    !
    ! CONFIG OVERRIDES (the oracle reads these from the CORE2 namelists; the driver matches):
    !   - &ice_dyn (as M3c): cd_oce_ice=0.0055 / delta_min=1.0e-11 / ice_diff=0.0 /
    !     ice_gamma_fct=0.5 (namelist values, NOT t_ice single->WP defaults).
    !   - &ice_therm: the oracle reads these as namelist DOUBLES; the t_ice_thermo defaults are
    !     single->WP, so they MUST be re-set as _WP doubles to byte-match (con=2.1656_WP,
    !     hmin/Armin=0.01_WP, emiss=0.97_WP, albedos, albw=0.1_WP [not 0.066], Sice=4.0_WP,
    !     h0/h0_s=0.5_WP, c_melt=0.5_WP, con/consn). h_ml=2.5_WP / cc / cl / rho* / clh* /
    !     tmelt / boltzmann / cpair are NOT in the namelist -> the MOD_ICE single->WP defaults
    !     match on both sides (cc=rhowat*4190 / cl=rhoice*3.34e5 are exactly representable).
    !   - atmflux scalars: Ch_atm_ice=Ce_atm_ice=0.00175_WP (namelist.forcing doubles, NOT
    !     1.75e-3 single->WP); ref_sss_local=.true. (=> rsss=S_oc); use_virt_salt=.true.
    !     (linfs); l_snow=.true. (=> rain=prec_rain, snow=prec_snow).
    !
    !   FESOM3_MESH_DIR  mesh dir   (default: CORE2)
    !   FESOM3_IC_FILE   IC netcdf  (default: pool phc3.0_winter.nc — same file the oracle reads)
    !   FESOM3_WHICHEVP  0=EVP / 1=mEVP (default 0)
    !   FESOM3_THERMO_OUT out path  (default: icethermo_f3.bin)
    use mpi
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use mod_precision,      only: WP, MP
    use mod_param_phys,     only: alpha, theta
    use mod_mesh,           only: t_mesh
    use mod_partit,         only: t_partit
    use mod_partitioning,   only: par_init, par_ex
    use mod_mesh_read,      only: read_mesh
    use mod_mesh_areas,     only: compute_geometry
    use mod_dyn,            only: t_dyn
    use mod_tracer,         only: t_tracer
    use mod_ice,            only: t_ice
    use mod_ice_setup,      only: ice_setup
    use mod_ice_dyn,        only: ocean2ice, EVPdynamics_solve
    use mod_ice_fct,        only: ice_TG_rhs, ice_fct_solve
    use mod_ice_thermo,     only: t_atmflux, cut_off, thermodynamics
    use oce_initial_state,  only: t_ic3d_config, do_ic3d
    use oce_ssh_rhs,        only: init_stiff_mat_ale
    use mod_advhor_dump,    only: advhor_dump_open, advhor_dump_close, wr_r1
    implicit none

    ! CORE2 namelist timestep: dt = 86400/48 = 1800 s (computed like FESOM2, not a literal).
    real(kind=WP), parameter :: dt = 86400.0_WP / real(48, WP)

    character(len=512) :: mesh_dir, out_path, ic_file, whichevp_str
    type(t_partit)     :: partit
    type(t_mesh)       :: mesh
    type(t_dyn)        :: dyn
    type(t_tracer)     :: tracers
    type(t_ice)        :: ice
    type(t_atmflux)    :: atm
    type(t_ic3d_config):: ic
    integer :: nsw, n, nl, nzmin, nzmax, e, u, whichevp, ios
    real(kind=WP) :: lon, lat
    real(kind=MP) :: zbar_srf, zbar_bot

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2'
    call get_environment_variable('FESOM3_IC_FILE', ic_file)
    if (len_trim(ic_file) == 0) &
        ic_file = '/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc'
    call get_environment_variable('FESOM3_THERMO_OUT', out_path)
    if (len_trim(out_path) == 0) out_path = 'icethermo_f3.bin'
    whichevp = 0
    call get_environment_variable('FESOM3_WHICHEVP', whichevp_str)
    if (len_trim(whichevp_str) > 0) read(whichevp_str, *, iostat=ios) whichevp

    call par_init(partit)
    if (partit%npes /= 1) then
        if (partit%mype == 0) write(*,'(a)') 'fesom_icethermodump: requires 1 rank'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype); error stop 1
    end if

    !===========================================================================
    ! mesh + geometry (CORE2 rotation 50/15/-90, identical to the lifecycle/IC/ice gate).
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    nl = mesh%nl
    write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') 'fesom_icethermodump: nod2D=', mesh%nod2D, &
        ' elem2D=', mesh%elem2D, ' edge2D=', mesh%edge2D, ' nl=', nl, ' CW swaps=', nsw

    !===========================================================================
    ! ALE depths Z_3d_n (do_ic3d) + zbar_e_bot (init_stiff_mat_ale). linfs full cells,
    ! eta=0 — identical to fesom_icefctdump.
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
    allocate(mesh%zbar_e_bot(mesh%elem2D)); mesh%zbar_e_bot = 0.0_MP
    do e = 1, mesh%elem2D
        mesh%zbar_e_bot(e) = mesh%zbar(mesh%nlevels(e))
    end do
    allocate(mesh%hbar(mesh%nod2D)); mesh%hbar = 0.0_MP

    !===========================================================================
    ! 2-tracer ocean state + do_ic3d phc3.0 IC (the surface T/S the cold start reads + the
    ! surface T/S ocean2ice copies into srfoce_temp/salt -> the thermo's T_oc/S_oc).
    tracers%num_tracers = 2
    allocate(tracers%data(2))
    allocate(tracers%data(1)%values(nl-1, mesh%nod2D), tracers%data(2)%values(nl-1, mesh%nod2D))
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
    call do_ic3d(tracers, ic, mesh)
    write(*,'(a,2es12.4)') 'fesom_icethermodump: surface T range = ', &
        minval(tracers%data(1)%values(1,:)), maxval(tracers%data(1)%values(1,:))

    !===========================================================================
    ! ssh_stiff CSR (reduced-M2: alpha=theta=1) — the ice mass matrix rides its sparsity.
    alpha = 1.0_WP; theta = 1.0_WP
    call init_stiff_mat_ale(mesh, dt)

    !===========================================================================
    ! ice setup: allocate + ice_mass_matrix_fill + cold-start ice_initial_state. Then the
    ! &ice_dyn override (M3c) AND the &ice_therm override (M3d, namelist doubles).
    call ice_setup(ice, tracers, mesh, dt)
    ice%whichEVP      = whichevp        ! 0 = standard EVP, 1 = mEVP (FESOM3_WHICHEVP)
    ice%cd_oce_ice    = 0.0055_WP       ! &ice_dyn namelist (double)
    ice%delta_min     = 1.0e-11_WP      ! &ice_dyn namelist (double)
    ice%ice_diff      = 0.0_WP          ! &ice_dyn namelist (no artificial diffusion)
    ice%ice_gamma_fct = 0.5_WP          ! &ice_dyn namelist
    ! &ice_therm namelist DOUBLES (the oracle's namelist read overwrites the type defaults;
    ! these single->WP defaults would differ ~1 ULP, so re-set as _WP doubles to byte-match).
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
    ice%thermo%albw              = 0.1_WP    ! CORE2 namelist (NOT the 0.066 LY2004 default)
    ice%thermo%open_water_albedo = 0
    ice%thermo%con               = 2.1656_WP
    ice%thermo%consn             = 0.31_WP
    ice%thermo%snowdist          = .true.
    ice%thermo%c_melt            = 0.5_WP
    write(*,'(a,i0,a,3es12.4)') 'fesom_icethermodump: whichEVP=', whichevp, &
        '  max a_ice/m_ice/m_snow (IC) = ', &
        maxval(ice%data(1)%values), maxval(ice%data(2)%values), maxval(ice%data(3)%values)

    !===========================================================================
    ! atmflux: allocate (nod2D) + scalar config. The atm scalars match the CORE2 namelists.
    call alloc_atm(atm, mesh%nod2D)
    atm%Ch_atm_ice    = 0.00175_WP      ! namelist.forcing double (NOT 1.75e-3 single->WP)
    atm%Ce_atm_ice    = 0.00175_WP
    atm%ref_sss       = 34.0_WP         ! namelist.tra (dead — ref_sss_local=.true. => rsss=S_oc)
    atm%ref_sss_local = .true.
    atm%use_virt_salt = .true.          ! which_ALE='linfs'
    atm%l_snow        = .true.          ! namelist.forcing

    !===========================================================================
    ! prescribe the analytic surface ocean velocity (level 1; ocean2ice reads UV(:,1,:)).
    allocate(dyn%uv(2, nl-1, mesh%elem2D)); dyn%uv = 0.0_WP
    do e = 1, mesh%elem2D
        lon = mesh%coord_nod2D(1, mesh%elem2D_nodes(1,e))
        lat = mesh%coord_nod2D(2, mesh%elem2D_nodes(1,e))
        dyn%uv(1,1,e) =  0.20_WP*cos(lat)*sin(lon)
        dyn%uv(2,1,e) = -0.15_WP*sin(lat)*cos(2.0_WP*lon)
    end do

    !===========================================================================
    ! prescribe elevation (hbar), wind-on-ice stress, AND the atmospheric forcing arrays.
    ! ALL formulas MUST equal the FESOM2 oracle src/fesom_ice_dump.F90::ice_thermo_dump_write.
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        mesh%hbar(n)              =  0.5_WP*sin(lon)*cos(lat) - 0.2_WP*cos(2.0_WP*lat)
        ice%stress_atmice_x(n)    =  0.10_WP*cos(lat)*sin(lon)
        ice%stress_atmice_y(n)    = -0.08_WP*sin(lat)*cos(2.0_WP*lon)
        atm%shortwave(n)      = 150.0_WP + 100.0_WP*cos(lat)*cos(lon)
        atm%longwave(n)       = 250.0_WP +  50.0_WP*sin(lat)
        atm%Tair(n)           = -10.0_WP +  20.0_WP*cos(lat) + 5.0_WP*sin(2.0_WP*lon)
        atm%shum(n)           =   0.003_WP + 0.002_WP*cos(lat)
        atm%prec_rain(n)      =   1.0e-8_WP*(1.0_WP + 0.5_WP*sin(lon))
        atm%prec_snow(n)      =   2.0e-8_WP*(1.0_WP + 0.5_WP*cos(lat))
        atm%runoff(n)         =   0.0_WP
        atm%u_wind(n)         =   6.0_WP*cos(lat)*sin(lon)
        atm%v_wind(n)         =  -5.0_WP*sin(lat)*cos(2.0_WP*lon)
        atm%Ch_atm_oce_arr(n) =   1.2e-3_WP + 1.0e-4_WP*cos(lat)
        atm%Ce_atm_oce_arr(n) =   1.2e-3_WP + 1.0e-4_WP*sin(lat)
    end do

    !===========================================================================
    ! ocean -> ice coupling + EVP dynamics (M3b) + ice FCT advection (M3c) -> advected ice
    ! state + the friction-velocity inputs (uice/vice, srfoce_u/v/temp/salt).
    call ocean2ice(ice, dyn, tracers, mesh)
    call EVPdynamics_solve(ice, mesh)
    call ice_TG_rhs(ice, mesh)
    call ice_fct_solve(ice, mesh)
    write(*,'(a,3es12.4)') 'fesom_icethermodump: post-adv max a_ice/m_ice/m_snow = ', &
        maxval(ice%data(1)%values(1:mesh%nod2D)), &
        maxval(ice%data(2)%values(1:mesh%nod2D)), &
        maxval(ice%data(3)%values(1:mesh%nod2D))

    !===========================================================================
    ! M3d kernels: cut_off (hmin/Armin clamp) + thermodynamics (growth/melt + air-sea fluxes).
    call cut_off(ice, mesh)
    call thermodynamics(ice, mesh, atm)
    write(*,'(a,3es12.4)') 'fesom_icethermodump: post-thermo max m_ice/flx_h/t_skin = ', &
        maxval(ice%data(2)%values(1:mesh%nod2D)), &
        maxval(abs(ice%flx_h(1:mesh%nod2D))), &
        maxval(abs(ice%thermo%t_skin(1:mesh%nod2D)))

    !===========================================================================
    ! dump (FADVHDMP, 1-rank identity ids). a/m_ice/m_snow + flx_h/flx_fw + t_skin (all nod2D).
    call advhor_dump_open(u, trim(out_path), mesh%nod2D, mesh%elem2D, mesh%edge2D, nl)
    call wr_r1(u, 'ice_a_ice',  ice%data(1)%values(1:mesh%nod2D))
    call wr_r1(u, 'ice_m_ice',  ice%data(2)%values(1:mesh%nod2D))
    call wr_r1(u, 'ice_m_snow', ice%data(3)%values(1:mesh%nod2D))
    call wr_r1(u, 'ice_flx_h',  ice%flx_h(1:mesh%nod2D))
    call wr_r1(u, 'ice_flx_fw', ice%flx_fw(1:mesh%nod2D))
    call wr_r1(u, 'ice_t_skin', ice%thermo%t_skin(1:mesh%nod2D))
    call advhor_dump_close(u)
    write(*,'(a)') 'fesom_icethermodump: wrote '//trim(out_path)

    call par_ex(partit%MPI_COMM_FESOM, partit%mype)

contains

    ! allocate + zero all t_atmflux per-node arrays (sized nod2D for the 1-rank gate).
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
        a%shortwave = 0.0_WP; a%longwave = 0.0_WP; a%Tair = 0.0_WP; a%shum = 0.0_WP
        a%prec_rain = 0.0_WP; a%prec_snow = 0.0_WP; a%runoff = 0.0_WP
        a%u_wind = 0.0_WP; a%v_wind = 0.0_WP; a%Ch_atm_oce_arr = 0.0_WP; a%Ce_atm_oce_arr = 0.0_WP
        a%evaporation = 0.0_WP; a%ice_sublimation = 0.0_WP; a%flice = 0.0_WP; a%real_salt_flux = 0.0_WP
        a%fw_ice = 0.0_WP; a%fw_snw = 0.0_WP
        a%hf_Qlat = 0.0_WP; a%hf_Qsen = 0.0_WP; a%hf_Qradtot = 0.0_WP
        a%hf_Qswr = 0.0_WP; a%hf_Qlwr = 0.0_WP; a%hf_Qlwrout = 0.0_WP
    end subroutine alloc_atm

end program fesom_icethermodump
