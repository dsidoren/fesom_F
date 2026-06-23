program fesom_icefluxdump
    ! M3e sea-ice -> ocean coupling-out byte-gate driver (THE PAYOFF). Extends
    ! fesom_icethermodump (M3d): same init (mesh + do_ic3d IC + ssh_stiff + ice_setup +
    ! cold-start ice), then ocean2ice + EVPdynamics_solve + ice_TG_rhs + ice_fct_solve +
    ! cut_off + thermodynamics (the M3b/M3c/M3d-proven chain) to get the advected + thermo'd
    ! ice state + the air-sea fluxes (ice%flx_h/flx_fw, atm%evaporation/...). Then PRESCRIBES
    ! the analytic atm-ocean momentum stress (stress_atmoce_x/y) + the SSS-restoring
    ! climatology (Ssurf), runs the M3e kernels oce_fluxes_mom + oce_fluxes, and dumps
    ! heat_flux/water_flux/virtual_salt/relax_salt (nod2D) + stress_surf (2,elem2D).
    ! tools/run_iceflux_gate_core2.sh compares vs the FESOM2 oracle
    ! (src/fesom_ice_dump.F90::ice_flux_dump_write) for max|delta|=0, CORE2 1-rank, BOTH
    ! whichEVP=0 and 1 (the advection upstream). The 5 fields == the proven M2.11c-2
    ! fesom_flux_dump field set (heat_flux/water_flux/virtual_salt/relax_salt/stress_surf):
    ! M3e produces NATIVELY what the lifecycle has been prescribing.
    !
    ! CONFIG OVERRIDES (the oracle reads these from the CORE2 namelists; the driver matches):
    !   - &ice_dyn / &ice_therm: as M3c/M3d (cd_oce_ice=0.0055, ice_diff=0.0, ice_gamma_fct=0.5,
    !     &ice_therm namelist doubles, albw=0.1, ...).
    !   - atmflux scalars: Ch_atm_ice=Ce_atm_ice=0.00175, ref_sss_local=.true. (=> rsss=S_oc),
    !     use_virt_salt=.true. (which_ALE='linfs'), l_snow=.true., ref_sss=34.0 (dead).
    !   - M3e adds: surf_relax_S=1.929e-06 (namelist.tra double, NOT the o_PARAM default
    !     10/(60*3600*24)); density_0=1030 (mod_constants, == FESOM2 o_PARAM). ocean_area is
    !     geometry (mod_mesh_areas faithful loop == FESOM2 oce_mesh.F90:2385).
    !
    !   FESOM3_MESH_DIR  mesh dir   (default: CORE2)
    !   FESOM3_IC_FILE   IC netcdf  (default: pool phc3.0_winter.nc — same file the oracle reads)
    !   FESOM3_WHICHEVP  0=EVP / 1=mEVP (default 0)
    !   FESOM3_FLUX_OUT  out path   (default: iceflux_f3.bin)
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
    use mod_ice_dyn,        only: ocean2ice
    use mod_ice_step,       only: ice_timestep
    use mod_ice_thermo,     only: t_atmflux
    use mod_ice_oce_coupling, only: oce_fluxes_mom, oce_fluxes
    use oce_initial_state,  only: t_ic3d_config, do_ic3d
    use oce_ssh_rhs,        only: init_stiff_mat_ale
    use mod_advhor_dump,    only: advhor_dump_open, advhor_dump_close, wr_r1, wr_r2
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
    real(kind=WP), allocatable :: stress_surf(:,:)
    integer :: nsw, n, nl, nzmin, nzmax, e, u, whichevp, ios
    real(kind=WP) :: lon, lat
    real(kind=MP) :: zbar_srf, zbar_bot

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2'
    call get_environment_variable('FESOM3_IC_FILE', ic_file)
    if (len_trim(ic_file) == 0) &
        ic_file = '/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc'
    call get_environment_variable('FESOM3_FLUX_OUT', out_path)
    if (len_trim(out_path) == 0) out_path = 'iceflux_f3.bin'
    whichevp = 0
    call get_environment_variable('FESOM3_WHICHEVP', whichevp_str)
    if (len_trim(whichevp_str) > 0) read(whichevp_str, *, iostat=ios) whichevp

    call par_init(partit)
    if (partit%npes /= 1) then
        if (partit%mype == 0) write(*,'(a)') 'fesom_icefluxdump: requires 1 rank'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype); error stop 1
    end if

    !===========================================================================
    ! mesh + geometry (CORE2 rotation 50/15/-90, identical to the lifecycle/IC/ice gate).
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    nl = mesh%nl
    write(*,'(a,i0,a,i0,a,i0,a,i0,a,es16.8)') 'fesom_icefluxdump: nod2D=', mesh%nod2D, &
        ' elem2D=', mesh%elem2D, ' nl=', nl, ' CW swaps=', nsw, ' ocean_area=', mesh%ocean_area

    !===========================================================================
    ! ALE depths Z_3d_n (do_ic3d) + zbar_e_bot (init_stiff_mat_ale). linfs full cells,
    ! eta=0 — identical to fesom_icethermodump.
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
    ! surface T/S ocean2ice copies into srfoce_temp/salt -> the thermo's T_oc/S_oc, AND the
    ! surface S the oce_fluxes virtual_salt / relax_salt read).
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
    write(*,'(a,2es12.4)') 'fesom_icefluxdump: surface S range = ', &
        minval(tracers%data(2)%values(1,:)), maxval(tracers%data(2)%values(1,:))

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
    ! &ice_therm namelist DOUBLES (the M3d override; byte-match the oracle's namelist read).
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

    !===========================================================================
    ! atmflux: allocate (nod2D) + scalar config. The atm scalars match the CORE2 namelists.
    call alloc_atm(atm, mesh%nod2D)
    atm%Ch_atm_ice    = 0.00175_WP      ! namelist.forcing double (NOT 1.75e-3 single->WP)
    atm%Ce_atm_ice    = 0.00175_WP
    atm%ref_sss       = 34.0_WP         ! namelist.tra (dead — ref_sss_local=.true. => rsss=S_oc)
    atm%ref_sss_local = .true.
    atm%use_virt_salt = .true.          ! which_ALE='linfs'
    atm%l_snow        = .true.          ! namelist.forcing
    atm%surf_relax_S  = 1.929e-06_WP    ! namelist.tra double (NOT the o_PARAM 10/(60*3600*24) default)

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
    ! prescribe elevation (hbar), wind-on-ice stress, the atmospheric forcing arrays, AND
    ! (M3e) the atm-ocean momentum stress + the SSS-restoring climatology. ALL formulas MUST
    ! equal the FESOM2 oracle src/fesom_ice_dump.F90::ice_flux_dump_write.
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
        ! M3e atm-ocean stress + SSS-restoring climatology (distinct from the ice forcing).
        atm%stress_atmoce_x(n) =  0.06_WP*cos(lat)*sin(2.0_WP*lon)
        atm%stress_atmoce_y(n) = -0.05_WP*sin(2.0_WP*lat)*cos(lon)
        atm%Ssurf(n)           =  34.0_WP + 0.5_WP*cos(lat) + 0.2_WP*sin(lon)
    end do

    !===========================================================================
    ! ocean -> ice coupling + EVP dynamics (M3b) + ice FCT advection (M3c) + thermodynamics
    ! (M3d) -> advected + thermo'd ice state + the air-sea fluxes ice%flx_h/flx_fw +
    ! atm%evaporation/ice_sublimation (the oce_fluxes inputs).
    call ocean2ice(ice, dyn, tracers, mesh)
    call ice_timestep(ice, mesh, atm)
    write(*,'(a,2es12.4)') 'fesom_icefluxdump: post-thermo max|flx_h|/|flx_fw| = ', &
        maxval(abs(ice%flx_h(1:mesh%nod2D))), maxval(abs(ice%flx_fw(1:mesh%nod2D)))

    !===========================================================================
    ! M3e kernels: oce_fluxes_mom (momentum stress -> stress_surf) + oce_fluxes (heat/water/
    ! virtual_salt/relax_salt). The PAYOFF — the ocean step's surface BCs, native at last.
    allocate(stress_surf(2, mesh%elem2D)); stress_surf = 0.0_WP
    call oce_fluxes_mom(ice, atm, stress_surf, mesh)
    call oce_fluxes(ice, tracers, atm, mesh)
    write(*,'(a,4es12.4)') 'fesom_icefluxdump: max|hf|/|wf|/|vs|/|rs| = ', &
        maxval(abs(atm%heat_flux(1:mesh%nod2D))),  maxval(abs(atm%water_flux(1:mesh%nod2D))), &
        maxval(abs(atm%virtual_salt(1:mesh%nod2D))), maxval(abs(atm%relax_salt(1:mesh%nod2D)))

    !===========================================================================
    ! dump (FADVHDMP, 1-rank identity ids). The 5 proven flux-dump fields: heat_flux/
    ! water_flux/virtual_salt/relax_salt (nod2D) + stress_surf (2,elem2D, wr_r2).
    call advhor_dump_open(u, trim(out_path), mesh%nod2D, mesh%elem2D, mesh%edge2D, nl)
    call wr_r1(u, 'heat_flux',    atm%heat_flux(1:mesh%nod2D))
    call wr_r1(u, 'water_flux',   atm%water_flux(1:mesh%nod2D))
    call wr_r1(u, 'virtual_salt', atm%virtual_salt(1:mesh%nod2D))
    call wr_r1(u, 'relax_salt',   atm%relax_salt(1:mesh%nod2D))
    call wr_r2(u, 'stress_surf',  stress_surf(:, 1:mesh%elem2D))
    call advhor_dump_close(u)
    write(*,'(a)') 'fesom_icefluxdump: wrote '//trim(out_path)

    call par_ex(partit%MPI_COMM_FESOM, partit%mype)

contains

    ! allocate + zero all t_atmflux per-node arrays (sized nod2D for the 1-rank gate),
    ! including the M3e oce_fluxes I/O arrays (heat_flux/water_flux/.../stress_node_surf/
    ! stress_atmoce_x/y/Ssurf).
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

end program fesom_icefluxdump
