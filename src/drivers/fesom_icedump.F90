program fesom_icedump
    ! M3a sea-ice foundation byte-gate driver. Mirrors fesom_icdump (M2.11b) /
    ! fesom_lifecycle (M2.11c) init up to the ssh_stiff CSR, then runs the transcribed
    ! ice_setup (allocate + ice_mass_matrix_fill + cold-start ice_initial_state) on the
    ! do_ic3d ocean state, and dumps the ice cold-start IC (a_ice/m_ice/m_snow) plus the
    ! FCT mass matrix (fct_massmatrix) in the oracle's FADVHDMP format. tools/ice_diff.py
    ! compares vs the FESOM2 oracle's live ice%data/fct_massmatrix (src/fesom_ice_dump.F90)
    ! for max|delta|=0, CORE2 1-rank.
    !
    !   FESOM3_MESH_DIR  mesh dir   (default: CORE2)
    !   FESOM3_IC_FILE   IC netcdf  (default: pool phc3.0_winter.nc — same file the oracle reads)
    !   FESOM3_ICE_OUT   out path   (default: ice_f3.bin)
    use mpi
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use mod_precision,      only: WP, MP
    use mod_param_phys,     only: alpha, theta
    use mod_mesh,           only: t_mesh
    use mod_partit,         only: t_partit
    use mod_partitioning,   only: par_init, par_ex
    use mod_mesh_read,      only: read_mesh
    use mod_mesh_areas,     only: compute_geometry
    use mod_tracer,         only: t_tracer
    use mod_ice,            only: t_ice
    use mod_ice_setup,      only: ice_setup
    use oce_initial_state,  only: t_ic3d_config, do_ic3d
    use oce_ssh_rhs,        only: init_stiff_mat_ale
    use mod_advhor_dump,    only: advhor_dump_open, advhor_dump_close, wr_r1
    implicit none

    ! CORE2 namelist timestep: dt = 86400/48 = 1800 s (computed like FESOM2, not a literal).
    real(kind=WP), parameter :: dt = 86400.0_WP / real(48, WP)

    character(len=512) :: mesh_dir, out_path, ic_file
    type(t_partit)     :: partit
    type(t_mesh)       :: mesh
    type(t_tracer)     :: tracers
    type(t_ice)        :: ice
    type(t_ic3d_config):: ic
    integer :: nsw, n, nz, nl, nzmin, nzmax, e, u
    real(kind=MP) :: zbar_srf, zbar_bot

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2'
    call get_environment_variable('FESOM3_IC_FILE', ic_file)
    if (len_trim(ic_file) == 0) &
        ic_file = '/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc'
    call get_environment_variable('FESOM3_ICE_OUT', out_path)
    if (len_trim(out_path) == 0) out_path = 'ice_f3.bin'

    call par_init(partit)
    if (partit%npes /= 1) then
        if (partit%mype == 0) write(*,'(a)') 'fesom_icedump: requires 1 rank'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype); error stop 1
    end if

    !===========================================================================
    ! mesh + geometry (CORE2 rotation 50/15/-90, identical to the lifecycle/IC gate).
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    nl = mesh%nl
    write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') 'fesom_icedump: nod2D=', mesh%nod2D, &
        ' elem2D=', mesh%elem2D, ' edge2D=', mesh%edge2D, ' nl=', nl, ' CW swaps=', nsw

    !===========================================================================
    ! ALE depths Z_3d_n (do_ic3d) + zbar_e_bot (init_stiff_mat_ale). linfs full cells,
    ! eta=0 — identical to fesom_lifecycle.
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

    !===========================================================================
    ! 2-tracer ocean state + do_ic3d phc3.0 IC (gives the surface T the cold start reads).
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
    write(*,'(a,2es12.4)') 'fesom_icedump: surface T range = ', &
        minval(tracers%data(1)%values(1,:)), maxval(tracers%data(1)%values(1,:))

    !===========================================================================
    ! ssh_stiff CSR (reduced-M2: alpha=theta=1) — the ice mass matrix rides its sparsity.
    alpha = 1.0_WP; theta = 1.0_WP
    call init_stiff_mat_ale(mesh, dt)

    !===========================================================================
    ! ice setup: allocate + ice_mass_matrix_fill + cold-start ice_initial_state.
    call ice_setup(ice, tracers, mesh, dt)
    write(*,'(a,3es12.4)') 'fesom_icedump: max a_ice/m_ice/m_snow = ', &
        maxval(ice%data(1)%values), maxval(ice%data(2)%values), maxval(ice%data(3)%values)
    write(*,'(a,i0,a,es12.4)') 'fesom_icedump: nza=', mesh%ssh_stiff%nza, &
        '  sum(massmatrix)=', sum(ice%work%fct_massmatrix)

    !===========================================================================
    ! dump (FADVHDMP, 1-rank identity ids). a_ice/m_ice/m_snow on nod2D; massmatrix on nza.
    call advhor_dump_open(u, trim(out_path), mesh%nod2D, mesh%elem2D, mesh%edge2D, nl)
    call wr_r1(u, 'ice_a_ice',      ice%data(1)%values(1:mesh%nod2D))
    call wr_r1(u, 'ice_m_ice',      ice%data(2)%values(1:mesh%nod2D))
    call wr_r1(u, 'ice_m_snow',     ice%data(3)%values(1:mesh%nod2D))
    call wr_r1(u, 'ice_massmatrix', ice%work%fct_massmatrix(1:mesh%ssh_stiff%nza))
    call advhor_dump_close(u)
    write(*,'(a)') 'fesom_icedump: wrote '//trim(out_path)

    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_icedump
