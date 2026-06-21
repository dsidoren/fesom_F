program fesom_icdump
    ! M2.11b initial-conditions (do_ic3d) byte-gate driver. Loads the CORE2 mesh
    ! 1-rank with the SAME rotation as the FESOM2 CORE2 run (alpha/beta/gamma =
    ! 50/15/-90, cyclic 360, force_rotation; identical to pi — work_core/namelist.config),
    ! builds the ALE depths Z_3d_n/zbar_3d_n (linfs full cells = reference mid-depth at
    ! init, M1.2/pressure-gate proven), allocates a 2-tracer state (data(1)%ID=1 temp,
    ! data(2)%ID=2 salt), and runs the transcribed oce_initial_state::do_ic3d on the
    ! phc3.0_winter.nc climatology (read+bilinear+vertical-interp+extrapolate+insitu2pot).
    ! Dumps Z_3d_n (input check), ic_temp (=data(1), potential T) and ic_salt (=data(2))
    ! in the oracle's FADVHDMP format. tools/pressure_diff.py compares vs the FESOM2
    ! oracle's live Tclim/Sclim (src/fesom_ic_dump.F90) for max|delta|=0.
    !
    !   FESOM3_MESH_DIR  mesh dir   (default: CORE2)
    !   FESOM3_IC_FILE   IC netcdf  (default: pool phc3.0_winter.nc — same file the oracle reads)
    !   FESOM3_IC_OUT    out path   (default: ic_f3.bin)
    use mpi
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use mod_precision,      only: WP, MP
    use mod_mesh,           only: t_mesh
    use mod_partit,         only: t_partit
    use mod_partitioning,   only: par_init, par_ex
    use mod_mesh_read,      only: read_mesh
    use mod_mesh_areas,     only: compute_geometry
    use mod_tracer,         only: t_tracer
    use oce_initial_state,  only: t_ic3d_config, do_ic3d
    use mod_advhor_dump,    only: advhor_dump_open, advhor_dump_close, wr_r2
    implicit none

    character(len=512) :: mesh_dir, out_path, ic_file
    type(t_partit)     :: partit
    type(t_mesh)       :: mesh
    type(t_tracer)     :: tracers
    type(t_ic3d_config):: ic
    integer :: nsw, n, nz, nl, nzmin, nzmax, u
    real(kind=MP) :: zbar_srf, zbar_bot

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2'
    call get_environment_variable('FESOM3_IC_FILE', ic_file)
    if (len_trim(ic_file) == 0) &
        ic_file = '/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc'
    call get_environment_variable('FESOM3_IC_OUT', out_path)
    if (len_trim(out_path) == 0) out_path = 'ic_f3.bin'

    call par_init(partit)
    if (partit%npes /= 1) then
        if (partit%mype == 0) write(*,'(a)') 'fesom_icdump: requires 1 rank'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype); error stop 1
    end if

    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    nl = mesh%nl
    write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') 'fesom_icdump: nod2D=', mesh%nod2D, &
        ' elem2D=', mesh%elem2D, ' edge2D=', mesh%edge2D, ' nl=', nl, ' CW swaps=', nsw

    ! zbar_3d_n / Z_3d_n: per-node ALE interface/mid depths, built exactly as FESOM2
    ! init_ale (oce_ale.F90:531-566) at the initial state (no cavity, full cells, eta=0).
    ! Identical to the pressure-gate build; getcoeffld interpolates onto Z_3d_n and
    ! insitu2pot reads mesh%Z. Dumped as an input so any divergence is localised.
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

    ! 2-tracer ocean state: data(1)=temperature (ID 1), data(2)=salinity (ID 2),
    ! matching work_core/namelist.tra (nml_tracer_list 1,2). values init irrelevant
    ! (do_ic3d overwrites with dummy then the interpolation).
    tracers%num_tracers = 2
    allocate(tracers%data(2))
    allocate(tracers%data(1)%values(nl-1, mesh%nod2D), tracers%data(2)%values(nl-1, mesh%nod2D))
    tracers%data(1)%values = 0.0_WP;  tracers%data(1)%ID = 1
    tracers%data(2)%values = 0.0_WP;  tracers%data(2)%ID = 2

    ! IC config = work_core/namelist.tra &tracer_init3d: n_ic3d=2, idlist=2,1,
    ! filelist=2x phc3.0_winter.nc, varlist='salt','temp', t_insitu=.true.
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

    write(*,'(a,2(es12.4))') 'fesom_icdump: T range = ', &
        minval(tracers%data(1)%values), maxval(tracers%data(1)%values)
    write(*,'(a,2(es12.4))') 'fesom_icdump: S range = ', &
        minval(tracers%data(2)%values), maxval(tracers%data(2)%values)

    call advhor_dump_open(u, trim(out_path), mesh%nod2D, mesh%elem2D, mesh%edge2D, nl)
    call wr_r2(u, 'Z_3d_n',  mesh%Z_3d_n(1:nl-1, 1:mesh%nod2D))
    call wr_r2(u, 'ic_temp', tracers%data(1)%values(1:nl-1, 1:mesh%nod2D))
    call wr_r2(u, 'ic_salt', tracers%data(2)%values(1:nl-1, 1:mesh%nod2D))
    call advhor_dump_close(u)
    write(*,'(a)') 'fesom_icdump: wrote '//trim(out_path)

    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_icdump
