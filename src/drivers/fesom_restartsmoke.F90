program fesom_restartsmoke
    ! Restart Stage 3 (Task 3.2) GATE driver: exercise the mod_io_restart WRITE path end-to-end on a
    ! real pi mesh WITHOUT the ocean (fast, MR-capable) — the restart analog of fesom_outputsmoke.
    !
    ! Registers ONE node field whose owned value at CANONICAL node g is the partition-independent
    ! formula value(g) = real(g), associates it as a live pointer, and calls restart_write to produce
    ! a checkpoint folder  <dir>/fesom.2000.001.03600/  holding eta_n.zarr (a single-variable,
    ! single-entity SNAPSHOT store — NO time dim) + checkpoint.json. tools/zarr_diff.py --restart then
    ! asserts the store opens in xarray with dims (nod2,), embedded finite lon/lat, the data var's
    ! _ARRAY_DIMENSIONS, value(g)==g (writer self-consistency; same at any rank count), and that
    ! checkpoint.json parses with the expected keys. No FESOM2 oracle needed.
    !
    !   FESOM3_MESH_DIR          mesh dir (default: pi)
    !   FESOM3_RESTARTSMOKE_DIR  output dir (default: ./restartsmoke)
    !   FESOM3_CHUNK_HORIZ       horizontal chunk (default 1000 here => multiple chunks + writer subset)
    use mpi
    use, intrinsic :: iso_fortran_env, only: real64
    use mod_precision,    only: WP
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex, set_partition
    use mod_mesh_read,    only: read_mesh
    use mod_mesh_areas,   only: compute_geometry
    use mod_part_bounds,  only: is_multirank, local_dims
    use mod_io_decomp,    only: DECOMP_NODE
    use mod_io_restart
    implicit none
    character(len=512) :: mesh_dir, out_dir
    type(t_partit) :: partit
    type(t_mesh)   :: mesh
    type(t_restart) :: R
    real(WP), allocatable, target :: eta(:)
    real(WP), pointer :: peta(:)
    logical :: mr
    integer :: nsw, i, gid
    integer :: nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call get_environment_variable('FESOM3_RESTARTSMOKE_DIR', out_dir)
    if (len_trim(out_dir) == 0) out_dir = 'restartsmoke'

    call par_init(partit)
    call set_partition(partit, trim(mesh_dir))
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)

    mr = is_multirank(partit)
    call local_dims(mesh, partit, nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF)

    call restart_init(R, mesh, partit, chunk_horiz=1000)

    ! synthetic node field: owned value at canonical node g == real(g) (partition-independent).
    allocate(eta(max(1,nNodO)))
    do i = 1, nNodO
        gid = i
        if (mr) gid = partit%myList_nod2D(i)
        eta(i) = real(gid, WP)
    end do
    peta => eta
    call restart_register_field(R, 'eta_n', 'm', DECOMP_NODE, p2d=peta)

    ! write one checkpoint: year 2000, day-of-year 1, sec-of-day 3600 -> fesom.2000.001.03600/
    call restart_write(R, trim(out_dir), 2000, 1, 3600.0_real64, globalstep=1)
    call restart_finalize(R)

    if (partit%mype == 0) then
        write(*,'(a,i0,a)') 'fesom_restartsmoke: nod2D=', mesh%nod2D, ' OK'
        write(*,'(a)') 'RESTARTSMOKE OK'
    end if
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_restartsmoke
