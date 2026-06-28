program fesom_outputsmoke
    ! M9 Stage 2 (Task 2.1) GATE driver: exercise mod_io_means end-to-end on a real mesh WITHOUT the
    ! ocean (fast, MR-capable) — the field-output analog of fesom_zarrsmoke / fesom_meshdiagdump.
    !
    ! Registers two synthetic node-scalar fields whose value at (record k, CANONICAL node g) follows a
    ! partition-independent generator formula, writes NREC snapshot records, then tools/zarr_diff.py
    ! --output verifies the per-variable-per-year stores: dims (time,nod2), CF time coord, embedded
    ! lon/lat, and every value == the formula. Same formula at any rank count => dist_2 == dist_8
    ! (partition independence) AND == 1-rank, with no FESOM2 oracle needed (writer self-consistency).
    !
    !   FESOM3_MESH_DIR      mesh dir (default: pi)
    !   FESOM3_OUTSMOKE_DIR  output dir (default: ./outsmoke)
    !   FESOM3_OUTSMOKE_NREC number of records (default 5)
    !   FESOM3_CHUNK_HORIZ   horizontal chunk (default 1000 here => multiple chunks + writer subset)
    use mpi
    use, intrinsic :: iso_fortran_env, only: real64
    use mod_precision,    only: WP
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex, set_partition
    use mod_mesh_read,    only: read_mesh
    use mod_mesh_areas,   only: compute_geometry
    use mod_part_bounds,  only: is_multirank, local_dims
    use, intrinsic :: iso_fortran_env, only: real64
    use mod_io_means
    implicit none
    character(len=512) :: mesh_dir, out_dir, buf
    type(t_partit) :: partit
    type(t_mesh)   :: mesh
    type(t_io_means) :: io
    type(t_means_clock) :: clk
    logical :: mr
    integer :: nsw, nrec, k, i, j, L, nl, ios, gid
    integer :: nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF, ge
    real(WP), allocatable :: fld_a(:), fld_b(:), fld_m(:), fld_3(:,:), fld_u(:,:), fld_v(:,:), fld_f2(:)
    ! Task 2.7 ELEMENT fields: fld_e2 (2-D scalar), fld_e3 (3-D scalar, full levels nz), fld_eu/fld_ev
    ! (3-D vector pair, nz1 layers) — exercise the element decomp + elem-centroid r2g + element mask.
    real(WP), allocatable :: fld_e2(:), fld_e3(:,:), fld_eu(:,:), fld_ev(:,:)

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call get_environment_variable('FESOM3_OUTSMOKE_DIR', out_dir)
    if (len_trim(out_dir) == 0) out_dir = 'outsmoke'
    nrec = 5
    call get_environment_variable('FESOM3_OUTSMOKE_NREC', buf, status=ios)
    if (ios == 0 .and. len_trim(buf) > 0) read(buf, *, iostat=ios) nrec

    call par_init(partit)
    call set_partition(partit, trim(mesh_dir))
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)

    mr = is_multirank(partit)
    call local_dims(mesh, partit, nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF)
    nl = mesh%nl

    call means_init(io, trim(out_dir), mesh, partit, chunk_horiz=1000)
    call means_define_node2d(io, 'fld_a', 'synthetic snapshot a', 'm')
    call means_define_node2d(io, 'fld_b', 'synthetic snapshot b', 'psu')
    ! fld_m is a MEAN stream: 3 sub-steps per record with values g-1, g, g+1 => mean == g (tests the
    ! running-sum accumulate + divide-by-count in the output precision, not just an overwrite).
    call means_define_node2d(io, 'fld_m', 'synthetic mean m', '1', mean=.true.)
    ! fld_3 is a 3-D tracer-like field (nl-1 layers): value(g,L) = g+L at valid levels, below-bottom
    ! (L >= nlevels_nod2D) masked to NC_FILL (-> NaN). Tests the 3-D chunk + nlevels mask + vert coord.
    call means_define_node3d(io, 'fld_3', 'synthetic 3d', 'X', on_full_levels=.false.)
    ! fld_u/fld_v: a 3-D node VECTOR pair (nl-1 layers), SNAPSHOT, f8. Exact polynomial values per
    ! (canonical g, layer L) so the native frame is bit-exact AND numpy can reproduce the r2g rotation.
    ! FESOM3_VEC_FRAME (set by the gate) picks native (raw) vs geographic (r2g-rotated to geo). Tests the
    ! vector pairing + per-node rotation + that it stays partition-independent (canonical node coords).
    call means_define_vector3d(io, 'fld_u', 'fld_v', 'synthetic vector u', 'synthetic vector v', 'm/s', &
                               on_full_levels=.false., precision='double')
    ! fld_f2: a 2-D snapshot on a freq=2 STEP cadence (Task 2.6) — written only when mod(istep,2)==0,
    ! so it ends up with floor(nrec/2) records. Value = gid+1000 (k-independent) so the gate checks the
    ! per-field step_event fired (record count) without needing the k->istep mapping.
    call means_define_node2d(io, 'fld_f2', 'synthetic freq-2 snapshot', '1', freq=2, unit='s')
    ! Task 2.7 ELEMENT streams (canonical elem id g): fld_e2 2-D scalar; fld_e3 3-D scalar on FULL levels
    ! (nz, voff=0 — the element nlevels mask); fld_eu/fld_ev 3-D vector pair on nz1 layers, r2g-rotated at
    ! the elem centroid (geographic) or raw (native), f8 so native is bit-exact + numpy reproduces the rot.
    call means_define_elem2d(io, 'fld_e2', 'synthetic element scalar', '1')
    call means_define_elem3d(io, 'fld_e3', 'synthetic element 3d (full levels)', 'X', on_full_levels=.true.)
    call means_define_vector3d_elem(io, 'fld_eu', 'fld_ev', 'synthetic elem vector u', &
                                    'synthetic elem vector v', 'm/s', on_full_levels=.false., precision='double')

    allocate(fld_a(max(1,nNodO)), fld_b(max(1,nNodO)), fld_m(max(1,nNodO)), fld_3(nl-1, max(1,nNodO)), &
             fld_u(nl-1, max(1,nNodO)), fld_v(nl-1, max(1,nNodO)), fld_f2(max(1,nNodO)), &
             fld_e2(max(1,nElemO)), fld_e3(nl, max(1,nElemO)), &
             fld_eu(nl-1, max(1,nElemO)), fld_ev(nl-1, max(1,nElemO)))
    do k = 0, nrec - 1
        ! mean stream: accumulate 3 sub-steps (g-1, g, g+1) -> mean = g
        do j = 1, 3
            do i = 1, nNodO
                gid = i
                if (mr) gid = partit%myList_nod2D(i)
                fld_m(i) = real(gid, WP) + real(j - 2, WP)
            end do
            call means_accumulate(io, 'fld_m', fld_m(1:nNodO))
        end do
        ! snapshot streams: accumulate once (overwrite) -> value == formula
        do i = 1, nNodO
            gid = i
            if (mr) gid = partit%myList_nod2D(i)
            fld_a(i) = real(gid, WP)        + 0.25_WP * real(k, WP)
            fld_b(i) = real(gid, WP)*0.5_WP - real(k, WP)
        end do
        call means_accumulate(io, 'fld_a', fld_a(1:nNodO))
        call means_accumulate(io, 'fld_b', fld_b(1:nNodO))
        ! 3-D snapshot: value(g, L) = g + L (filled at ALL layers; means_write masks below-bottom)
        do i = 1, nNodO
            gid = i
            if (mr) gid = partit%myList_nod2D(i)
            do L = 1, nl-1
                fld_3(L,i) = real(gid + L, WP)
            end do
        end do
        call means_accumulate(io, 'fld_3', fld_3(1:nl-1, 1:nNodO))
        ! vector pair (snapshot): u(g,L)=0.001g+0.5L, v(g,L)=-0.002g+0.25L+1 — exact in f8, bounded O(1-10)
        do i = 1, nNodO
            gid = i
            if (mr) gid = partit%myList_nod2D(i)
            do L = 1, nl-1
                fld_u(L,i) = real(gid,WP)*0.001_WP  + real(L,WP)*0.5_WP
                fld_v(L,i) = real(gid,WP)*(-0.002_WP) + real(L,WP)*0.25_WP + 1.0_WP
            end do
        end do
        call means_accumulate(io, 'fld_u', fld_u(1:nl-1, 1:nNodO))
        call means_accumulate(io, 'fld_v', fld_v(1:nl-1, 1:nNodO))
        ! freq-2 field: snapshot value gid+1000 (constant in k)
        do i = 1, nNodO
            gid = i
            if (mr) gid = partit%myList_nod2D(i)
            fld_f2(i) = real(gid, WP) + 1000.0_WP
        end do
        call means_accumulate(io, 'fld_f2', fld_f2(1:nNodO))
        ! ELEMENT fields (canonical elem id ge): e2=0.1*ge+500; e3(L)=ge+100*L (full levels nz);
        ! eu(L)=0.001*ge+0.5*L, ev(L)=-0.002*ge+0.25*L+1 (same form as the node vector, at elem coords).
        do i = 1, nElemO
            ge = i
            if (mr) ge = partit%myList_elem2D(i)
            fld_e2(i) = real(ge, WP) + 500.0_WP        ! integer-valued => float32-exact (no round-off)
            do L = 1, nl
                fld_e3(L,i) = real(ge, WP) + real(L, WP)*100.0_WP
            end do
            do L = 1, nl-1
                fld_eu(L,i) = real(ge,WP)*0.001_WP    + real(L,WP)*0.5_WP
                fld_ev(L,i) = real(ge,WP)*(-0.002_WP) + real(L,WP)*0.25_WP + 1.0_WP
            end do
        end do
        call means_accumulate(io, 'fld_e2', fld_e2(1:nElemO))
        call means_accumulate(io, 'fld_e3', fld_e3(1:nl,   1:nElemO))
        call means_accumulate(io, 'fld_eu', fld_eu(1:nl-1, 1:nElemO))
        call means_accumulate(io, 'fld_ev', fld_ev(1:nl-1, 1:nElemO))

        ! synthetic clock: year 2000, day 1, timenew = k*3600 -> time_sec = k*3600 (matches the gate).
        ! istep = k+1 drives the per-field STEP events (fld_f2 freq=2 => due at even istep).
        clk%year = 2000; clk%yearstart = 2000; clk%daynew = 1; clk%ndpyr = 365
        clk%month = 1;   clk%day_in_month = 1; clk%ndim_month = 31
        clk%timenew = real(k, real64) * 3600.0_real64
        call means_output(io, k + 1, clk)
    end do
    call means_finalize(io)

    if (partit%mype == 0) then
        write(*,'(a,i0,a,i0,a)') 'fesom_outputsmoke: nod2D=', mesh%nod2D, ' nrec=', nrec, ' OK'
        write(*,'(a)') 'OUTPUTSMOKE OK'
    end if
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_outputsmoke
