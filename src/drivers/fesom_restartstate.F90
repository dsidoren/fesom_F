program fesom_restartstate
    ! Restart Stage 3 (Task 3.4) GATE driver: exercise restart_register_state — the FULL oce+ice
    ! prognostic field set incl. EVP sigma and the real(MP) mesh%hbar/hnode — by writing one complete
    ! checkpoint on a real pi mesh WITHOUT the ocean (fast, MR-capable). The state analog of
    ! fesom_restartsmoke: instead of one synthetic node field it allocates correctly-shaped REAL
    ! t_dyn / t_tracer / t_ice instances + mesh%hbar/hnode, fills the OWNED slots with a
    ! partition-independent synthetic formula, associates them via restart_register_state (the exact
    ! accessors Stage 5 uses: dyn%uv(1,:,:), tracers%data(j)%valuesold(1,:,:), ice%work%sigma11, ...),
    ! and calls restart_write. tools/zarr_diff.py --restart-state then asserts EVERY expected store
    ! exists with the right entity (nod2/elem) x level-kind (nz / nz1 / 2-D) shape, embedded finite
    ! lon/lat, finite data, that element stores are ~2x the node count, and value-checks each store
    ! against its class formula (so the MP->WP staging of hbar/hnode is proven lossless). No FESOM2 oracle.
    !
    ! Synthetic formulas (canonical id g, 1-based; level L 1-based) — MUST match zarr_diff.py:
    !   node 2-D : value = g                         (eta_n d_eta hbar ssh_rhs_old area hice hsnow uice vice t_skin)
    !   node 3-D : value = g + 0.5*L                 (hnode temp* salt* w w_expl w_impl tke)
    !   elem 2-D : value = g                         (sigma11 sigma12 sigma22)
    !   elem 3-D : value = g + 0.5*L                 (u v urhs_AB vrhs_AB [+_AB3])
    !
    !   FESOM3_MESH_DIR          mesh dir (default: pi)
    !   FESOM3_RESTARTSTATE_DIR  output dir (default: ./restartstate)
    !   FESOM3_AB_ORDER          1|2|3 Adams-Bashforth order (default 2; 3 adds *_AB3 + <tr>_M2 stores)
    !   FESOM3_CHUNK_HORIZ       horizontal chunk (default 1000 => multiple chunks + writer subset)
    use mpi
    use, intrinsic :: iso_fortran_env, only: real64
    use mod_precision,    only: WP, MP
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex, set_partition
    use mod_mesh_read,    only: read_mesh
    use mod_mesh_areas,   only: compute_geometry
    use mod_part_bounds,  only: is_multirank, local_dims
    use mod_dyn,          only: t_dyn
    use mod_tracer,       only: t_tracer
    use mod_ice,          only: t_ice
    use mod_io_restart
    implicit none
    character(len=512) :: mesh_dir, out_dir, abbuf
    type(t_partit)         :: partit
    type(t_mesh),   target :: mesh
    type(t_dyn),    target :: dyn
    type(t_tracer), target :: tracers
    type(t_ice),    target :: ice
    type(t_restart) :: R
    logical :: mr
    integer :: nsw, i, j, L, gid, nl, ab, ios
    real(WP) :: vn, ve, v3
    integer :: nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call get_environment_variable('FESOM3_RESTARTSTATE_DIR', out_dir)
    if (len_trim(out_dir) == 0) out_dir = 'restartstate'
    ab = 2
    call get_environment_variable('FESOM3_AB_ORDER', abbuf)
    if (len_trim(abbuf) > 0) then
        read(abbuf, *, iostat=ios) i
        if (ios == 0 .and. (i == 1 .or. i == 2 .or. i == 3)) ab = i
    end if

    call par_init(partit)
    call set_partition(partit, trim(mesh_dir))
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)

    mr = is_multirank(partit)
    call local_dims(mesh, partit, nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF)
    nl = mesh%nl

    call restart_init(R, mesh, partit, chunk_horiz=1000)

    ! ---- allocate correctly-shaped REAL model state (node arrays nNodL, element arrays nElemF) -------
    dyn%AB_order = ab
    allocate(dyn%eta_n(nNodL), dyn%d_eta(nNodL), dyn%ssh_rhs_old(nNodL))
    allocate(dyn%w(nl,nNodL), dyn%w_e(nl,nNodL), dyn%w_i(nl,nNodL))
    allocate(dyn%uv(2, nl-1, nElemF))
    allocate(dyn%uv_rhsAB(ab-1, 2, nl-1, nElemF))
    allocate(dyn%work%tke(nl, nNodL))                      ! exercise the optional TKE path (mix_scheme=5)
    if (.not. allocated(mesh%hbar))  allocate(mesh%hbar(nNodL))        ! real(MP) — precision gotcha
    if (.not. allocated(mesh%hnode)) allocate(mesh%hnode(nl-1, nNodL)) ! real(MP)

    tracers%num_tracers = 2
    allocate(tracers%data(2))
    do j = 1, 2
        tracers%data(j)%ID       = j                       ! 1=temp, 2=salt
        tracers%data(j)%AB_order = ab
        allocate(tracers%data(j)%values  (nl-1, nNodL))
        allocate(tracers%data(j)%valuesAB(nl-1, nNodL))
        allocate(tracers%data(j)%valuesold(max(1,ab-1), nl-1, nNodL))  ! M1 at (1,:,:); M2 at (2,:,:) if ab3
    end do

    ice%num_itracers = 3
    allocate(ice%data(3))
    do j = 1, 3
        allocate(ice%data(j)%values(nNodL))
    end do
    allocate(ice%uice(nNodL), ice%vice(nNodL))
    allocate(ice%thermo%t_skin(nNodL))                     ! carried thermo skin temp (NODE 2-D)
    allocate(ice%work%sigma11(nElemF), ice%work%sigma12(nElemF), ice%work%sigma22(nElemF))

    ! ---- zero-init (halo + unused-by-writer slots stay finite) ----------------------------------------
    dyn%eta_n = 0; dyn%d_eta = 0; dyn%ssh_rhs_old = 0
    dyn%w = 0; dyn%w_e = 0; dyn%w_i = 0; dyn%uv = 0; dyn%uv_rhsAB = 0; dyn%work%tke = 0
    mesh%hbar = 0; mesh%hnode = 0
    do j = 1, 2
        tracers%data(j)%values = 0; tracers%data(j)%valuesAB = 0; tracers%data(j)%valuesold = 0
    end do
    do j = 1, 3
        ice%data(j)%values = 0
    end do
    ice%uice = 0; ice%vice = 0; ice%thermo%t_skin = 0
    ice%work%sigma11 = 0; ice%work%sigma12 = 0; ice%work%sigma22 = 0

    ! ---- fill OWNED nodes (value = g for 2-D ; g + 0.5*L for 3-D) -------------------------------------
    do i = 1, nNodO
        gid = i; if (mr) gid = partit%myList_nod2D(i)
        vn  = real(gid, WP)
        dyn%eta_n(i)          = vn
        dyn%d_eta(i)          = vn
        dyn%ssh_rhs_old(i)    = vn
        mesh%hbar(i)          = real(gid, MP)
        ice%data(1)%values(i) = vn; ice%data(2)%values(i) = vn; ice%data(3)%values(i) = vn
        ice%uice(i)           = vn; ice%vice(i)           = vn
        ice%thermo%t_skin(i)  = vn
        do L = 1, nl-1
            v3 = vn + 0.5_WP*real(L, WP)
            mesh%hnode(L,i)               = real(gid, MP) + 0.5_MP*real(L, MP)
            tracers%data(1)%values(L,i)   = v3; tracers%data(2)%values(L,i)   = v3
            tracers%data(1)%valuesAB(L,i) = v3; tracers%data(2)%valuesAB(L,i) = v3
            tracers%data(1)%valuesold(:,L,i) = v3; tracers%data(2)%valuesold(:,L,i) = v3
        end do
        do L = 1, nl
            v3 = vn + 0.5_WP*real(L, WP)
            dyn%w(L,i) = v3; dyn%w_e(L,i) = v3; dyn%w_i(L,i) = v3; dyn%work%tke(L,i) = v3
        end do
    end do

    ! ---- fill OWNED elements (value = g for 2-D ; g + 0.5*L for 3-D) ----------------------------------
    do i = 1, nElemO
        gid = i; if (mr) gid = partit%myList_elem2D(i)
        ve  = real(gid, WP)
        ice%work%sigma11(i) = ve; ice%work%sigma12(i) = ve; ice%work%sigma22(i) = ve
        do L = 1, nl-1
            v3 = ve + 0.5_WP*real(L, WP)
            dyn%uv(:,L,i)        = v3
            dyn%uv_rhsAB(:,:,L,i) = v3
        end do
    end do

    ! ---- register the full state + write one checkpoint (mix_scheme=5 => tke registered) --------------
    call restart_register_state(R, dyn, tracers, ice, mesh, mix_scheme=5)
    call restart_write(R, trim(out_dir), 2000, 1, 3600.0_real64, globalstep=1)
    call restart_finalize(R)

    if (partit%mype == 0) then
        write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') 'fesom_restartstate: nod2D=', mesh%nod2D, &
            ' elem2D=', mesh%elem2D, ' nl=', nl, ' AB_order=', ab, ' fields=', R%nf
        write(*,'(a)') 'RESTARTSTATE OK'
    end if
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_restartstate
