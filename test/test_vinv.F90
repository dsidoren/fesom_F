program test_vinv
    ! Vector-invariant momentum advection (momadv_opt==1, src/oce/oce_dyn_vinv.F90).
    !
    ! WRITTEN BEFORE THE OPERATORS (TDD), because there is no oracle: the FESOM2 v2.7.3
    ! momadv_opt==1 branch aborts ("not adapted mom_adv advection typ for ALE"), so these
    ! assertions ARE the specification. See
    ! docs/plans/2026-10-01-momadv-vector-invariant.md.
    !
    ! WHY cartesian=.true.
    ! --------------------
    ! On the SPHERICAL metric the median-dual contour does not close: edge_cross_dxdy(1:2)
    ! carries elem_cos(el1) and (3:4) carries elem_cos(el2) at the SAME edge midpoint
    ! (src/mesh/mod_mesh_areas.F90:336-347), so a uniform field leaves a residual of
    ! relative size (L/R)*tan(lat) ~ 3% on pi -- about 1e9 x round-off. No exactness
    ! assertion is possible there. compute_geometry(..., cartesian=.true.) sets
    ! elem_cos = 1 and metric_factor = 0 (mod_mesh_areas.F90:239-242), the contour closes,
    ! and the discrete circulation integral becomes EXACT for a velocity field linear in
    ! (x,y): for straight segments the quadrature of the closed line integral is exact, so
    ! zeta = circulation/area reproduces curl(u) to round-off. That is what makes V2 below
    ! a real test rather than a tolerance guess.
    !
    ! WHAT IS ASSERTED
    ! ----------------
    !   V1 uniform flow     -> zeta == 0 to round-off. Catches a sign error or a missing
    !                          one-sided level range: a constant field must cancel edge by
    !                          edge.
    !   V2 linear shear     -> u = alpha*y, v = 0  =>  zeta == -alpha to ROUND-OFF.
    !                          THE load-bearing test. Getting -alpha exactly requires
    !                          area(n) == sum over wet adjacent elements of (1/3)*elem_area,
    !                          so this is what pins the area(nz,n) -> area(n) deviation.
    !                          The field depends on LATITUDE ONLY, so it is continuous
    !                          across pi's cyclic seam (a field linear in longitude would
    !                          not be, and would give garbage there).
    !   V3 linearity        -> zeta(-u) == -zeta(u) exactly.
    !   V4 non-degeneracy   -> zeta /= 0 somewhere. Catches a stub or a no-op (the whole
    !                          point of writing this test before the implementation).
    !   V5 finiteness       -> no NaN/Inf anywhere, including at bathymetry steps.
    !
    ! RESTRICTIONS on every assertion: INTERIOR nodes only (the dual contour is open at a
    ! boundary node, so zeta there is not curl(u) -- pi has ~455 boundary edges), and only
    ! levels nz < nlevels_nod2D_min(n) where ALL adjacent elements are wet (at a bathymetry
    ! step the one-sided ranges deliberately drop the dry element's contribution while
    ! area(n) keeps its full weight, so zeta is not curl(u) there either -- that is the
    ! documented zero-with-full-weight rule, not an error).
    use mpi
    use mod_precision,     only: WP, MP
    use mod_constants,     only: r_earth
    use mod_mesh,          only: t_mesh
    use mod_dyn,           only: t_dyn
    use mod_partit,        only: t_partit
    use mod_partitioning,  only: par_init, par_ex, set_partition
    use mod_mesh_read,     only: read_mesh
    use mod_mesh_areas,    only: compute_geometry
    use mod_part_bounds,   only: owned_bounds
    use oce_dyn_vinv,      only: relative_vorticity, momentum_adv_vinv
    implicit none

    character(len=512) :: mesh_dir
    type(t_partit) :: partit
    type(t_mesh)   :: mesh
    type(t_dyn)    :: dyn
    integer :: nfail, nsw, ierr
    integer :: nNodO, nNodL, nEdgeO, nElemO, nElemF, nl, n, nz
    real(kind=WP) :: zmax
    real(kind=WP), allocatable :: zsave(:,:)
    logical,       allocatable :: interior(:), elem_int(:)
    real(kind=WP) :: scal
    integer       :: nzall
    real(kind=WP), parameter   :: ALPHA = 1.0e-5_WP   ! du/dy [1/s]

    nfail = 0
    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call par_init(partit)

    if (partit%npes > 1) call set_partition(partit, trim(mesh_dir))
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    ! cartesian: makes the discrete curl exact for a linear field (see header)
    call compute_geometry(mesh, partit, cartesian=.true.)
    call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
    if (partit%npes == 1) then
        nElemF = mesh%elem2D
    else
        nElemF = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
    end if
    nl = mesh%nl

    allocate(dyn%uv(2, nl-1, nElemF))
    allocate(dyn%work%vorticity(nl-1, nNodL))
    allocate(zsave(nl-1, nNodL))
    allocate(interior(nNodL))
    dyn%uv = 0.0_WP; dyn%work%vorticity = 0.0_WP

    ! interior = touched by no boundary edge (edge_tri(2,ed) <= 0 is the partition-robust
    ! boundary test; the ed > edge2D_in form needs myList_edge2D at npes>1)
    interior = .true.
    do n = 1, nEdgeO                 ! LOCAL owned edges; mesh%edge2D is the GLOBAL count
        if (mesh%edge_tri(2,n) <= 0) then
            if (mesh%edges(1,n) <= nNodL) interior(mesh%edges(1,n)) = .false.
            if (mesh%edges(2,n) <= nNodL) interior(mesh%edges(2,n)) = .false.
        end if
    end do
    write(*,'(a,i7,a,i7)') '  interior nodes: ', count(interior(1:nNodO)), ' of ', nNodO

    !=========================================================================
    ! V1 - uniform flow must give zeta == 0 EXACTLY
    !=========================================================================
    dyn%uv(1,:,:) =  0.17_WP
    dyn%uv(2,:,:) = -0.43_WP
    call relative_vorticity(dyn, mesh, partit)
    zmax = 0.0_WP
    do n = 1, nNodO
        if (.not. interior(n)) cycle
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D_min(n)-1
            zmax = max(zmax, abs(dyn%work%vorticity(nz,n)))
        end do
    end do
    call max_all(zmax)
    write(*,'(a,es12.4)') '  V1 uniform flow   : max|zeta| = ', zmax
    call check_true('V1 uniform flow -> zeta == 0 (interior, all-wet)', zmax < 1.0e-16_WP)

    !=========================================================================
    ! V2 / V4 / V5 - linear shear u = alpha*y : zeta must be EXACTLY -alpha
    !=========================================================================
    call set_shear(1.0_WP)
    call relative_vorticity(dyn, mesh, partit)
    zsave = dyn%work%vorticity

    zmax = owned_absmax(dyn%work%vorticity)
    write(*,'(a,es12.4)') '  V4 non-degeneracy : max|zeta| = ', zmax
    call check_true('V4 zeta /= 0 somewhere (operator is not a no-op)', zmax > 1.0e-12_WP)

    zmax = 0.0_WP
    do n = 1, nNodO
        if (.not. interior(n)) cycle
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D_min(n)-1
            zmax = max(zmax, abs(dyn%work%vorticity(nz,n) + ALPHA))
        end do
    end do
    call max_all(zmax)
    write(*,'(a,es12.4,a,es12.4)') '  V2 linear shear   : max|zeta-(-alpha)| = ', zmax, &
                                   '   relative = ', zmax/ALPHA
    ! DIAGNOSTIC: is the error systematic (a factor) or localised (a few bad nodes)?
    block
        integer :: nok, nbad, nworst
        real(kind=WP) :: r, racc, rworst, latw
        nok = 0; nbad = 0; racc = 0.0_WP; rworst = 0.0_WP; nworst = 0
        do n = 1, nNodO
            if (.not. interior(n)) cycle
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D_min(n)-1
                r = dyn%work%vorticity(nz,n)/(-ALPHA)
                racc = racc + r; nok = nok + 1
                if (abs(r - 1.0_WP) > 0.01_WP) nbad = nbad + 1
                if (abs(r - 1.0_WP) > rworst) then
                    rworst = abs(r - 1.0_WP); nworst = n
                end if
            end do
        end do
        if (nok > 0) then
            latw = real(mesh%coord_nod2D(2, max(nworst,1)), WP)*57.29578_WP
            write(*,'(a,i7,a,f9.4,a,i7,a,f7.2,a,f8.2)') &
                '     levels=', nok, '  mean zeta/(-alpha)=', racc/real(nok,WP), &
                '  off-by->1%=', nbad, '  worst ratio off=', rworst, ' at lat=', latw
        end if
    end block
    call check_true('V2 u=alpha*y -> zeta == -alpha to round-off (interior, all-wet)', &
                    zmax/ALPHA < 1.0e-11_WP)

    call check_true('V5 finite', all(dyn%work%vorticity == dyn%work%vorticity) .and. &
                                 all(abs(dyn%work%vorticity) < huge(1.0_WP)))

    !=========================================================================
    ! V3 - linearity: reversing the flow must flip zeta exactly
    !=========================================================================
    call set_shear(-1.0_WP)
    call relative_vorticity(dyn, mesh, partit)
    zmax = owned_absmax(dyn%work%vorticity + zsave)
    write(*,'(a,es12.4)') '  V3 linearity      : max|zeta(-u)+zeta(u)| = ', zmax
    call check_true('V3 zeta(-u) == -zeta(u)', zmax < 1.0e-16_WP)

    !=========================================================================
    ! V6 - full operator: uniform u with NONZERO w must give ZERO tendency.
    ! This is the Block C net. For a uniform u: zeta == 0 (V1), grad(KE) == 0 at interior
    ! elements, and the vertical flux must telescope EXACTLY -- including across the bottom
    ! face, which is where qq's uvert(nl1+1)=0 breaks under bottom-at-vertices. A `w == 0`
    ! test cannot see that, since every vertical term is proportional to w.
    !=========================================================================
    allocate(mesh%helem(nl-1, nElemF))
    allocate(dyn%uv_rhsAB(1, 2, nl-1, nElemF))
    allocate(dyn%w_e(nl, nNodL))
    allocate(elem_int(nElemO))
    do nz = 1, nl-1
        mesh%helem(nz,:) = 10.0_MP + 2.0_MP*real(nz, MP)   ! any positive ALE thickness
    end do
    do n = 1, nNodL
        do nz = 1, nl
            dyn%w_e(nz,n) = 1.0e-5_WP*real(nl-nz, WP)      ! nonzero AND depth-varying
        end do
    end do
    do n = 1, nElemO
        elem_int(n) = all(interior(mesh%elem2D_nodes(1:3,n)))
    end do

    dyn%uv(1,:,:) =  0.17_WP
    dyn%uv(2,:,:) = -0.43_WP
    dyn%uv_rhsAB  = 0.0_WP
    call momentum_adv_vinv(dyn, mesh, partit)

    ! Restricted to levels where ALL elements adjacent to ALL THREE nodes are wet. Below
    ! that, KE(nz,n) sums only over the elements wet at nz while dividing by the FULL
    ! area(n) -- the zero-with-full-weight rule -- so KE is deliberately not uniform there
    ! and grad(KE) /= 0 even for a uniform u. That is intended, not an error.
    zmax = 0.0_WP; scal = 0.0_WP
    do n = 1, nElemO
        if (.not. elem_int(n)) cycle
        nzall = minval(mesh%nlevels_nod2D_min(mesh%elem2D_nodes(1:3,n))) - 1
        do nz = mesh%ulevels(n), min(mesh%nlevels(n)-1, nzall)
            zmax = max(zmax, abs(dyn%uv_rhsAB(1,1,nz,n)), abs(dyn%uv_rhsAB(1,2,nz,n)))
            ! scale of the individual vertical terms, so the tolerance is relative
            scal = max(scal, abs(dyn%w_e(nz,mesh%elem2D_nodes(1,n))*0.43_WP) &
                             *real(mesh%elem_area(n),WP)/real(mesh%helem(nz,n),WP))
        end do
    end do
    call max_all(zmax); call max_all(scal)
    write(*,'(a,es12.4,a,es12.4)') '  V6 uniform u, w/=0: max|UV_rhsAB| = ', zmax, &
                                   '   relative = ', zmax/max(scal, tiny(1.0_WP))
    call check_true('V6 uniform u + nonzero w -> zero tendency (interior elems, incl bottom face)', &
                    zmax/max(scal, tiny(1.0_WP)) < 1.0e-12_WP)

    !=========================================================================
    ! V7 - the full operator is not a no-op
    !=========================================================================
    call set_shear(1.0_WP)
    dyn%uv_rhsAB = 0.0_WP
    call momentum_adv_vinv(dyn, mesh, partit)
    zmax = 0.0_WP
    do n = 1, nElemO
        do nz = mesh%ulevels(n), mesh%nlevels(n)-1
            zmax = max(zmax, abs(dyn%uv_rhsAB(1,1,nz,n)), abs(dyn%uv_rhsAB(1,2,nz,n)))
        end do
    end do
    call max_all(zmax)
    write(*,'(a,es12.4)') '  V7 full operator  : max|UV_rhsAB| = ', zmax
    call check_true('V7 full operator is not a no-op', zmax > 0.0_WP)

    if (nfail == 0) then
        write(*,'(a)') 'test_vinv: OK'
    else
        write(*,'(a,i0,a)') 'test_vinv: ', nfail, ' FAILURE(S)'
    end if
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)   ! NO third arg: abort is presence-based
    if (nfail /= 0) error stop 1

contains

    subroutine set_shear(sgn)
        ! u = sgn*ALPHA*y, v = 0 with y = r_earth*lat : LINEAR in the cartesian metric, so
        ! the discrete curl is exact and zeta == -sgn*ALPHA. Depends on latitude only, so
        ! it is continuous across pi's cyclic seam. Element value from the mean of its 3
        ! corner latitudes (exact for a linear field).
        !
        ! coord_nod2D, NOT geo_coord_nod2D: the mesh is read with force_rotation=.true., so
        ! the geometry (elem_area, edge_cross_dxdy, element centres) lives in the ROTATED
        ! computational frame. A field linear in GEOGRAPHIC latitude is not linear in that
        ! frame, and the test then fails worst near the rotated pole -- which is exactly how
        ! this was found (mean zeta/(-alpha) = 0.836, worst at rotated lat 84 deg).
        real(kind=WP), intent(in) :: sgn
        integer :: el, k
        real(kind=WP) :: latm
        ! 1..nElemO, NOT nElemF: elem2D_nodes is owned-only, and both triangles of an
        ! owned edge are owned (LESSONS L31), so the edge loop never reads a halo element.
        do el = 1, nElemO
            latm = 0.0_WP
            do k = 1, 3
                latm = latm + real(mesh%coord_nod2D(2, mesh%elem2D_nodes(k,el)), WP)
            end do
            latm = latm/3.0_WP
            dyn%uv(1,:,el) = sgn*ALPHA*r_earth*latm
            dyn%uv(2,:,el) = 0.0_WP
        end do
    end subroutine set_shear

    real(kind=WP) function owned_absmax(f)
        real(kind=WP), intent(in) :: f(:,:)
        integer :: nn, kk
        owned_absmax = 0.0_WP
        do nn = 1, nNodO
            do kk = mesh%ulevels_nod2D(nn), mesh%nlevels_nod2D(nn)-1
                owned_absmax = max(owned_absmax, abs(f(kk,nn)))
            end do
        end do
        call max_all(owned_absmax)
    end function owned_absmax

    subroutine sum_all(x)
        real(kind=WP), intent(inout) :: x
        real(kind=WP) :: t
        if (partit%npes > 1) then
            call MPI_Allreduce(x, t, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                               partit%MPI_COMM_FESOM, ierr)
            x = t
        end if
    end subroutine sum_all

    subroutine max_all(x)
        real(kind=WP), intent(inout) :: x
        real(kind=WP) :: t
        if (partit%npes > 1) then
            call MPI_Allreduce(x, t, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                               partit%MPI_COMM_FESOM, ierr)
            x = t
        end if
    end subroutine max_all

    subroutine sum_int(k)
        integer, intent(inout) :: k
        integer :: t
        if (partit%npes > 1) then
            call MPI_Allreduce(k, t, 1, MPI_INTEGER, MPI_SUM, &
                               partit%MPI_COMM_FESOM, ierr)
            k = t
        end if
    end subroutine sum_int

    subroutine check_true(name, cond)
        character(len=*), intent(in) :: name
        logical,          intent(in) :: cond
        if (.not. cond) then
            nfail = nfail + 1
            if (partit%mype == 0) write(*,'(a)') '  FAIL: '//name
        end if
    end subroutine check_true

end program test_vinv
