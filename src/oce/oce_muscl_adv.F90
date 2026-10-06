module oce_muscl_adv
    ! MUSCL-type advection setup, transcribed from FESOM2 v2.7.3 oce_muscl_adv.F90.
    !   muscl_adv_init            -> nboundary_lay (per-node bottom-boundary layer)
    !   find_up_downwind_triangles-> edge_up_dn_tri (upwind/downwind tri per edge)
    !   fill_up_dn_grad           -> per-edge up/dn gradient (ORACLE only, see the routine)
    !   muscl_node_grad           -> Miura node-averaged gradient (once per node)
    !
    ! Reference: Abalakin, Dervieux, Kozubskaya (2002), INRIA RR-4459; the concept
    ! of upwind/downwind triangles to a given edge (sergey.danilov@awi.de 2012).
    !
    ! SCOPE / clean-architecture deviations (D7, plan "USE-globals -> explicit args"):
    !  * 1-rank only (myDim_* == global; eDim==eXDim==0). The FESOM2 halo exchanges
    !    (exchange_elem of coord_elem/e_nodes in find_up_downwind_triangles) are
    !    no-ops at 1-rank and are dropped here; multi-rank is deferred to M1.5/M2.12.
    !    coord_elem(:,n,el)==coord_nod2D(:,elem2D_nodes(n,el)) and the FESOM2 global
    !    id myList_nod2D(node)==node, so the e_nodes(n,el)==myList_nod2D(...) tests
    !    become elem2D_nodes(n,el)==node.
    !  * muscl_adv_init's nn_num/nn_pos block (FESOM2 oce_muscl_adv.F90:59-90,106-127)
    !    is OMITTED: those node-neighbour arrays size off SSH_stiff (not built until
    !    M2.6) and are used only by quadratic reconstruction, NOT by horizontal
    !    advection. Only the nboundary_lay loop is kept.
    !  * tr_xy is an explicit argument to fill_up_dn_grad (FESOM2 reads global o_ARRAYS).
    use mod_precision,   only: WP
    use mod_mesh,        only: t_mesh
    use mod_partit,      only: t_partit
    use mod_tracer,      only: t_tracer_work
    use mod_mesh_rotate, only: get_cyclic_length
    use mod_part_bounds, only: owned_bounds, is_multirank
    use mod_halo,        only: exchange_elem_full
    implicit none
    private
    public :: muscl_adv_init, find_up_downwind_triangles, fill_up_dn_grad, muscl_node_grad
    public :: muscl_node_ranges

contains

    !---------------------------------------------------------------------------
    subroutine muscl_adv_init(twork, mesh, partit)
        ! oce_muscl_adv.F90:32-158 (nboundary_lay part only — see module header).
        ! M2.12b: optional partit. nboundary_lay is sized OWNED+HALO and built from
        ! OWNED edges with NO exchange — exactly as FESOM2 (muscl_adv_init:92-155); the
        ! halo entries are the partition-local min, identical on both codes for the same
        ! partition, so the owned-edge MUSCL flux byte-matches.
        type(t_mesh),        intent(in)    :: mesh
        type(t_tracer_work), intent(inout) :: twork
        type(t_partit),      intent(in), optional :: partit
        integer :: n
        integer :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        ! find upwind and downwind triangle for each local edge
        call find_up_downwind_triangles(twork, mesh, partit)

        ! node n becomes a boundary node below layer twork%nboundary_lay(n)
        if (.not. allocated(twork%nboundary_lay)) allocate(twork%nboundary_lay(nNodL))
        twork%nboundary_lay = mesh%nl - 1
        do n = 1, nEdgeO
            if (any(mesh%edge_tri(:, n) <= 0)) then
                ! edge nodes already at the surface boundary: sign(1, nboundary_lay-nz)
                ! must be negative at nz=1, so set nboundary_lay(edge nodes)=0.
                twork%nboundary_lay(mesh%edges(:, n)) = 0
            else
                ! edge becomes a boundary edge with depth (bottom topography): at depth
                ! nboundary_lay the edge still has two valid ocean triangles; below it
                ! the edge is a boundary edge.
                twork%nboundary_lay(mesh%edges(1, n)) = min(twork%nboundary_lay(mesh%edges(1, n)), &
                    minval(mesh%nlevels(mesh%edge_tri(:, n))) - 1)
                twork%nboundary_lay(mesh%edges(2, n)) = min(twork%nboundary_lay(mesh%edges(2, n)), &
                    minval(mesh%nlevels(mesh%edge_tri(:, n))) - 1)
            end if
        end do

        ! levels on which the MUSCL kernels read the node-averaged gradient (static)
        call muscl_node_ranges(twork, mesh, partit)
    end subroutine muscl_adv_init

    !---------------------------------------------------------------------------
    subroutine muscl_node_ranges(twork, mesh, partit)
        ! For every node, the levels on which an OWNED edge reads the Miura node average
        ! gnod: an edge with both up/downwind triangles reads the triangles' gradient on
        ! its shared range [nzmin_e, nzmax_e) (nzmin_e = maxval(ulevels_nod2D_max),
        ! nzmax_e = minval(nlevels_nod2D_min) of its nodes) and gnod on the node's other
        ! wet levels; an edge without them reads gnod on all of them. So node n needs gnod on
        !     [ulevels_nod2D(n), gnod_lo(n))  and  [gnod_hi(n), nlevels_nod2D(n))
        ! with gnod_lo = max_e nzmin_e, gnod_hi = min_e nzmax_e over its owned edges (the
        ! whole wet range if one of them lacks a triangle; nothing if none is owned).
        ! nzmin_e >= ulevels_nod2D(n) and nzmax_e <= nlevels_nod2D(n), so the two ranges
        ! lie in the node's wet range. muscl_node_grad computes only these levels -- on a
        ! flat-bottomed interior the shared ranges cover almost the whole column and the
        ! node average is needed on a few levels only. Levels are static: call once (and
        ! again only if the level structure is changed, as test_muscl_onthefly does).
        type(t_tracer_work), intent(inout) :: twork
        type(t_mesh),        intent(in)    :: mesh
        type(t_partit),      intent(in), optional :: partit
        integer :: edge, k, n, nzmin, nzmax, ednodes(2)
        logical :: both
        integer :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        if (.not. allocated(twork%gnod_lo)) allocate(twork%gnod_lo(nNodL), twork%gnod_hi(nNodL))
        ! no owned edge: nothing needed ([ulev, ulev) and [nlev, nlev) are empty)
        twork%gnod_lo = mesh%ulevels_nod2D(1:nNodL)
        twork%gnod_hi = mesh%nlevels_nod2D(1:nNodL)
        do edge = 1, nEdgeO
            ednodes = mesh%edges(:, edge)
            both  = (twork%edge_up_dn_tri(1, edge) /= 0) .and. (twork%edge_up_dn_tri(2, edge) /= 0)
            nzmin = maxval(mesh%ulevels_nod2D_max(ednodes))
            nzmax = minval(mesh%nlevels_nod2D_min(ednodes))
            do k = 1, 2
                n = ednodes(k)
                if (both) then
                    twork%gnod_lo(n) = max(twork%gnod_lo(n), nzmin)
                    twork%gnod_hi(n) = min(twork%gnod_hi(n), nzmax)
                else
                    twork%gnod_lo(n) = mesh%nlevels_nod2D(n)          ! the whole wet range
                end if
            end do
        end do
    end subroutine muscl_node_ranges

    !---------------------------------------------------------------------------
    subroutine find_up_downwind_triangles(twork, mesh, partit)
        ! oce_muscl_adv.F90:162-352. For each edge, find the triangle around node 1
        ! that contains the direction -edge (upwind), and around node 2 the triangle
        ! containing +edge (downwind). Decompose b (a triangle side) and the search
        ! direction x along c (the other side) and the 90deg-CCW normal, then compare
        ! the atan2 angles.
        !
        ! M2.12b (optional partit): a halo node's element list (nod_in_elem2D) reaches
        ! eXDim elements whose elem2D_nodes are NOT stored locally and whose nodes may
        ! be beyond the local node list. So at npes>1 we build, exactly as FESOM2
        ! (oce_muscl_adv.F90:192-235), coord_elem(2,3,*) and e_nodes(3,*) over the FULL
        ! element halo via exchange_elem_full, and pick the triangle sides from those +
        ! the GLOBAL anchor id (myList_nod2D). At 1-rank (partit absent) the proven path
        ! through tri_sides(elem2D_nodes/coord_nod2D) is used unchanged.
        type(t_mesh),        intent(in)    :: mesh
        type(t_tracer_work), intent(inout) :: twork
        type(t_partit),      intent(in), optional :: partit
        integer       :: n, k, ednodes(2), elem
        real(kind=WP) :: x(2), b(2), c(2), cr, bx, by, xx, xy, ab, ax, cl
        integer       :: nNodO, nNodL, nEdgeO, nElemO, nElemF, kk, nn
        logical       :: lmr
        real(kind=WP), allocatable :: coord_elem(:,:,:), tmp(:)
        integer,       allocatable :: e_nodes(:,:), tmp_i(:)

        lmr = is_multirank(partit)
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        if (.not. allocated(twork%edge_up_dn_tri))  allocate(twork%edge_up_dn_tri(2, nEdgeO))
        twork%edge_up_dn_tri = 0
        cl = get_cyclic_length()

        if (lmr) then
            ! element vertex coords + GLOBAL node ids over the full element halo
            nElemF = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
            allocate(coord_elem(2, 3, nElemF), e_nodes(3, nElemF))
            allocate(tmp(nElemF), tmp_i(nElemF))
            do nn = 1, 3
                do kk = 1, 2
                    tmp = 0.0_WP
                    do elem = 1, nElemO
                        tmp(elem) = mesh%coord_nod2D(kk, mesh%elem2D_nodes(nn, elem))
                    end do
                    call exchange_elem_full(tmp, partit)
                    coord_elem(kk, nn, :) = tmp
                end do
                tmp_i = 0
                do elem = 1, nElemO
                    tmp_i(elem) = partit%myList_nod2D(mesh%elem2D_nodes(nn, elem))
                end do
                call exchange_elem_full(tmp_i, partit)
                e_nodes(nn, :) = tmp_i
            end do
            deallocate(tmp, tmp_i)
        end if

        do n = 1, nEdgeO
            ednodes = mesh%edges(:, n)
            x = mesh%coord_nod2D(:, ednodes(2)) - mesh%coord_nod2D(:, ednodes(1))
            if (x(1) >  cl/2.0_WP) x(1) = x(1) - cl
            if (x(1) < -cl/2.0_WP) x(1) = x(1) + cl

            ! Find upwind (in the sense of x) triangle, i.e. which contains -x:
            x = -x
            do k = 1, mesh%nod_in_elem2D_num(ednodes(1))
                elem = mesh%nod_in_elem2D(k, ednodes(1))
                call get_bc(elem, ednodes(1), b, c)
                if (b(1) >  cl/2.0_WP) b(1) = b(1) - cl
                if (b(1) < -cl/2.0_WP) b(1) = b(1) + cl
                if (c(1) >  cl/2.0_WP) c(1) = c(1) - cl
                if (c(1) < -cl/2.0_WP) c(1) = c(1) + cl
                ! Decompose b and x into parts along c and along (-cy,cx) (90deg CCW)
                cr = sum(c*c)
                bx = sum(b*c)/cr
                by = (-b(1)*c(2)+b(2)*c(1))/cr
                xx = sum(x*c)/cr
                xy = (-x(1)*c(2)+x(2)*c(1))/cr
                ab = atan2(by, bx)
                ax = atan2(xy, xx)
                ! Since b,c are triangle sides, |ab|<pi, so atan2 is what is needed
                if ((ab > 0.0_WP) .and. (ax > 0.0_WP) .and. (ax < ab)) then
                    twork%edge_up_dn_tri(1, n) = elem; cycle
                end if
                if ((ab < 0.0_WP) .and. (ax < 0.0_WP) .and. (ax > ab)) then
                    twork%edge_up_dn_tri(1, n) = elem; cycle
                end if
                if ((ab == ax) .or. (ax == 0.0_WP)) then
                    twork%edge_up_dn_tri(1, n) = elem; cycle
                end if
            end do

            ! Find downwind element (direction +x)
            x = -x
            do k = 1, mesh%nod_in_elem2D_num(ednodes(2))
                elem = mesh%nod_in_elem2D(k, ednodes(2))
                call get_bc(elem, ednodes(2), b, c)
                if (b(1) >  cl/2.0_WP) b(1) = b(1) - cl
                if (b(1) < -cl/2.0_WP) b(1) = b(1) + cl
                if (c(1) >  cl/2.0_WP) c(1) = c(1) - cl
                if (c(1) < -cl/2.0_WP) c(1) = c(1) + cl
                cr = sum(c*c)
                bx = sum(b*c)/cr
                by = (-b(1)*c(2)+b(2)*c(1))/cr
                xx = sum(x*c)/cr
                xy = (-x(1)*c(2)+x(2)*c(1))/cr
                ab = atan2(by, bx)
                ax = atan2(xy, xx)
                if ((ab > 0.0_WP) .and. (ax > 0.0_WP) .and. (ax < ab)) then
                    twork%edge_up_dn_tri(2, n) = elem; cycle
                end if
                if ((ab < 0.0_WP) .and. (ax < 0.0_WP) .and. (ax > ab)) then
                    twork%edge_up_dn_tri(2, n) = elem; cycle
                end if
                if ((ab == ax) .or. (ax == 0.0_WP)) then
                    twork%edge_up_dn_tri(2, n) = elem; cycle
                end if
            end do
        end do

        if (lmr) deallocate(coord_elem, e_nodes)

        ! For edges touching the boundary, up/downwind elements may be absent; the
        ! MUSCL kernels then return to standard Miura (the node average gnod).

    contains
        subroutine get_bc(el, anchor_loc, bb, cc)
            ! Triangle sides from vertex `anchor_loc`. 1-rank: via elem2D_nodes/
            ! coord_nod2D (tri_sides). Multi-rank: via coord_elem + the GLOBAL anchor
            ! id matched against e_nodes (FESOM2 oce_muscl_adv.F90:251-303).
            integer,       intent(in)  :: el, anchor_loc
            real(kind=WP), intent(out) :: bb(2), cc(2)
            integer :: anchor_g
            if (.not. lmr) then
                call tri_sides(mesh, el, anchor_loc, bb, cc)
                return
            end if
            anchor_g = partit%myList_nod2D(anchor_loc)
            if (e_nodes(1, el) == anchor_g) then
                bb = coord_elem(:, 2, el) - coord_elem(:, 1, el)
                cc = coord_elem(:, 3, el) - coord_elem(:, 1, el)
            elseif (e_nodes(2, el) == anchor_g) then
                bb = coord_elem(:, 1, el) - coord_elem(:, 2, el)
                cc = coord_elem(:, 3, el) - coord_elem(:, 2, el)
            else
                bb = coord_elem(:, 1, el) - coord_elem(:, 3, el)
                cc = coord_elem(:, 2, el) - coord_elem(:, 3, el)
            end if
        end subroutine get_bc
    end subroutine find_up_downwind_triangles

    pure subroutine tri_sides(mesh, elem, anchor, b, c)
        ! The two triangle sides emanating from the vertex 'anchor' of 'elem', in
        ! FESOM2's b/c convention (oce_muscl_adv.F90:251-260): if anchor is vertex 1
        ! then b=v2-v1, c=v3-v1; if vertex 2 then b=v1-v2, c=v3-v2; else b=v1-v3,
        ! c=v2-v3. (1-rank: coord_elem == coord_nod2D(:,elem2D_nodes(:,elem)).)
        type(t_mesh),  intent(in)  :: mesh
        integer,       intent(in)  :: elem, anchor
        real(kind=WP), intent(out) :: b(2), c(2)
        integer :: e1, e2, e3
        e1 = mesh%elem2D_nodes(1, elem)
        e2 = mesh%elem2D_nodes(2, elem)
        e3 = mesh%elem2D_nodes(3, elem)
        if (e1 == anchor) then
            b = mesh%coord_nod2D(:, e2) - mesh%coord_nod2D(:, e1)
            c = mesh%coord_nod2D(:, e3) - mesh%coord_nod2D(:, e1)
        elseif (e2 == anchor) then
            b = mesh%coord_nod2D(:, e1) - mesh%coord_nod2D(:, e2)
            c = mesh%coord_nod2D(:, e3) - mesh%coord_nod2D(:, e2)
        else
            b = mesh%coord_nod2D(:, e1) - mesh%coord_nod2D(:, e3)
            c = mesh%coord_nod2D(:, e2) - mesh%coord_nod2D(:, e3)
        end if
    end subroutine tri_sides

    !---------------------------------------------------------------------------
    subroutine fill_up_dn_grad(eudg, twork, tr_xy, mesh, partit)
        ! oce_muscl_adv.F90:356-525. Per edge, build the up/downwind elemental tracer
        ! gradient eudg(1:4,nz,edge): (1,3)=upwind (x,y), (2,4)=downwind.
        ! NOT USED BY THE MODEL: FESOM2 stores this array in twork%edge_up_dn_grad; FESOM3
        ! looks the same values up on the fly in the MUSCL kernels (muscl_node_grad +
        ! the kernels' per-level lookup). Kept as the ORACLE of that equivalence (test_muscl_onthefly)
        ! and for the legacy FESOM2 dump drivers (fesom_advhordump*). eudg is zeroed here,
        ! so levels the fill does not write are 0, as in FESOM2's once-zeroed array.
        ! On shared levels take the gradient straight from edge_up_dn_tri; on
        ! not-shared levels (and on boundary edges) area-weighted-average tr_xy over
        ! the triangles around each edge node (standard Miura).
        ! M2.12b: optional partit -> loop OWNED edges. The body reads nod_in_elem2D
        ! (completed for halo nodes by the find_neighbors dance), tr_xy and elem_area
        ! at halo elements (both halo-exchanged by the caller / compute_geometry).
        type(t_mesh),        intent(in)    :: mesh
        real(kind=WP),       intent(out)   :: eudg(4, mesh%nl-1, *)   ! (4, nl-1, nEdgeO)
        type(t_tracer_work), intent(in)    :: twork
        real(kind=WP),       intent(in)    :: tr_xy(2, mesh%nl-1, mesh%elem2D)
        type(t_partit),      intent(in), optional :: partit
        integer       :: edge, nz, elem, k, ednodes(2), nzmin, nzmax
        real(kind=WP) :: tvol, tx, ty
        integer       :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        eudg(:, :, 1:nEdgeO) = 0.0_WP
        do edge = 1, nEdgeO
            ednodes = mesh%edges(:, edge)
            !___ edge has both upwind and downwind triangle on the surface __________
            if ((twork%edge_up_dn_tri(1, edge) /= 0) .and. (twork%edge_up_dn_tri(2, edge) /= 0)) then
                nzmin = maxval(mesh%ulevels_nod2D_max(ednodes))
                nzmax = minval(mesh%nlevels_nod2D_min(ednodes))

                ! not-shared upper levels of edge node 1
                do nz = mesh%ulevels_nod2D(ednodes(1)), nzmin-1
                    tvol = 0.0_WP; tx = 0.0_WP; ty = 0.0_WP
                    do k = 1, mesh%nod_in_elem2D_num(ednodes(1))
                        elem = mesh%nod_in_elem2D(k, ednodes(1))
                        if (mesh%nlevels(elem)-1 < nz .or. nz < mesh%ulevels(elem)) cycle
                        tvol = tvol + mesh%elem_area(elem)
                        tx = tx + tr_xy(1, nz, elem)*mesh%elem_area(elem)
                        ty = ty + tr_xy(2, nz, elem)*mesh%elem_area(elem)
                    end do
                    eudg(1, nz, edge) = tx/tvol
                    eudg(3, nz, edge) = ty/tvol
                end do
                ! not-shared upper levels of edge node 2
                do nz = mesh%ulevels_nod2D(ednodes(2)), nzmin-1
                    tvol = 0.0_WP; tx = 0.0_WP; ty = 0.0_WP
                    do k = 1, mesh%nod_in_elem2D_num(ednodes(2))
                        elem = mesh%nod_in_elem2D(k, ednodes(2))
                        if (mesh%nlevels(elem)-1 < nz .or. nz < mesh%ulevels(elem)) cycle
                        tvol = tvol + mesh%elem_area(elem)
                        tx = tx + tr_xy(1, nz, elem)*mesh%elem_area(elem)
                        ty = ty + tr_xy(2, nz, elem)*mesh%elem_area(elem)
                    end do
                    eudg(2, nz, edge) = tx/tvol
                    eudg(4, nz, edge) = ty/tvol
                end do
                ! shared levels: take gradient straight from up/downwind triangle
                do nz = nzmin, nzmax-1
                    eudg(1:2, nz, edge) = tr_xy(1, nz, twork%edge_up_dn_tri(:, edge))
                    eudg(3:4, nz, edge) = tr_xy(2, nz, twork%edge_up_dn_tri(:, edge))
                end do
                ! not-shared lower levels of edge node 1
                do nz = nzmax, mesh%nlevels_nod2D(ednodes(1))-1
                    tvol = 0.0_WP; tx = 0.0_WP; ty = 0.0_WP
                    do k = 1, mesh%nod_in_elem2D_num(ednodes(1))
                        elem = mesh%nod_in_elem2D(k, ednodes(1))
                        if (mesh%nlevels(elem)-1 < nz .or. nz < mesh%ulevels(elem)) cycle
                        tvol = tvol + mesh%elem_area(elem)
                        tx = tx + tr_xy(1, nz, elem)*mesh%elem_area(elem)
                        ty = ty + tr_xy(2, nz, elem)*mesh%elem_area(elem)
                    end do
                    eudg(1, nz, edge) = tx/tvol
                    eudg(3, nz, edge) = ty/tvol
                end do
                ! not-shared lower levels of edge node 2
                do nz = nzmax, mesh%nlevels_nod2D(ednodes(2))-1
                    tvol = 0.0_WP; tx = 0.0_WP; ty = 0.0_WP
                    do k = 1, mesh%nod_in_elem2D_num(ednodes(2))
                        elem = mesh%nod_in_elem2D(k, ednodes(2))
                        if (mesh%nlevels(elem)-1 < nz .or. nz < mesh%ulevels(elem)) cycle
                        tvol = tvol + mesh%elem_area(elem)
                        tx = tx + tr_xy(1, nz, elem)*mesh%elem_area(elem)
                        ty = ty + tr_xy(2, nz, elem)*mesh%elem_area(elem)
                    end do
                    eudg(2, nz, edge) = tx/tvol
                    eudg(4, nz, edge) = ty/tvol
                end do
            !___ edge has only one triangle on the surface (boundary edge) __________
            else
                ! Only linear reconstruction part (standard Miura at both nodes)
                nzmin = mesh%ulevels_nod2D(ednodes(1))
                nzmax = mesh%nlevels_nod2D(ednodes(1))
                do nz = nzmin, nzmax-1
                    tvol = 0.0_WP; tx = 0.0_WP; ty = 0.0_WP
                    do k = 1, mesh%nod_in_elem2D_num(ednodes(1))
                        elem = mesh%nod_in_elem2D(k, ednodes(1))
                        if (mesh%nlevels(elem)-1 < nz .or. nz < mesh%ulevels(elem)) cycle
                        tvol = tvol + mesh%elem_area(elem)
                        tx = tx + tr_xy(1, nz, elem)*mesh%elem_area(elem)
                        ty = ty + tr_xy(2, nz, elem)*mesh%elem_area(elem)
                    end do
                    eudg(1, nz, edge) = tx/tvol
                    eudg(3, nz, edge) = ty/tvol
                end do
                nzmin = mesh%ulevels_nod2D(ednodes(2))
                nzmax = mesh%nlevels_nod2D(ednodes(2))
                do nz = nzmin, nzmax-1
                    tvol = 0.0_WP; tx = 0.0_WP; ty = 0.0_WP
                    do k = 1, mesh%nod_in_elem2D_num(ednodes(2))
                        elem = mesh%nod_in_elem2D(k, ednodes(2))
                        if (mesh%nlevels(elem)-1 < nz .or. nz < mesh%ulevels(elem)) cycle
                        tvol = tvol + mesh%elem_area(elem)
                        tx = tx + tr_xy(1, nz, elem)*mesh%elem_area(elem)
                        ty = ty + tr_xy(2, nz, elem)*mesh%elem_area(elem)
                    end do
                    eudg(2, nz, edge) = tx/tvol
                    eudg(4, nz, edge) = ty/tvol
                end do
            end if
        end do
    end subroutine fill_up_dn_grad

    !---------------------------------------------------------------------------
    subroutine muscl_node_grad(gnod, tr_xy, gnod_lo, gnod_hi, mesh, partit)
        ! Miura node-averaged tracer gradient, ONCE PER NODE:
        !   gnod(:,nz,n) = sum_{elem around n, wet at nz} tr_xy(:,nz,elem)*elem_area(elem)
        !                / sum_{same elems} elem_area(elem)
        ! for nz in [ulevels_nod2D(n), gnod_lo(n)) and [gnod_hi(n), nlevels_nod2D(n)),
        ! n = 1..nNodL (owned + halo): the levels on which an owned edge reads it
        ! (muscl_node_ranges). The other entries are NOT written: the caller zeroes gnod
        ! once at allocation, so the kernels read 0 outside a node's wet range (FESOM2's
        ! never-written entries), and the shared levels are never read.
        !
        ! WHY: fill_up_dn_grad uses exactly this average on the levels of an edge that are
        ! not shared by both up/downwind triangles and on edges without them, but evaluates
        ! it per EDGE, i.e. ~6 times per node, and stores the result in the 4-component
        ! edge array edge_up_dn_grad. The average depends on (node, level) only, so it is
        ! computed here once; the MUSCL kernels look it up on the fly together with tr_xy
        ! of the up/downwind triangle (docs/plans/completed/2026-10-06-muscl-onthefly.md).
        ! The loop order over nod_in_elem2D, the wet test and the tx/tvol division are
        ! fill_up_dn_grad's, so every value is bit-identical to the one the fill stores.
        ! A level with no wet element around the node (a node deeper than all its
        ! elements under bottom-at-vertices) gets 0 here; the fill stores 0/0 there.
        ! Neither is ever read by the kernels: they only visit levels wet in one of the
        ! edge's elements, which are elements around both edge nodes.
        ! Halo nodes: nod_in_elem2D is complete for halo nodes and tr_xy/elem_area are
        ! full-halo valid (exchange_elem_full), so no exchange of gnod is needed.
        real(kind=WP),  intent(inout) :: gnod(:,:,:)        ! (2, nl-1, nNodL)
        real(kind=WP),  intent(in)    :: tr_xy(:,:,:)       ! (2, nl-1, nElemF), full halo
        integer,        intent(in)    :: gnod_lo(:), gnod_hi(:)   ! (nNodL), muscl_node_ranges
        type(t_mesh),   intent(in)    :: mesh
        type(t_partit), intent(in), optional :: partit
        integer       :: n, nz, k, elem
        real(kind=WP) :: tvol, tx, ty
        integer       :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        do n = 1, nNodL
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n)-1
                if (nz >= gnod_lo(n) .and. nz < gnod_hi(n)) cycle      ! not read
                tvol = 0.0_WP; tx = 0.0_WP; ty = 0.0_WP
                do k = 1, mesh%nod_in_elem2D_num(n)
                    elem = mesh%nod_in_elem2D(k, n)
                    if (mesh%nlevels(elem)-1 < nz .or. nz < mesh%ulevels(elem)) cycle
                    tvol = tvol + mesh%elem_area(elem)
                    tx = tx + tr_xy(1, nz, elem)*mesh%elem_area(elem)
                    ty = ty + tr_xy(2, nz, elem)*mesh%elem_area(elem)
                end do
                if (tvol > 0.0_WP) then
                    gnod(1, nz, n) = tx/tvol
                    gnod(2, nz, n) = ty/tvol
                end if
            end do
        end do
    end subroutine muscl_node_grad

end module oce_muscl_adv
