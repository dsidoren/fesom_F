module mod_mesh_areas
    ! Core 2D mesh geometry, transcribed from FESOM2 v2.7.3 oce_mesh.F90 (cites
    ! inline). Arity-general: every FESOM2 `/3.0` vertex average becomes `/nv`
    ! (nv = elem2D_nnodes(elem)); for triangles nv==3 so the anchor is unchanged.
    !
    ! Provides what M1 (tracer advection) needs: elem_cos, metric_factor, elem_area,
    ! gradient_sca, edge_dxdy, edge_cross_dxdy, area/areasvol(+inv). gradient_vec
    ! (M2 momentum), coriolis (M2) and mesh_resolution smoothing (M4 GM) are deferred.
    !
    ! BYTE-FAITHFULNESS (M1 geometry byte-gate). FESOM2 splits this work across two
    ! routines whose ORDER of operations the bits depend on:
    !   mesh_areas (oce_mesh.F90:2157)         -> elem_area (UNSCALED), area accumulate
    !                                             (still UNSCALED), THEN one *r_earth^2
    !                                             on elem_area/area/areasvol together,
    !                                             then area_inv/areasvol_inv.
    !   mesh_auxiliary_arrays (oce_mesh.F90:2453) -> elem_cos/metric_factor (via
    !                                             elem_center), edge_dxdy,
    !                                             edge_cross_dxdy, gradient_sca
    !                                             (the latter uses the *scaled*
    !                                             elem_area). We reproduce that exact
    !                                             ordering and the elem_center /
    !                                             edge_center wrap arithmetic verbatim.
    use mod_precision,   only: WP, MP
    use mod_constants,   only: r_earth
    use mod_mesh,        only: t_mesh
    use mod_partit,      only: t_partit
    use mod_mesh_rotate, only: trim_cyclic, get_cyclic_length
    implicit none
    private
    public :: compute_geometry

contains

    subroutine compute_geometry(mesh, partit, cartesian)
        ! Same data-dependency order as FESOM2 mesh_areas + mesh_auxiliary_arrays:
        ! elem_cos before edge/gradient; elem_area UNSCALED before area accumulation;
        ! elem_area scaled (in compute_node_areas) before gradient_sca.
        type(t_mesh),   intent(inout) :: mesh
        type(t_partit), intent(in)    :: partit
        logical,        intent(in)    :: cartesian
        call compute_elem_metric(mesh, cartesian)   ! elem_cos, metric_factor
        call compute_elem_area(mesh, cartesian)     ! elem_area (UNSCALED radians^2)
        call compute_node_areas(mesh)               ! accumulate area, then scale *r_earth^2
        call compute_edge_geometry(mesh)            ! edge_dxdy, edge_cross_dxdy (uses elem_cos)
        call compute_gradient_sca(mesh)             ! uses SCALED elem_area + elem_cos
    end subroutine compute_geometry

    !--------------------------------------------------------------------------
    pure subroutine elem_center(mesh, n, cx, cy)
        ! FESOM2 oce_mesh.F90:2135-2155: wrap each vertex longitude relative to the
        ! element's MINIMUM longitude (amin), then cx = sum(lon)/nv ; cy = sum(lat)/nv.
        ! The wrap uses >= / < against cyclic_length/2 (NOT the symmetric trim_cyclic).
        type(t_mesh), intent(in)  :: mesh
        integer,      intent(in)  :: n
        real(kind=WP), intent(out) :: cx, cy
        integer :: j, nv
        real(kind=WP) :: ax(size(mesh%elem2D_nodes,1)), amin, cl
        nv = mesh%elem2D_nnodes(n)
        do j = 1, nv
            ax(j) = mesh%coord_nod2D(1, mesh%elem2D_nodes(j, n))
        end do
        amin = minval(ax(1:nv))
        cl = get_cyclic_length()
        do j = 1, nv
            if (ax(j) - amin >=  cl/2.0_WP) ax(j) = ax(j) - cl
            if (ax(j) - amin <  -cl/2.0_WP) ax(j) = ax(j) + cl
        end do
        ! Divide by the LITERAL 3.0_WP, NOT real(nv,WP): under -no-prec-div a runtime
        ! divisor uses a reciprocal approximation that differs from the compile-time
        ! 1/3 by 1 ULP -> breaks the geometry byte-gate. (Arity caveat, per the plan:
        ! the /3 -> /n_vert generalization is NOT bit-safe here. Quads must branch on
        ! nv with the matching literal; v1 is triangles only, nv==3 everywhere.)
        cx = sum(ax(1:nv)) / 3.0_WP
        cy = sum(mesh%coord_nod2D(2, mesh%elem2D_nodes(1:nv, n))) / 3.0_WP
    end subroutine elem_center

    pure subroutine edge_center(mesh, e, cx, cy)
        ! FESOM2 oce_mesh.F90:2115-2133: a = coord(n1), b = coord(n2); shift a(1) DOWN
        ! if a-b too positive, shift b(1) DOWN if a-b too negative; cx = 0.5*(a1+b1).
        type(t_mesh), intent(in)  :: mesh
        integer,      intent(in)  :: e
        real(kind=WP), intent(out) :: cx, cy
        real(kind=WP) :: a1, a2, b1, b2, cl
        a1 = mesh%coord_nod2D(1, mesh%edges(1, e)); a2 = mesh%coord_nod2D(2, mesh%edges(1, e))
        b1 = mesh%coord_nod2D(1, mesh%edges(2, e)); b2 = mesh%coord_nod2D(2, mesh%edges(2, e))
        cl = get_cyclic_length()
        if (a1 - b1 >  cl/2.0_WP) a1 = a1 - cl
        if (a1 - b1 < -cl/2.0_WP) b1 = b1 - cl
        cx = 0.5_WP * (a1 + b1)
        cy = 0.5_WP * (a2 + b2)
    end subroutine edge_center

    !--------------------------------------------------------------------------
    subroutine compute_elem_metric(mesh, cartesian)
        ! elem_cos = cos(lat_center); metric_factor = tan(lat_center)/r_earth.
        ! (oce_mesh.F90:2508-2528)
        type(t_mesh), intent(inout) :: mesh
        logical,      intent(in)    :: cartesian
        integer :: n
        real(kind=WP) :: cx, cy
        allocate(mesh%elem_cos(mesh%elem2D), mesh%metric_factor(mesh%elem2D))
        do n = 1, mesh%elem2D
            call elem_center(mesh, n, cx, cy)
            mesh%elem_cos(n)      = cos(cy)
            mesh%metric_factor(n) = tan(cy) / r_earth
        end do
        if (cartesian) then
            mesh%elem_cos = 1.0_MP; mesh%metric_factor = 0.0_MP
        end if
    end subroutine compute_elem_metric

    subroutine compute_elem_area(mesh, cartesian)
        ! oce_mesh.F90:2202-2214. ay = cos(sum(lat)/nv). elem_area is left UNSCALED
        ! here (radians^2); the *r_earth^2 happens in compute_node_areas, matching
        ! FESOM2's mesh_areas where the scaling is deferred to after area accumulation.
        type(t_mesh), intent(inout) :: mesh
        logical,      intent(in)    :: cartesian
        integer :: n, nv, n1
        real(kind=WP) :: ay, a1, a2, b1, b2
        allocate(mesh%elem_area(mesh%elem2D))
        do n = 1, mesh%elem2D
            nv = mesh%elem2D_nnodes(n)
            n1 = mesh%elem2D_nodes(1, n)
            ! literal 3.0_WP divisor (see elem_center: -no-prec-div arity caveat)
            ay = sum(mesh%coord_nod2D(2, mesh%elem2D_nodes(1:nv, n))) / 3.0_WP
            ay = cos(ay)
            if (cartesian) ay = 1.0_WP
            a1 = mesh%coord_nod2D(1, mesh%elem2D_nodes(2, n)) - mesh%coord_nod2D(1, n1)
            a2 = mesh%coord_nod2D(2, mesh%elem2D_nodes(2, n)) - mesh%coord_nod2D(2, n1)
            b1 = mesh%coord_nod2D(1, mesh%elem2D_nodes(3, n)) - mesh%coord_nod2D(1, n1)
            b2 = mesh%coord_nod2D(2, mesh%elem2D_nodes(3, n)) - mesh%coord_nod2D(2, n1)
            call trim_cyclic(a1); call trim_cyclic(b1)
            a1 = a1 * ay; b1 = b1 * ay
            mesh%elem_area(n) = 0.5_WP * abs(a1*b2 - b1*a2)
        end do
    end subroutine compute_elem_area

    subroutine compute_gradient_sca(mesh)
        ! Linear shape-function gradient coefficients (oce_mesh.F90:2619-2641).
        ! Uses the SCALED elem_area (dfactor = -0.5*r_earth/elem_area).
        type(t_mesh), intent(inout) :: mesh
        integer :: e, n1, n2, n3
        real(kind=WP) :: dX31, dX21, dY31, dY21, dfac
        allocate(mesh%gradient_sca(2*size(mesh%elem2D_nodes,1), mesh%elem2D))  ! (2*MAX_NV, elem2D)
        mesh%gradient_sca = 0.0_MP
        do e = 1, mesh%elem2D
            n1 = mesh%elem2D_nodes(1, e); n2 = mesh%elem2D_nodes(2, e); n3 = mesh%elem2D_nodes(3, e)
            dX31 = mesh%coord_nod2D(1, n3) - mesh%coord_nod2D(1, n1); call trim_cyclic(dX31)
            dX31 = mesh%elem_cos(e) * dX31
            dX21 = mesh%coord_nod2D(1, n2) - mesh%coord_nod2D(1, n1); call trim_cyclic(dX21)
            dX21 = mesh%elem_cos(e) * dX21
            dY31 = mesh%coord_nod2D(2, n3) - mesh%coord_nod2D(2, n1)
            dY21 = mesh%coord_nod2D(2, n2) - mesh%coord_nod2D(2, n1)
            dfac = -0.5_WP * r_earth / mesh%elem_area(e)
            mesh%gradient_sca(1, e) = (-dY31 + dY21) * dfac
            mesh%gradient_sca(2, e) = dY31 * dfac
            mesh%gradient_sca(3, e) = -dY21 * dfac
            mesh%gradient_sca(4, e) = (dX31 - dX21) * dfac
            mesh%gradient_sca(5, e) = -dX31 * dfac
            mesh%gradient_sca(6, e) = dX21 * dfac
        end do
    end subroutine compute_gradient_sca

    subroutine compute_edge_geometry(mesh)
        ! edge_dxdy (along-edge, radians) + edge_cross_dxdy (edge-center to elem
        ! centers, metres). oce_mesh.F90:2534-2573.
        type(t_mesh), intent(inout) :: mesh
        integer :: n, el1, el2
        real(kind=WP) :: a1, a2, ecx, ecy, cx, cy, b1, b2
        allocate(mesh%edge_dxdy(2, mesh%edge2D), mesh%edge_cross_dxdy(4, mesh%edge2D))
        do n = 1, mesh%edge2D
            a1 = mesh%coord_nod2D(1, mesh%edges(2, n)) - mesh%coord_nod2D(1, mesh%edges(1, n))
            a2 = mesh%coord_nod2D(2, mesh%edges(2, n)) - mesh%coord_nod2D(2, mesh%edges(1, n))
            call trim_cyclic(a1)
            mesh%edge_dxdy(1, n) = a1; mesh%edge_dxdy(2, n) = a2
        end do
        do n = 1, mesh%edge2D
            call edge_center(mesh, n, ecx, ecy)
            el1 = mesh%edge_tri(1, n); el2 = mesh%edge_tri(2, n)
            call elem_center(mesh, el1, cx, cy)
            b1 = cx - ecx; b2 = cy - ecy; call trim_cyclic(b1)
            b1 = b1 * mesh%elem_cos(el1)
            mesh%edge_cross_dxdy(1, n) = b1 * r_earth
            mesh%edge_cross_dxdy(2, n) = b2 * r_earth
            if (el2 > 0) then
                call elem_center(mesh, el2, cx, cy)
                b1 = cx - ecx; b2 = cy - ecy; call trim_cyclic(b1)
                b1 = b1 * mesh%elem_cos(el2)
                mesh%edge_cross_dxdy(3, n) = b1 * r_earth
                mesh%edge_cross_dxdy(4, n) = b2 * r_earth
            else
                mesh%edge_cross_dxdy(3, n) = 0.0_WP; mesh%edge_cross_dxdy(4, n) = 0.0_WP
            end if
        end do
    end subroutine compute_edge_geometry

    subroutine compute_node_areas(mesh)
        ! Control-volume area per level (oce_mesh.F90:2252-2351). area(nz,n) gathers
        ! elem_area/nv from adjacent elements deep enough to reach level nz. The
        ! accumulation runs on UNSCALED elem_area; then elem_area, area and areasvol
        ! are all multiplied by r_earth^2 together (a single deferred scaling, as in
        ! FESOM2 mesh_areas:2313-2315) so the per-node sums round identically.
        type(t_mesh), intent(inout) :: mesh
        integer :: n, j, elem, nz, nzmin, nzmax
        allocate(mesh%area(mesh%nl, mesh%nod2D), mesh%area_inv(mesh%nl, mesh%nod2D))
        allocate(mesh%areasvol(mesh%nl, mesh%nod2D), mesh%areasvol_inv(mesh%nl, mesh%nod2D))
        mesh%area = 0.0_MP
        do n = 1, mesh%nod2D
            do j = 1, mesh%nod_in_elem2D_num(n)
                elem = mesh%nod_in_elem2D(j, n)
                nzmin = mesh%ulevels(elem)
                nzmax = mesh%nlevels(elem) - 1
                do nz = nzmin, nzmax
                    ! literal 3.0_MP divisor (FESOM2 mesh_areas:2266; -no-prec-div
                    ! arity caveat — see elem_center). Triangles only at the anchor.
                    mesh%area(nz, n) = mesh%area(nz, n) &
                        + mesh%elem_area(elem) / 3.0_MP
                end do
            end do
        end do
        ! non-cavity: "mid" cell area == upper-edge area
        mesh%areasvol = 0.0_MP
        do n = 1, mesh%nod2D
            nzmin = mesh%ulevels_nod2D(n)
            nzmax = mesh%nlevels_nod2D(n) - 1
            do nz = nzmin, nzmax
                mesh%areasvol(nz, n) = mesh%area(nz, n)
            end do
        end do
        ! deferred single scaling to physical m^2 (mesh_areas:2313-2315)
        mesh%elem_area = mesh%elem_area * r_earth * r_earth
        mesh%area      = mesh%area      * r_earth * r_earth
        mesh%areasvol  = mesh%areasvol  * r_earth * r_earth
        ! inverse areas (mesh_areas:2321-2351); non-cavity areasvol_inv == area_inv
        mesh%area_inv = 0.0_MP
        do n = 1, mesh%nod2D
            nzmin = mesh%ulevels_nod2D(n)
            nzmax = mesh%nlevels_nod2D(n)
            do nz = nzmin, nzmax
                if (mesh%area(nz, n) > 0.0_MP) then
                    mesh%area_inv(nz, n) = 1.0_MP / mesh%area(nz, n)
                else
                    mesh%area_inv(nz, n) = 0.0_MP
                end if
            end do
        end do
        mesh%areasvol_inv = mesh%area_inv
        mesh%ocean_area = sum(mesh%area(1, 1:mesh%nod2D))
    end subroutine compute_node_areas

end module mod_mesh_areas
