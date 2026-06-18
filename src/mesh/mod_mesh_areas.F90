module mod_mesh_areas
    ! Core 2D mesh geometry, transcribed from FESOM2 v2.7.3 oce_mesh.F90 (cites
    ! inline). Arity-general: every FESOM2 `/3.0` vertex average becomes `/nv`
    ! (nv = elem2D_nnodes(elem)); for triangles nv==3 so the anchor is unchanged.
    !
    ! Provides what M1 (tracer advection) needs: elem_cos, metric_factor, elem_area,
    ! gradient_sca, edge_dxdy, edge_cross_dxdy, area/areasvol(+inv). gradient_vec
    ! (M2 momentum) and mesh_resolution smoothing (M4 GM) are deferred.
    use mod_precision,   only: WP, MP
    use mod_constants,   only: r_earth, pi
    use mod_mesh,        only: t_mesh
    use mod_partit,      only: t_partit
    use mod_mesh_rotate, only: trim_cyclic
    implicit none
    private
    public :: compute_geometry

contains

    subroutine compute_geometry(mesh, partit, cartesian)
        type(t_mesh),   intent(inout) :: mesh
        type(t_partit), intent(in)    :: partit
        logical,        intent(in)    :: cartesian
        call compute_elem_metric(mesh, cartesian)
        call compute_elem_area(mesh, cartesian)
        call compute_gradient_sca(mesh)
        call compute_edge_geometry(mesh)
        call compute_node_areas(mesh)
    end subroutine compute_geometry

    !--------------------------------------------------------------------------
    pure subroutine elem_center(mesh, n, cx, cy)
        type(t_mesh), intent(in)  :: mesh
        integer,      intent(in)  :: n
        real(kind=WP), intent(out) :: cx, cy
        integer :: j, nv, n1
        real(kind=WP) :: dx, sx, sy
        nv = mesh%elem2D_nnodes(n)
        n1 = mesh%elem2D_nodes(1, n)
        sx = 0.0_WP; sy = 0.0_WP
        do j = 1, nv
            dx = mesh%coord_nod2D(1, mesh%elem2D_nodes(j, n)) - mesh%coord_nod2D(1, n1)
            call trim_cyclic(dx)
            sx = sx + dx
            sy = sy + mesh%coord_nod2D(2, mesh%elem2D_nodes(j, n))
        end do
        cx = mesh%coord_nod2D(1, n1) + sx / real(nv, WP)
        cy = sy / real(nv, WP)
    end subroutine elem_center

    pure subroutine edge_center(mesh, e, cx, cy)
        type(t_mesh), intent(in)  :: mesh
        integer,      intent(in)  :: e
        real(kind=WP), intent(out) :: cx, cy
        real(kind=WP) :: dx
        dx = mesh%coord_nod2D(1, mesh%edges(2, e)) - mesh%coord_nod2D(1, mesh%edges(1, e))
        call trim_cyclic(dx)
        cx = mesh%coord_nod2D(1, mesh%edges(1, e)) + 0.5_WP * dx
        cy = 0.5_WP * (mesh%coord_nod2D(2, mesh%edges(1, e)) + mesh%coord_nod2D(2, mesh%edges(2, e)))
    end subroutine edge_center

    !--------------------------------------------------------------------------
    subroutine compute_elem_metric(mesh, cartesian)
        ! elem_cos = cos(lat_center); metric_factor = tan(lat_center)/r_earth.
        ! (oce_mesh.F90:2505-2528)
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
        ! oce_mesh.F90:2201-2213, 2313. ay = sum(lat)/nv (the arity site).
        type(t_mesh), intent(inout) :: mesh
        logical,      intent(in)    :: cartesian
        integer :: n, j, nv, n1
        real(kind=WP) :: ay, a1, a2, b1, b2
        allocate(mesh%elem_area(mesh%elem2D))
        do n = 1, mesh%elem2D
            nv = mesh%elem2D_nnodes(n)
            n1 = mesh%elem2D_nodes(1, n)
            ay = 0.0_WP
            do j = 1, nv
                ay = ay + mesh%coord_nod2D(2, mesh%elem2D_nodes(j, n))
            end do
            ay = ay / real(nv, WP)
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
        mesh%elem_area = mesh%elem_area * r_earth * r_earth
    end subroutine compute_elem_area

    subroutine compute_gradient_sca(mesh)
        ! Linear shape-function gradient coefficients (oce_mesh.F90:2619-2640).
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
        ! centers, metres). oce_mesh.F90:2540-2570.
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
        ! Control-volume area per level (oce_mesh.F90:2252-2328). area(nz,n) gathers
        ! elem_area/nv from adjacent elements deep enough to reach level nz.
        type(t_mesh), intent(inout) :: mesh
        integer :: n, j, elem, nz, nzmax
        allocate(mesh%area(mesh%nl, mesh%nod2D), mesh%area_inv(mesh%nl, mesh%nod2D))
        allocate(mesh%areasvol(mesh%nl, mesh%nod2D), mesh%areasvol_inv(mesh%nl, mesh%nod2D))
        mesh%area = 0.0_MP
        do n = 1, mesh%nod2D
            do j = 1, mesh%nod_in_elem2D_num(n)
                elem = mesh%nod_in_elem2D(j, n)
                do nz = mesh%ulevels(elem), mesh%nlevels(elem) - 1
                    mesh%area(nz, n) = mesh%area(nz, n) &
                        + mesh%elem_area(elem) / real(mesh%elem2D_nnodes(elem), MP)
                end do
            end do
        end do
        mesh%areasvol = mesh%area                       ! linfs / no cavity
        mesh%area_inv = 0.0_MP; mesh%areasvol_inv = 0.0_MP
        do n = 1, mesh%nod2D
            nzmax = mesh%nlevels_nod2D(n) - 1
            do nz = mesh%ulevels_nod2D(n), nzmax
                if (mesh%area(nz, n) > 0.0_MP)     mesh%area_inv(nz, n)     = 1.0_MP / mesh%area(nz, n)
                if (mesh%areasvol(nz, n) > 0.0_MP) mesh%areasvol_inv(nz, n) = 1.0_MP / mesh%areasvol(nz, n)
            end do
        end do
        mesh%ocean_area = sum(mesh%area(1, 1:mesh%nod2D))
    end subroutine compute_node_areas

end module mod_mesh_areas
