module mod_mesh_analytic
    ! In-memory analytic mesh generator (no file I/O), for fast tests + the
    ! fesom_analytic driver. A regular nx*ny Cartesian grid, 2 triangles per cell,
    ! nl uniform levels. cartesian=.true. so the metric is off. Synthesizes a
    ! 1-rank identity partition too. (Doubly-periodic variant deferred; closed only.)
    use mod_precision, only: WP, MP
    use mod_constants, only: r_earth
    use mod_mesh,      only: t_mesh, MAX_NV, MAX_ADJACENT
    use mod_partit,    only: t_partit
    use mod_mesh_areas,  only: compute_geometry
    use mod_mesh_rotate, only: init_mesh_rotation
    use mod_mesh_read,   only: enforce_cw_orientation
    implicit none
    private
    public :: generate_analytic_mesh, build_edges

contains

    subroutine generate_analytic_mesh(mesh, partit, nx, ny, nl, Lx, Ly, max_depth)
        type(t_mesh),   intent(inout) :: mesh
        type(t_partit), intent(inout) :: partit
        integer,        intent(in)    :: nx, ny, nl
        real(kind=WP),  intent(in)    :: Lx, Ly, max_depth
        integer :: i, j, nod, el, k
        real(kind=WP) :: dx, dy

        mesh%nod2D  = nx * ny
        mesh%elem2D = 2 * (nx-1) * (ny-1)
        mesh%nl     = nl
        ! Coords stored radians-like (physical metres / r_earth) so the shared
        ! geometry pipeline's `*r_earth^2` (sphere convention) yields physical m^2.
        dx = Lx / real(nx-1, WP) / r_earth; dy = Ly / real(ny-1, WP) / r_earth

        ! ---- nodes (Cartesian) ----
        allocate(mesh%coord_nod2D(2, mesh%nod2D), mesh%geo_coord_nod2D(2, mesh%nod2D))
        allocate(mesh%bc_index_nod2D(mesh%nod2D))
        mesh%bc_index_nod2D = 0
        do j = 1, ny
            do i = 1, nx
                nod = (j-1)*nx + i
                mesh%coord_nod2D(1, nod) = real(i-1, WP) * dx
                mesh%coord_nod2D(2, nod) = real(j-1, WP) * dy
            end do
        end do
        mesh%geo_coord_nod2D = mesh%coord_nod2D

        ! ---- elements (2 triangles per cell) ----
        allocate(mesh%elem2D_nodes(MAX_NV, mesh%elem2D), mesh%elem2D_nnodes(mesh%elem2D))
        mesh%elem2D_nodes = 0; mesh%elem2D_nnodes = 3
        el = 0
        do j = 1, ny-1
            do i = 1, nx-1
                el = el + 1
                mesh%elem2D_nodes(1, el) = (j-1)*nx + i
                mesh%elem2D_nodes(2, el) = (j-1)*nx + i + 1
                mesh%elem2D_nodes(3, el) = j*nx + i
                el = el + 1
                mesh%elem2D_nodes(1, el) = (j-1)*nx + i + 1
                mesh%elem2D_nodes(2, el) = j*nx + i + 1
                mesh%elem2D_nodes(3, el) = j*nx + i
            end do
        end do

        ! ---- vertical: uniform levels 0 .. -max_depth ----
        allocate(mesh%zbar(nl), mesh%Z(nl-1))
        do k = 1, nl
            mesh%zbar(k) = -max_depth * real(k-1, MP) / real(nl-1, MP)
        end do
        mesh%Z = 0.5_MP * (mesh%zbar(1:nl-1) + mesh%zbar(2:nl))
        allocate(mesh%nlevels(mesh%elem2D), mesh%nlevels_nod2D(mesh%nod2D))
        allocate(mesh%ulevels(mesh%elem2D), mesh%ulevels_nod2D(mesh%nod2D))
        allocate(mesh%ulevels_nod2D_max(mesh%nod2D), mesh%nlevels_nod2D_min(mesh%nod2D))
        allocate(mesh%elem_depth(mesh%elem2D), mesh%depth(mesh%nod2D))
        mesh%nlevels = nl; mesh%nlevels_nod2D = nl
        mesh%ulevels = 1; mesh%ulevels_nod2D = 1; mesh%ulevels_nod2D_max = 1
        mesh%nlevels_nod2D_min = nl; mesh%elem_depth = -max_depth; mesh%depth = -max_depth

        ! ---- orientation, topology, 1-rank partition, geometry ----
        ! Cartesian coords are in metres; disable cyclic wrapping (huge cyclic
        ! length) so trim_cyclic inside enforce/compute is a no-op here. CW
        ! enforcement must precede geometry (gradient_sca assumes CW node order).
        call init_mesh_rotation(0.0_WP, 0.0_WP, 0.0_WP, 1.0e12_WP)
        call enforce_cw_orientation(mesh, el)
        call build_nod_in_elem_local(mesh)
        call build_edges(mesh)
        call synth_partit(mesh, partit)
        call compute_geometry(mesh, partit, cartesian=.true.)
    end subroutine generate_analytic_mesh

    subroutine build_nod_in_elem_local(mesh)
        type(t_mesh), intent(inout) :: mesh
        integer :: n, j, node
        allocate(mesh%nod_in_elem2D_num(mesh%nod2D))
        allocate(mesh%nod_in_elem2D(MAX_ADJACENT, mesh%nod2D))
        mesh%nod_in_elem2D_num = 0; mesh%nod_in_elem2D = 0
        do n = 1, mesh%elem2D
            do j = 1, mesh%elem2D_nnodes(n)
                node = mesh%elem2D_nodes(j, n)
                mesh%nod_in_elem2D_num(node) = mesh%nod_in_elem2D_num(node) + 1
                mesh%nod_in_elem2D(mesh%nod_in_elem2D_num(node), node) = n
            end do
        end do
    end subroutine build_nod_in_elem_local

    subroutine build_edges(mesh)
        ! Generic edge + edge_tri builder from elem2D_nodes (triangles). edges(:,e)
        ! directed lo->hi; edge_tri(1)=first elem seen, edge_tri(2)=second (or 0 at
        ! boundary). Exact left/right orientation is not FESOM2's (analytic only).
        type(t_mesh), intent(inout) :: mesh
        integer, allocatable :: head_hi(:,:), head_eid(:,:), head_n(:)
        integer :: n, j, a, b, lo, hi, k, eid, nedge
        integer, parameter :: PAIR(2,3) = reshape([1,2, 2,3, 3,1], [2,3])
        allocate(head_hi(MAX_ADJACENT, mesh%nod2D), head_eid(MAX_ADJACENT, mesh%nod2D))
        allocate(head_n(mesh%nod2D))
        head_n = 0
        allocate(mesh%edges(2, 3*mesh%elem2D), mesh%edge_tri(2, 3*mesh%elem2D))
        mesh%edge_tri = 0
        nedge = 0
        do n = 1, mesh%elem2D
            do j = 1, 3
                a = mesh%elem2D_nodes(PAIR(1,j), n); b = mesh%elem2D_nodes(PAIR(2,j), n)
                lo = min(a,b); hi = max(a,b)
                eid = 0
                do k = 1, head_n(lo)
                    if (head_hi(k, lo) == hi) then; eid = head_eid(k, lo); exit; end if
                end do
                if (eid == 0) then
                    nedge = nedge + 1
                    mesh%edges(1, nedge) = lo; mesh%edges(2, nedge) = hi
                    mesh%edge_tri(1, nedge) = n
                    head_n(lo) = head_n(lo) + 1
                    head_hi(head_n(lo), lo) = hi; head_eid(head_n(lo), lo) = nedge
                else
                    mesh%edge_tri(2, eid) = n
                end if
            end do
        end do
        mesh%edge2D = nedge; mesh%edge2D_in = count(mesh%edge_tri(2, 1:nedge) > 0)
        ! shrink to actual edge count
        mesh%edges     = mesh%edges(:, 1:nedge)
        mesh%edge_tri  = mesh%edge_tri(:, 1:nedge)
        deallocate(head_hi, head_eid, head_n)
    end subroutine build_edges

    subroutine synth_partit(mesh, partit)
        type(t_mesh),   intent(in)    :: mesh
        type(t_partit), intent(inout) :: partit
        integer :: i
        partit%myDim_nod2D = mesh%nod2D; partit%eDim_nod2D = 0
        partit%myDim_elem2D = mesh%elem2D; partit%eDim_elem2D = 0; partit%eXDim_elem2D = 0
        partit%myDim_edge2D = mesh%edge2D; partit%eDim_edge2D = 0
        if (allocated(partit%myList_nod2D)) deallocate(partit%myList_nod2D)
        allocate(partit%myList_nod2D(mesh%nod2D)); partit%myList_nod2D = [(i, i=1, mesh%nod2D)]
        partit%com_nod2D%rPEnum = 0; partit%com_nod2D%sPEnum = 0
        partit%com_elem2D%rPEnum = 0; partit%com_elem2D%sPEnum = 0
        partit%com_elem2D_full%rPEnum = 0; partit%com_elem2D_full%sPEnum = 0
    end subroutine synth_partit

end module mod_mesh_analytic
