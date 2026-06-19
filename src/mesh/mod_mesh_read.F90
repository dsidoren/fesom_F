module mod_mesh_read
    ! Read a partitioned FESOM mesh into t_mesh (1-rank: local == global; the
    ! multi-rank remap is deferred to M2.12). File formats + topology transcribed
    ! from FESOM2 v2.7.3 oce_mesh.F90 (file:line cites inline). Builds nod_in_elem2D,
    ! enforces clockwise element orientation (test_tri), sets up the vertical levels.
    use mod_precision,   only: WP, MP
    use mod_constants,   only: rad
    use mod_mesh,        only: t_mesh, MAX_NV, MAX_ADJACENT
    use mod_partit,      only: t_partit
    use mod_mesh_rotate, only: init_mesh_rotation, g2r, r2g, trim_cyclic
    implicit none
    private
    public :: read_mesh, enforce_cw_orientation

contains

    subroutine read_mesh(mesh, partit, mesh_dir, alpha_deg, beta_deg, gamma_deg, &
                         cyclic_deg, force_rotation, n_cw_swaps)
        type(t_mesh),     intent(inout) :: mesh
        type(t_partit),   intent(in)    :: partit
        character(len=*), intent(in)    :: mesh_dir
        real(kind=WP),    intent(in)    :: alpha_deg, beta_deg, gamma_deg, cyclic_deg
        logical,          intent(in)    :: force_rotation
        integer,          intent(out)   :: n_cw_swaps
        integer :: u, ios, i, n, id, n1, n2, n3, bnd
        real(kind=WP) :: lon, lat, gx, gy

        if (partit%npes /= 1) then
            write(*,'(a)') 'read_mesh: only 1-rank (global) reading implemented (M0.7)'
            error stop 1
        end if
        call init_mesh_rotation(alpha_deg, beta_deg, gamma_deg, cyclic_deg)

        ! ---- nodes: nod2d.out (id lon_deg lat_deg bnd) ----
        open(newunit=u, file=trim(mesh_dir)//'/nod2d.out', status='old', action='read', iostat=ios)
        read(u,*) mesh%nod2D
        allocate(mesh%coord_nod2D(2, mesh%nod2D), mesh%geo_coord_nod2D(2, mesh%nod2D))
        allocate(mesh%bc_index_nod2D(mesh%nod2D))
        do i = 1, mesh%nod2D
            read(u,*) id, lon, lat, bnd
            gx = lon * rad; gy = lat * rad           ! geographic, radians
            if (force_rotation) then
                call g2r(gx, gy, mesh%coord_nod2D(1, id), mesh%coord_nod2D(2, id))
            else
                mesh%coord_nod2D(1, id) = gx; mesh%coord_nod2D(2, id) = gy
            end if
            mesh%bc_index_nod2D(id) = bnd
        end do
        close(u)
        ! geo_coord_nod2D = r2g(coord_nod2D)  (oce_mesh.F90:2480-2492)
        do i = 1, mesh%nod2D
            call r2g(mesh%geo_coord_nod2D(1,i), mesh%geo_coord_nod2D(2,i), &
                     mesh%coord_nod2D(1,i), mesh%coord_nod2D(2,i))
        end do

        ! ---- elements: elem2d.out (n1 n2 n3) ----
        open(newunit=u, file=trim(mesh_dir)//'/elem2d.out', status='old', action='read', iostat=ios)
        read(u,*) mesh%elem2D
        allocate(mesh%elem2D_nodes(MAX_NV, mesh%elem2D), mesh%elem2D_nnodes(mesh%elem2D))
        mesh%elem2D_nodes = 0
        do i = 1, mesh%elem2D
            read(u,*) n1, n2, n3
            mesh%elem2D_nodes(1, i) = n1; mesh%elem2D_nodes(2, i) = n2; mesh%elem2D_nodes(3, i) = n3
            mesh%elem2D_nnodes(i) = 3
        end do
        close(u)

        ! ---- edges: edgenum.out, edges.out, edge_tri.out ----
        open(newunit=u, file=trim(mesh_dir)//'/edgenum.out', status='old', action='read', iostat=ios)
        read(u,*) mesh%edge2D
        read(u,*, iostat=ios) mesh%edge2D_in
        if (ios /= 0) mesh%edge2D_in = mesh%edge2D
        close(u)
        allocate(mesh%edges(2, mesh%edge2D), mesh%edge_tri(2, mesh%edge2D))
        open(newunit=u, file=trim(mesh_dir)//'/edges.out', status='old', action='read', iostat=ios)
        do i = 1, mesh%edge2D
            read(u,*) mesh%edges(1, i), mesh%edges(2, i)
        end do
        close(u)
        open(newunit=u, file=trim(mesh_dir)//'/edge_tri.out', status='old', action='read', iostat=ios)
        do i = 1, mesh%edge2D
            read(u,*) mesh%edge_tri(1, i), mesh%edge_tri(2, i)
        end do
        close(u)
        ! Boundary edges store -999 (no 2nd triangle) on disk; FESOM2 load_edges
        ! (oce_mesh.F90:1885-1887) zeroes the negatives. el(2)==0 marks a boundary edge.
        where (mesh%edge_tri < 0) mesh%edge_tri = 0

        ! ---- vertical level counts: elvls.out (elem), nlvls.out (node) ----
        allocate(mesh%nlevels(mesh%elem2D), mesh%nlevels_nod2D(mesh%nod2D))
        open(newunit=u, file=trim(mesh_dir)//'/elvls.out', status='old', action='read', iostat=ios)
        do i = 1, mesh%elem2D
            read(u,*) mesh%nlevels(i)
        end do
        close(u)
        open(newunit=u, file=trim(mesh_dir)//'/nlvls.out', status='old', action='read', iostat=ios)
        do i = 1, mesh%nod2D
            read(u,*) mesh%nlevels_nod2D(i)
        end do
        close(u)

        ! ---- aux3d.out: nl, zbar(nl), per-node depth ----
        open(newunit=u, file=trim(mesh_dir)//'/aux3d.out', status='old', action='read', iostat=ios)
        read(u,*) mesh%nl
        allocate(mesh%zbar(mesh%nl), mesh%Z(mesh%nl-1), mesh%depth(mesh%nod2D))
        read(u,*) mesh%zbar
        if (mesh%zbar(2) > 0.0_MP) mesh%zbar = -mesh%zbar            ! oce_mesh.F90:578
        mesh%Z = 0.5_MP * (mesh%zbar(1:mesh%nl-1) + mesh%zbar(2:mesh%nl))  ! :580-581
        do i = 1, mesh%nod2D
            read(u,*, iostat=ios) mesh%depth(i)
            if (ios /= 0) then; mesh%depth(i) = 0.0_MP; end if
        end do
        close(u)

        ! ---- topology + orientation + vertical structure ----
        call build_nod_in_elem(mesh)
        call enforce_cw_orientation(mesh, n_cw_swaps)
        call setup_vertical(mesh)
    end subroutine read_mesh

    subroutine build_nod_in_elem(mesh)
        ! node->element adjacency, dense (oce_mesh.F90:2010-2044), 2-pass.
        type(t_mesh), intent(inout) :: mesh
        integer :: n, j, node, kmax
        allocate(mesh%nod_in_elem2D_num(mesh%nod2D))
        mesh%nod_in_elem2D_num = 0
        do n = 1, mesh%elem2D
            do j = 1, mesh%elem2D_nnodes(n)
                node = mesh%elem2D_nodes(j, n)
                mesh%nod_in_elem2D_num(node) = mesh%nod_in_elem2D_num(node) + 1
            end do
        end do
        kmax = maxval(mesh%nod_in_elem2D_num)
        if (kmax > MAX_ADJACENT) then
            write(*,'(a,i0)') 'build_nod_in_elem: MAX_ADJACENT exceeded: ', kmax
            error stop 1
        end if
        allocate(mesh%nod_in_elem2D(MAX_ADJACENT, mesh%nod2D))
        mesh%nod_in_elem2D = 0
        mesh%nod_in_elem2D_num = 0
        do n = 1, mesh%elem2D
            do j = 1, mesh%elem2D_nnodes(n)
                node = mesh%elem2D_nodes(j, n)
                mesh%nod_in_elem2D_num(node) = mesh%nod_in_elem2D_num(node) + 1
                mesh%nod_in_elem2D(mesh%nod_in_elem2D_num(node), node) = n
            end do
        end do
    end subroutine build_nod_in_elem

    subroutine enforce_cw_orientation(mesh, n_swaps)
        ! test_tri (oce_mesh.F90:1679-1729): r = b1*c2 - b2*c1; if r>0 swap nodes 2,3.
        type(t_mesh), intent(inout) :: mesh
        integer,      intent(out)   :: n_swaps
        integer :: n, tmp
        real(kind=WP) :: a1, a2, b1, b2, c1, c2, r
        n_swaps = 0
        do n = 1, mesh%elem2D
            a1 = mesh%coord_nod2D(1, mesh%elem2D_nodes(1, n))
            a2 = mesh%coord_nod2D(2, mesh%elem2D_nodes(1, n))
            b1 = mesh%coord_nod2D(1, mesh%elem2D_nodes(2, n)) - a1
            b2 = mesh%coord_nod2D(2, mesh%elem2D_nodes(2, n)) - a2
            c1 = mesh%coord_nod2D(1, mesh%elem2D_nodes(3, n)) - a1
            c2 = mesh%coord_nod2D(2, mesh%elem2D_nodes(3, n)) - a2
            call trim_cyclic(b1)
            call trim_cyclic(c1)
            r = b1*c2 - b2*c1
            if (r > 0.0_WP) then
                tmp = mesh%elem2D_nodes(2, n)
                mesh%elem2D_nodes(2, n) = mesh%elem2D_nodes(3, n)
                mesh%elem2D_nodes(3, n) = tmp
                n_swaps = n_swaps + 1
            end if
        end do
    end subroutine enforce_cw_orientation

    subroutine setup_vertical(mesh)
        ! ulevels (no cavity => 1); nlevels_nod2D_min = min over adjacent elements;
        ! elem_depth = zbar(nlevels). (oce_mesh.F90:1103-1106, 1656-1665)
        type(t_mesh), intent(inout) :: mesh
        integer :: n, k
        allocate(mesh%ulevels(mesh%elem2D), mesh%ulevels_nod2D(mesh%nod2D))
        allocate(mesh%ulevels_nod2D_max(mesh%nod2D), mesh%nlevels_nod2D_min(mesh%nod2D))
        allocate(mesh%elem_depth(mesh%elem2D))
        mesh%ulevels = 1; mesh%ulevels_nod2D = 1; mesh%ulevels_nod2D_max = 1
        do n = 1, mesh%nod2D
            k = mesh%nod_in_elem2D_num(n)
            mesh%nlevels_nod2D_min(n) = minval(mesh%nlevels(mesh%nod_in_elem2D(1:k, n)))
        end do
        do n = 1, mesh%elem2D
            mesh%elem_depth(n) = mesh%zbar(mesh%nlevels(n))
        end do
    end subroutine setup_vertical

end module mod_mesh_read
