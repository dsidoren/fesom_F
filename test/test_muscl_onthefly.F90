module muscl_oracle
    ! THE ORACLE KERNELS: FESOM3's MUSCL / MFCT horizontal-flux kernels as they were before
    ! the on-the-fly refactor (verbatim, commit 48cc272: they read the stored per-edge
    ! up/downwind gradient edge_up_dn_grad built by fill_up_dn_grad). The model no longer
    ! has them; test_muscl_onthefly part F compares the production kernels with these.
    use mod_precision,   only: WP
    use mod_mesh,        only: t_mesh
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: owned_bounds
    implicit none
    private
    public :: oracle_muscl, oracle_mfct
contains

    subroutine oracle_muscl(vel, ttf, mesh, num_ord, flux, edge_up_dn_grad, nboundary_lay, o_init_zero, partit)
        ! MUSCL horizontal flux (oce_adv_tra_hor.F90:261-542). num_ord = fraction of
        ! 4th-order (centered) contribution; (1-num_ord) is 3rd-order upwind. The
        ! per-node clamp c_lo = max(sign(1,nboundary_lay-nz),0) switches off the
        ! linear-reconstruction increment below a node's boundary layer.
        ! M2.12b: optional partit -> loop over OWNED edges (myDim_edge2D).
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: num_ord
        real(kind=WP), intent(in)    :: ttf(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: vel(2, mesh%nl-1, mesh%elem2D)
        real(kind=WP), intent(inout) :: flux(mesh%nl-1, mesh%edge2D)
        integer,       intent(in)    :: nboundary_lay(mesh%nod2D)
        real(kind=WP), intent(in)    :: edge_up_dn_grad(4, mesh%nl-1, mesh%edge2D)
        logical, optional, intent(in) :: o_init_zero
        type(t_partit), intent(in), optional :: partit
        logical       :: l_init_zero
        real(kind=WP) :: deltaX1, deltaY1, deltaX2, deltaY2, vflux
        integer       :: el(2), enodes(2), nz, edge, nu12, nl12, nl1, nl2, nu1, nu2
        integer       :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        l_init_zero = .true.
        if (present(o_init_zero)) l_init_zero = o_init_zero
        if (l_init_zero) then
            do edge = 1, nEdgeO
                flux(:, edge) = 0.0_WP
            end do
        end if

        do edge = 1, nEdgeO
            enodes = mesh%edges(:, edge)
            el     = mesh%edge_tri(:, edge)
            nl1    = mesh%nlevels(el(1)) - 1
            nu1    = mesh%ulevels(el(1))
            deltaX1 = mesh%edge_cross_dxdy(1, edge)
            deltaY1 = mesh%edge_cross_dxdy(2, edge)
            nl2 = 0; nu2 = 0
            if (el(2) > 0) then
                deltaX2 = mesh%edge_cross_dxdy(3, edge)
                deltaY2 = mesh%edge_cross_dxdy(4, edge)
                nl2 = mesh%nlevels(el(2)) - 1
                nu2 = mesh%ulevels(el(2))
            end if
            nl12 = min(nl1, nl2)
            nu12 = max(nu1, nu2)
            ! (A)
            do nz = nu1, nu12-1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                call flux_ho(nz, clof(enodes(1),nz), clof(enodes(2),nz), vflux)
            end do
            ! (B)
            if (nu2 > 0) then
                do nz = nu2, nu12-1
                    vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                    call flux_ho(nz, clof(enodes(1),nz), clof(enodes(2),nz), vflux)
                end do
            end if
            ! (C)
            do nz = nu12, nl12
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1)) &
                      + ( vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                call flux_ho(nz, clof(enodes(1),nz), clof(enodes(2),nz), vflux)
            end do
            ! (D)
            do nz = nl12+1, nl1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                call flux_ho(nz, clof(enodes(1),nz), clof(enodes(2),nz), vflux)
            end do
            ! (E)
            do nz = nl12+1, nl2
                vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                call flux_ho(nz, clof(enodes(1),nz), clof(enodes(2),nz), vflux)
            end do
        end do

    contains
        real(kind=WP) function clof(node, nz)
            integer, intent(in) :: node, nz
            clof = real(max(sign(1, nboundary_lay(node)-nz), 0), WP)
        end function clof
        subroutine flux_ho(nz, clo1, clo2, vflux)
            integer,       intent(in) :: nz
            real(kind=WP), intent(in) :: clo1, clo2, vflux
            real(kind=WP) :: Tmean1, Tmean2, cHO
            Tmean2 = ttf(nz, enodes(2)) - &
                     (2.0_WP*(ttf(nz, enodes(2))-ttf(nz, enodes(1))) + &
                      mesh%edge_dxdy(1,edge)*edge_up_dn_grad(2,nz,edge) + &
                      mesh%edge_dxdy(2,edge)*edge_up_dn_grad(4,nz,edge))/6.0_WP*clo2
            Tmean1 = ttf(nz, enodes(1)) + &
                     (2.0_WP*(ttf(nz, enodes(2))-ttf(nz, enodes(1))) + &
                      mesh%edge_dxdy(1,edge)*edge_up_dn_grad(1,nz,edge) + &
                      mesh%edge_dxdy(2,edge)*edge_up_dn_grad(3,nz,edge))/6.0_WP*clo1
            cHO = (vflux+abs(vflux))*Tmean1 + (vflux-abs(vflux))*Tmean2
            flux(nz,edge) = -0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(Tmean1+Tmean2) - flux(nz,edge)
        end subroutine flux_ho
    end subroutine oracle_muscl

    subroutine oracle_mfct(vel, ttf, mesh, num_ord, flux, edge_up_dn_grad, o_init_zero, partit)
        ! MUSCL for the FCT path (oce_adv_tra_hor.F90:546-834). Same as
        ! adv_tra_hor_muscl but WITHOUT the c_lo bottom-boundary clamp (the
        ! reconstruction near bottom topography is not upwind; runs with FCT only).
        ! M2.12b: optional partit -> loop over OWNED edges (myDim_edge2D).
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: num_ord
        real(kind=WP), intent(in)    :: ttf(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: vel(2, mesh%nl-1, mesh%elem2D)
        real(kind=WP), intent(inout) :: flux(mesh%nl-1, mesh%edge2D)
        real(kind=WP), intent(in)    :: edge_up_dn_grad(4, mesh%nl-1, mesh%edge2D)
        logical, optional, intent(in) :: o_init_zero
        type(t_partit), intent(in), optional :: partit
        logical       :: l_init_zero
        real(kind=WP) :: deltaX1, deltaY1, deltaX2, deltaY2, vflux
        integer       :: el(2), enodes(2), nz, edge, nu12, nl12, nl1, nl2, nu1, nu2
        integer       :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        l_init_zero = .true.
        if (present(o_init_zero)) l_init_zero = o_init_zero
        if (l_init_zero) then
            do edge = 1, nEdgeO
                flux(:, edge) = 0.0_WP
            end do
        end if

        do edge = 1, nEdgeO
            enodes = mesh%edges(:, edge)
            el     = mesh%edge_tri(:, edge)
            nl1    = mesh%nlevels(el(1)) - 1
            nu1    = mesh%ulevels(el(1))
            deltaX1 = mesh%edge_cross_dxdy(1, edge)
            deltaY1 = mesh%edge_cross_dxdy(2, edge)
            nl2 = 0; nu2 = 0
            if (el(2) > 0) then
                deltaX2 = mesh%edge_cross_dxdy(3, edge)
                deltaY2 = mesh%edge_cross_dxdy(4, edge)
                nl2 = mesh%nlevels(el(2)) - 1
                nu2 = mesh%ulevels(el(2))
            end if
            nl12 = min(nl1, nl2)
            nu12 = max(nu1, nu2)
            ! (A)
            do nz = nu1, nu12-1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                call flux_ho(nz, vflux)
            end do
            ! (B)
            if (nu2 > 0) then
                do nz = nu2, nu12-1
                    vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                    call flux_ho(nz, vflux)
                end do
            end if
            ! (C)
            do nz = nu12, nl12
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1)) &
                      + ( vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                call flux_ho(nz, vflux)
            end do
            ! (D)
            do nz = nl12+1, nl1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                call flux_ho(nz, vflux)
            end do
            ! (E)
            do nz = nl12+1, nl2
                vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                call flux_ho(nz, vflux)
            end do
        end do

    contains
        subroutine flux_ho(nz, vflux)
            integer,       intent(in) :: nz
            real(kind=WP), intent(in) :: vflux
            real(kind=WP) :: Tmean1, Tmean2, cHO
            Tmean2 = ttf(nz, enodes(2)) - &
                     (2.0_WP*(ttf(nz, enodes(2))-ttf(nz, enodes(1))) + &
                      mesh%edge_dxdy(1,edge)*edge_up_dn_grad(2,nz,edge) + &
                      mesh%edge_dxdy(2,edge)*edge_up_dn_grad(4,nz,edge))/6.0_WP
            Tmean1 = ttf(nz, enodes(1)) + &
                     (2.0_WP*(ttf(nz, enodes(2))-ttf(nz, enodes(1))) + &
                      mesh%edge_dxdy(1,edge)*edge_up_dn_grad(1,nz,edge) + &
                      mesh%edge_dxdy(2,edge)*edge_up_dn_grad(3,nz,edge))/6.0_WP
            cHO = (vflux+abs(vflux))*Tmean1 + (vflux-abs(vflux))*Tmean2
            flux(nz,edge) = -0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(Tmean1+Tmean2) - flux(nz,edge)
        end subroutine flux_ho
    end subroutine oracle_mfct

end module muscl_oracle

program test_muscl_onthefly
    ! MUSCL horizontal advection with on-the-fly up/downwind gradients
    ! (docs/plans/2026-10-06-muscl-onthefly.md).
    !
    ! WHY THIS TEST EXISTS
    ! --------------------
    ! FESOM2's stored per-edge array edge_up_dn_grad(4, nl-1, nEdgeO), built by
    ! fill_up_dn_grad, is replaced in FESOM3 by an on-the-fly lookup: the elemental gradient tr_xy of
    ! the upwind/downwind triangle on the levels shared by both edge nodes, the Miura
    ! node-averaged gradient gnod (muscl_node_grad, once per node) on the levels of one node
    ! only and on edges without both triangles, and 0 where the fill never writes. The
    ! refactor must be BIT-IDENTICAL; fill_up_dn_grad (kept in oce_muscl_adv for this
    ! purpose only, writing a caller-owned array) and the pre-refactor kernels (module
    ! muscl_oracle above) are the oracle.
    !
    ! THE LOOKUP RULE (side k = 1: upwind / edge node 1, components 1,3 of the old array;
    ! k = 2: downwind / edge node 2, components 2,4), with
    ! nzmin = maxval(ulevels_nod2D_max(edge nodes)), nzmax = minval(nlevels_nod2D_min(...)):
    !   both up/dn triangles exist .and. nzmin <= nz < nzmax -> tr_xy(:,nz,edge_up_dn_tri(k))
    !   else ulevels_nod2D(nk) <= nz < nlevels_nod2D(nk)      -> gnod(:,nz,nk)
    !   else                                                  -> 0
    ! (with both triangles the fill's node-only ranges [ulev(nk), nzmin) and
    ! [nzmax, nlev(nk)) are exactly the complement of the shared range inside the node's
    ! range, because nzmin >= ulev(nk) and nzmax <= nlev(nk); an empty shared range makes
    ! the two node-only ranges overlap, which the fill resolves to the same average).
    !
    ! PART G - the rule reproduces the whole old array, every owned edge, every level
    !   nz = 1..nl-1, both sides, both components, BITWISE, on
    !   G1 the pi mesh as read (ulevels = 1), and
    !   G2 pi with synthesised cavity columns (ulevels_nod2D = 3 on every 10th global node
    !      with >= 14 levels; element ulevels and ulevels_nod2D_max rebuilt consistently,
    !      exchanged at np 2), which populates the upper node-only ranges [ulev(nk), nzmin)
    !      that do not exist without cavities.
    !   Positions where the fill stores 0/0 (a level of a node with no wet element around
    !   it) are compared as "oracle NaN <=> gnod 0" and counted; the kernels never read
    !   them (part F, Task 2, pins that through the fluxes).
    !   Every lookup class (triangle / node average / zero) must be populated, and in G2
    !   the upper node-only range must occur (teeth).
    !   Positive control: the reference RK3 code's rule (0 instead of the node average)
    !   must produce mismatches.
    !
    ! PART F - the fluxes: the pre-refactor kernels oracle_muscl / oracle_mfct with the
    !   stored array vs the production adv_tra_hor_muscl / adv_tra_hor_mfct with
    !   (tr_xy, gnod, edge_up_dn_tri),
    !   every owned edge and level, both signs of the velocity, BITWISE, on both meshes
    !   of part G. A nonzero edge velocity on every wet level and a tracer with structure
    !   in x, y and z make every reconstruction branch contribute. Positive control: the
    !   on-the-fly kernels with gnod = 0 (the reference RK3 rule) must differ.
    use mpi
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use mod_precision,     only: WP, MP
    use mod_mesh,          only: t_mesh
    use mod_partit,        only: t_partit
    use mod_tracer,        only: t_tracer_work
    use mod_partitioning,  only: par_init, par_ex, set_partition
    use mod_mesh_read,     only: read_mesh
    use mod_mesh_areas,    only: compute_geometry
    use mod_part_bounds,   only: owned_bounds, is_multirank
    use mod_halo,          only: exchange_elem_full, exchange_nod
    use oce_tracer_grad,   only: tracer_gradient_elements
    use oce_muscl_adv,     only: muscl_adv_init, fill_up_dn_grad, muscl_node_grad
    use oce_adv_tra_hor,   only: adv_tra_hor_muscl, adv_tra_hor_mfct
    use muscl_oracle,      only: oracle_muscl, oracle_mfct
    implicit none

    character(len=512) :: mesh_dir
    type(t_partit)      :: partit
    type(t_mesh)        :: mesh
    type(t_tracer_work) :: twork
    integer :: nfail, nsw
    integer :: nNodO, nNodL, nEdgeO, nElemO, nElemF, nl
    real(kind=WP), allocatable :: ttf(:,:), tr_xy(:,:,:), gnod(:,:,:)
    real(kind=WP), allocatable :: vel(:,:,:), eudg(:,:,:), f_old(:,:), f_new(:,:), gzero(:,:,:)
    real(kind=WP), parameter   :: num_ord = 0.25_WP      ! any mix of the 3rd/4th-order parts
    ! cavity synthesis (the constants of test_wsplit / test_wimpl_tra: same column set)
    integer, parameter :: cav_every    = 10      ! every 10th global node ...
    integer, parameter :: cav_ulev     = 3       ! ... gets its surface at level 3 ...
    integer, parameter :: cav_nlev_min = 14      ! ... if it has >= 14 levels (>= 11 layers below)

    nfail = 0
    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call par_init(partit)
    if (partit%npes > 1) call set_partition(partit, trim(mesh_dir))
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
    nElemF = mesh%elem2D
    if (partit%npes > 1) nElemF = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
    nl = mesh%nl
    allocate(ttf(nl-1, nNodL), tr_xy(2, nl-1, nElemF), gnod(2, nl-1, nNodL))
    call muscl_adv_init(twork, mesh, partit)      ! edge_up_dn_tri, nboundary_lay
    allocate(vel(2, nl-1, nElemF), eudg(4, nl-1, nEdgeO), f_old(nl-1, nEdgeO), f_new(nl-1, nEdgeO))
    allocate(gzero(2, nl-1, nNodL)); gzero = 0.0_WP
    allocate(mesh%helem(nl-1, nElemF))
    call build_flow()

    call part_g('G1 pi as read      ')
    call part_f('F1 pi as read      ')
    call synth_cavity()
    call part_g('G2 pi with cavities')
    call part_f('F2 pi with cavities')

    if (partit%mype == 0) then
        if (nfail == 0) then
            write(*,'(a)') 'test_muscl_onthefly: OK'
        else
            write(*,'(a,i0,a)') 'test_muscl_onthefly: ', nfail, ' FAILURE(S)'
        end if
    end if
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)   ! NO third arg: abort is presence-based
    if (nfail /= 0) error stop 1

contains
    subroutine check_true(name, cond)
        character(len=*), intent(in) :: name
        logical,          intent(in) :: cond
        if (.not. cond) then
            nfail = nfail + 1
            if (partit%mype == 0) write(*,'(a)') '  FAIL: '//name
        end if
    end subroutine check_true

    function gsum(i) result(s)
        integer, intent(in) :: i
        integer :: s, ierr
        s = i
        if (partit%npes > 1) call MPI_Allreduce(i, s, 1, MPI_INTEGER, MPI_SUM, &
                                                partit%MPI_COMM_FESOM, ierr)
    end function gsum

    !=========================================================================
    ! tracer field and its elemental gradient (full halo), as init_tracers_AB builds it
    !=========================================================================
    subroutine build_gradient()
        ! smooth large-scale structure + a level- and node-dependent wiggle, all functions
        ! of the node COORDINATES (identical on owner and halo copies)
        integer :: n, nz
        real(kind=WP) :: x, y
        do n = 1, nNodL
            x = mesh%coord_nod2D(1, n); y = mesh%coord_nod2D(2, n)
            do nz = 1, nl-1
                ttf(nz, n) = 10.0_WP*cos(y)*sin(x + 0.1_WP*real(nz, WP)) &
                           + 0.37_WP*sin(17.0_WP*x*y + 1.3_WP*real(nz, WP)) - 0.05_WP*real(nz, WP)
            end do
        end do
        tr_xy = 0.0_WP
        call tracer_gradient_elements(ttf, tr_xy, mesh, partit)
        if (is_multirank(partit)) call exchange_elem_full(tr_xy, partit)
    end subroutine build_gradient

    !=========================================================================
    ! PART G
    !=========================================================================
    subroutine part_g(label)
        character(len=*), intent(in) :: label
        integer :: edge, nz, k, nk, nzmin, nzmax, comp, ednodes(2)
        integer :: n_tri, n_nod, n_zero, n_nan, n_bad, n_upper, n_ctl
        logical :: both
        real(kind=WP) :: g(2), old(2)

        call build_gradient()
        call fill_up_dn_grad(eudg, twork, tr_xy, mesh, partit)    ! the oracle (zeroes eudg first)
        call muscl_node_grad(gnod, tr_xy, mesh, partit)

        n_tri = 0; n_nod = 0; n_zero = 0; n_nan = 0; n_bad = 0; n_upper = 0; n_ctl = 0
        do edge = 1, nEdgeO
            ednodes = mesh%edges(:, edge)
            both  = (twork%edge_up_dn_tri(1, edge) /= 0) .and. (twork%edge_up_dn_tri(2, edge) /= 0)
            nzmin = maxval(mesh%ulevels_nod2D_max(ednodes))
            nzmax = minval(mesh%nlevels_nod2D_min(ednodes))
            do k = 1, 2
                nk = ednodes(k)
                do nz = 1, nl-1
                    ! the old array: side k stores x in component k, y in component k+2
                    old(1) = eudg(k,   nz, edge)
                    old(2) = eudg(k+2, nz, edge)
                    ! the lookup rule (header)
                    if (both .and. nz >= nzmin .and. nz < nzmax) then
                        g = tr_xy(:, nz, twork%edge_up_dn_tri(k, edge))
                        n_tri = n_tri + 1
                    else if (nz >= mesh%ulevels_nod2D(nk) .and. nz < mesh%nlevels_nod2D(nk)) then
                        g = gnod(:, nz, nk)
                        n_nod = n_nod + 1
                        if (both .and. nz < nzmin) n_upper = n_upper + 1
                        ! positive control: the reference's 0 instead of the node average
                        if (any(old /= 0.0_WP .and. .not. ieee_is_nan(old))) n_ctl = n_ctl + 1
                    else
                        g = 0.0_WP
                        n_zero = n_zero + 1
                    end if
                    do comp = 1, 2
                        if (ieee_is_nan(old(comp))) then
                            ! fill's 0/0 at a level with no wet element: gnod stores 0 there
                            n_nan = n_nan + 1
                            if (g(comp) /= 0.0_WP) n_bad = n_bad + 1
                        else if (transfer(g(comp), 0_8) /= transfer(old(comp), 0_8)) then
                            n_bad = n_bad + 1
                        end if
                    end do
                end do
            end do
        end do
        n_tri = gsum(n_tri); n_nod = gsum(n_nod); n_zero = gsum(n_zero); n_nan = gsum(n_nan)
        n_bad = gsum(n_bad); n_upper = gsum(n_upper); n_ctl = gsum(n_ctl)
        if (partit%mype == 0) then
            write(*,'(2a,i0,a,i0,a,i0,a,i0)') '  ', label//': entries from triangle ', n_tri, &
                '  node average ', n_nod, '  zero ', n_zero, '  (of which upper node-only ', n_upper
            write(*,'(a,i0,a,i0,a,i0)') '      mismatches (bitwise) ', n_bad, &
                '  oracle 0/0 positions ', n_nan, '  positive control (0 for node avg) differs at ', n_ctl
        end if
        call check_true(trim(label)//': rule == fill_up_dn_grad bitwise', n_bad == 0)
        call check_true(trim(label)//': every lookup class populated', &
                        n_tri > 0 .and. n_nod > 0 .and. n_zero > 0)
        call check_true(trim(label)//': positive control (reference zero) differs', n_ctl > 0)
        if (index(label, 'cavities') > 0) &
            call check_true(trim(label)//': upper node-only range populated', n_upper > 0)
    end subroutine part_g

    !=========================================================================
    ! PART F
    !=========================================================================
    subroutine build_flow()
        ! element velocity with structure in both components and every level; helem any
        ! positive thickness (the kernels only multiply the edge velocity by it)
        integer :: e, nz
        real(kind=WP) :: a
        do e = 1, nElemF
            a = 0.731_WP*real(mod(e, 97), WP)
            do nz = 1, nl-1
                vel(1, nz, e) = 0.3_WP*sin(a + 0.2_WP*real(nz, WP)) + 0.05_WP
                vel(2, nz, e) = 0.2_WP*cos(1.7_WP*a - 0.1_WP*real(nz, WP))
                mesh%helem(nz, e) = 10.0_MP + 2.0_MP*real(nz, MP)
            end do
        end do
    end subroutine build_flow

    subroutine part_f(label)
        character(len=*), intent(in) :: label
        integer :: isign, n_bad_m, n_bad_f, n_ctl, n_flux
        ! part_g left tr_xy, gnod and the oracle array eudg of this mesh in place
        n_bad_m = 0; n_bad_f = 0; n_ctl = 0; n_flux = 0
        do isign = 1, 2
            if (isign == 2) vel = -vel
            call oracle_muscl(vel, ttf, mesh, num_ord, f_old, eudg, twork%nboundary_lay, &
                                   o_init_zero=.true., partit=partit)
            call adv_tra_hor_muscl(vel, ttf, mesh, num_ord, f_new, tr_xy, gnod, twork%edge_up_dn_tri, &
                                       twork%nboundary_lay, o_init_zero=.true., partit=partit)
            n_bad_m = n_bad_m + count(transfer(f_old, 0_8, size(f_old)) /= transfer(f_new, 0_8, size(f_new)))
            n_flux  = n_flux  + count(f_old /= 0.0_WP)
            call oracle_mfct(vel, ttf, mesh, num_ord, f_old, eudg, o_init_zero=.true., partit=partit)
            call adv_tra_hor_mfct(vel, ttf, mesh, num_ord, f_new, tr_xy, gnod, twork%edge_up_dn_tri, &
                                      o_init_zero=.true., partit=partit)
            n_bad_f = n_bad_f + count(transfer(f_old, 0_8, size(f_old)) /= transfer(f_new, 0_8, size(f_new)))
            ! positive control: the reference rule (no node average)
            call adv_tra_hor_mfct(vel, ttf, mesh, num_ord, f_new, tr_xy, gzero, twork%edge_up_dn_tri, &
                                      o_init_zero=.true., partit=partit)
            n_ctl = n_ctl + count(f_old /= f_new)
        end do
        vel = -vel                                   ! restore
        n_bad_m = gsum(n_bad_m); n_bad_f = gsum(n_bad_f); n_ctl = gsum(n_ctl); n_flux = gsum(n_flux)
        if (partit%mype == 0) write(*,'(2a,i0,a,i0,a,i0,a,i0)') '  ', label//': nonzero fluxes ', n_flux, &
            '  bitwise mismatches muscl ', n_bad_m, '  mfct ', n_bad_f, '  positive control differs at ', n_ctl
        call check_true(trim(label)//': adv_tra_hor_muscl == pre-refactor kernel bitwise', n_bad_m == 0)
        call check_true(trim(label)//': adv_tra_hor_mfct == pre-refactor kernel bitwise', n_bad_f == 0)
        call check_true(trim(label)//': fluxes nonzero (not vacuous)', n_flux > 0)
        call check_true(trim(label)//': positive control (gnod = 0) differs', n_ctl > 0)
    end subroutine part_f

    !=========================================================================
    ! Cavity columns, element-consistent
    !=========================================================================
    subroutine synth_cavity()
        ! node ulevels by the GLOBAL index (owner and halo copies agree); element ulevels
        ! = maxval over its vertices (owned elements, then the full element halo);
        ! ulevels_nod2D = minval and ulevels_nod2D_max = maxval over the elements around
        ! the node (owned, then exchange_nod) -- FESOM2's cavity convention.
        integer :: n, g, e, k, ncav
        ncav = 0
        do n = 1, nNodL
            g = n
            if (partit%npes > 1) g = partit%myList_nod2D(n)
            if (mod(g, cav_every) == 0 .and. mesh%nlevels_nod2D(n) >= cav_nlev_min) then
                mesh%ulevels_nod2D(n) = cav_ulev
                if (n <= nNodO) ncav = ncav + 1
            end if
        end do
        do e = 1, nElemO
            mesh%ulevels(e) = maxval(mesh%ulevels_nod2D(mesh%elem2D_nodes(1:3, e)))
        end do
        if (is_multirank(partit)) call exchange_elem_full(mesh%ulevels, partit)
        ! FESOM2's cavity convention: a node's surface is the SHALLOWEST surface of its
        ! elements (ulevels_nod2D = minval(ulevels around)), so a node whose elements all
        ! touch a cavity vertex moves down too; without this its top levels would have no
        ! wet element (the fill would then divide 0/0)
        do n = 1, nNodO
            k = mesh%nod_in_elem2D_num(n)
            mesh%ulevels_nod2D(n)     = minval(mesh%ulevels(mesh%nod_in_elem2D(1:k, n)))
            mesh%ulevels_nod2D_max(n) = maxval(mesh%ulevels(mesh%nod_in_elem2D(1:k, n)))
        end do
        if (is_multirank(partit)) then
            call exchange_nod(mesh%ulevels_nod2D, partit)
            call exchange_nod(mesh%ulevels_nod2D_max, partit)
        end if
        ncav = gsum(ncav)
        if (partit%mype == 0) write(*,'(a,i0)') '  synthesised owned cavity columns (ulevels_nod2D = 3): ', ncav
        call check_true('cavity columns synthesised', ncav > 0)
        do e = 1, nElemO
            call check_e(mesh%nlevels(e)-1 >= mesh%ulevels(e))
        end do
    end subroutine synth_cavity

    subroutine check_e(ok)
        logical, intent(in) :: ok
        if (.not. ok) call check_true('synthesised element keeps a wet level', .false.)
    end subroutine check_e
end program test_muscl_onthefly
