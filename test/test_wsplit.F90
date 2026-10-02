program test_wsplit
    ! Smooth Courant-number-dependent explicit/implicit vertical-velocity split
    ! (oce_wsplit + compute_Wvel_split, docs/plans/2026-10-02-wsplit-smooth.md).
    !
    ! WHY THIS TEST EXISTS
    ! --------------------
    ! use_wsplit splits w = w_e + w_i and advects w_i implicitly. FESOM2 (oce_ale.F90:
    ! 3217-3265) does this with a hard switch at CFL_z = wsplit_maxcfl, so d(w_e)/d(CFL)
    ! jumps from 1 to 0 there and a face oscillating around the threshold flips between
    ! high-order explicit and partly first-order implicit treatment. FESOM3 replaces the
    ! switch by the C^1 limiting function of Shchepetkin (2015, Ocean Modelling 91, 38-69,
    ! Sec. 3.1 / Fig. 9; NEMO sshwzv.F90 wAimp). There is no oracle for it: FESOM2 has only
    ! the hard switch. So this test IS the specification (L54): every property the plan
    ! claims for the function is asserted here, with the closed forms derived in the
    ! comments. Part W needs no mesh; parts S and X run on the pi mesh at np 1/2.
    !
    ! THE FUNCTION. f = wsplit_implicit_fraction(Cu, Cu_min, Cu_max) is the IMPLICIT share,
    ! w_i = f*w, w_e = w - w_i. With D = Cu_max - Cu_min, F = 4*Cu_max*D and
    ! Cu_cut = 2*Cu_max - Cu_min:
    !     Cu <= Cu_min           f = 0                                   (fully explicit)
    !     Cu_min < Cu < Cu_cut   f = x^2/(F + x^2),   x = Cu - Cu_min   (smooth bend)
    !     Cu >= Cu_cut           f = (Cu - Cu_max)/Cu                   (Cu_e == Cu_max)
    ! Closed forms used below: at Cu_cut, x = 2D, so the middle branch gives
    ! 4D^2/(4*Cu_max*D + 4D^2) = D/Cu_cut and the top branch (Cu_cut - Cu_max)/Cu_cut =
    ! D/Cu_cut (C^0). Slopes: middle 2xF/(F+x^2)^2 = 16*Cu_max*D^2/(16*D^2*Cu_cut^2) =
    ! Cu_max/Cu_cut^2, top Cu_max/Cu^2 = Cu_max/Cu_cut^2 (C^1). At Cu_min the middle branch
    ! has f = 0 and f' = 2xF/(F+x^2)^2 = 0, matching the flat lower branch (C^1). f -> 1
    ! as Cu -> inf; Cu_e = Cu*(1-f) is non-decreasing (middle branch: d/dCu proportional to
    ! F - x^2 - 2x*Cu_min >= 0 for x <= 2D) and equals Cu_max exactly on the top branch.
    !
    ! FLOATING POINT. On the top branch 1-f = Cu_max/Cu is formed by cancellation, so
    ! Cu*(1-f) carries an ABSOLUTE error ~eps*Cu (relative ~eps*Cu/Cu_max); W2 and W6
    ! account for it.
    !
    ! PART W - the function, for the parameter sets (0.5,1), (0.9,1), (0,1), (0.25,0.5):
    !   W1 f == 0 exactly for Cu in {0, Cu_min/2, Cu_min}
    !   W2 cap: |Cu*(1-f) - Cu_max| <= 1e-14*Cu for Cu in {Cu_cut, 2, 5, 50, 1e4}
    !   W3 f(1e6*Cu_max) > 1 - 2e-6; 0 <= f <= 1 on a 10^4-point grid in [0, 20]
    !   W4 C^0 joints: the two branch formulas agree AT Cu_min and AT Cu_cut to 1e-14, and
    !      the module reproduces each branch formula on its own branch (grid, 4 eps)
    !   W5 C^1 joints, h = 1e-6: at Cu_min both one-sided difference quotients have
    !      |q| <= 2h/F (closed-form slope 0: the right quotient is h/(F+h^2)); at Cu_cut
    !      both agree with Cu_max/Cu_cut^2 to 1e-5 relative
    !   W6 monotone: f non-decreasing on the grid; Cu_e non-decreasing up to the
    !      saturated-branch noise, 4*eps*Cu (see the derivation at the check)
    !   W7 elemental: the array call equals the element-wise loop bitwise
    !   W8 wsplit_check_params accepts (0,1), (0.5,1), (1,1); rejects (-0.1,1), (1.1,1),
    !      (0.5,0)
    !
    ! PART S - compute_Wvel_split on the pi mesh, every OWNED+HALO face (the routine
    ! produces the halo, which the next step's momentum advection reads). Prescribed w
    ! (nonzero at every level incl. the surface, both signs) and a cfl_z field spanning
    ! [0, 5] in steps of 0.05, so it contains faces exactly AT 0, Cu_min, Cu_max and Cu_cut;
    ! parameter sets (0.5,1), (0.9,1), (0.25,0.5):
    !   S1 |w_e + w_i - w| <= 1 ulp(w) at every face (w_e = w - w_i; the consumers need no
    !      more than that)
    !   S2 w_i == f(cfl_z)*w BITWISE: the split IS the limiting function, nothing else
    !   S3 use_wsplit=.false.: w_e == w and w_i == +0 bitwise (sign of zero included) --
    !      the production drivers' path, byte-identical to the previous off path
    !   S4 cfl_z <= Cu_min -> w_i == 0 exactly; cfl_z >= Cu_cut -> |w_e|*cfl_z/|w| ==
    !      Cu_max to 1e-14*cfl_z (the explicit Courant number is capped EXACTLY, FESOM2's
    !      meaning of wsplit_maxcfl is kept); every class non-empty
    !
    ! PART X - transition smoothness on the NON-FCT tracer path: QR4C explicit vertical
    ! advection on w_e (adv_tra_ver_qr4c -> oce_tra_adv_flux2dtracer -> the ALE reconstruct
    ! of diff_tracers_ale) followed by the upwind TDMA on w_i (diff_ver_part_impl_ale,
    ! do_wimpl), Kv = 0, no surface fluxes, Redi off, tra_adv_lim /= 'FCT'.
    !
    ! WHY THIS PATH AND NOT THE FCT LOW-ORDER PAIR (plan review). Explicit upwind on w_e
    ! followed by implicit upwind on w_i across one face is EXACTLY split-independent: the
    ! donor cell gives T2*(h2 - dt*w_e)/(h2 - dt*w + dt*w_i) = T2 and the receiver then holds
    ! T1*h1 + dt*w*T2 whatever f, so a hard switch leaves no first-order kink there. The
    ! kink is first order only where the explicit and implicit operators differ at first
    ! order: QR4C (centred, 4th order) against the upwind TDMA. Hence this scan.
    !
    ! X0 CLOSED FORM OF THE KINK. Uniform layers h (hnode = hnode_new = h: a smoothness
    ! probe, not a conservation test), uniform interior w (faces 2..nzmax-1; the surface and
    ! bottom faces 0, so compute_CFLz gives one interior Courant number Cu = 2|w|dt/h --
    ! it counts the face flux against BOTH adjacent cells), a quadratic T(z) = T0 + bq*z +
    ! cq*z^2 sampled at the cell mid-depths. Per unit w the per-step increments at a cell p
    ! whose faces p and p+1 lie in QR4C's 4th-order range (nz in [nzmin+2, nzmax-2]) are
    !   QR4C (num_ord = 1): face value T(z_face) - cq*h^2/12 at BOTH faces (the h^2 term
    !                       cancels in the divergence), so
    !                       P_q T(p) = -dt*(T(z_p + h/2) - T(z_p - h/2))/h = -dt*T'(z_p)
    !   upwind, w > 0:      face value = the cell below -> P_u T(p) = -dt*(T'(z_p) - cq*h)
    !   upwind, w < 0:      face value = the cell above -> P_u T(p) = -dt*(T'(z_p) + cq*h)
    ! The hard switch (Cu_min = Cu_max = C, threshold velocity w_c) gives
    !   below C:  Delta = w*P_q T,                        dDelta/dw = P_q T
    !   above C:  w_e = w_c frozen, w_i = w - w_c,  Delta = (I + w_i*U)^-1 T* - T, where
    !             T* = T + w_c*P_q T is the tracer AFTER the explicit step (which has already
    !             acted when the TDMA runs) and (I + w_i*U)^-1 is the TDMA with U the upwind
    !             matrix per unit velocity; d/dw_i at w_i = 0 is -U = P_u, so
    !             dDelta/dw = P_u T* = P_u T + w_c*P_u(P_q T)
    ! P_q T is linear in z with slope -2*cq*dt, and P_u of a linear profile is -dt times
    ! its slope = 2*cq*dt^2 for either sign of w. The jump of dDelta/dw at the threshold is
    !   J_w = (P_u - P_q)T + 2*cq*dt^2*w_c = sgn(w)*cq*dt*h*(1 + Cu_c),   Cu_c = 2|w_c|dt/h
    ! and in units of the scan variable Cu (dw/dCu = sgn(w)*h/(2dt))
    !   J = cq*h^2*(1 + Cu_c)/2                              (both signs of w)
    ! The plan's shorthand J = [UPW - QR4C](T) is the Cu_c -> 0 limit cq*h^2/2; the factor
    ! (1 + Cu_c) is the pre-advection of T to T* at the threshold and X2 pins it to 5 %.
    ! A function with a derivative jump J has the second difference J*dCu at a grid point
    ! ON the kink (the scan puts Cu = Cu_max exactly on a grid point; the explicit branch
    ! is exactly linear in Cu, the implicit branch adds O(dCu/(1+Cu_c)) curvature) and
    ! J*dCu split over two neighbours otherwise -- hence the plan's >= 0.5*|J|*dCu floor.
    ! bq is chosen so that T'(z_p) = 4*cq*h at the probe: the kink then exceeds the X1
    ! bound by ~8x for the hard switch while the smooth function stays ~20x below it.
    !
    !   X1 smooth (0.5,1), 601 points in [0, 3] (dCu = 0.005), probe cell p = 5 of every
    !      owned column with >= 10 levels (its stencil, and the T* cells the closed form
    !      needs, stay in the 4th-order range), both signs of w: the first differences D1
    !      of Delta(Cu) are continuous, i.e. the second differences obey
    !      |D2| <= 10*dCu*max|D1| everywhere
    !   X2 hard switch (1,1), the same scan: VIOLATES the X1 criterion (proves X1 has
    !      teeth), and at the kink D2 >= 0.5*|J|*dCu and |D2 - J*dCu| <= 0.05*|J|*dCu
    !      (pins J, both signs of w)
    !   X3 limits: for Cu <= Cu_min the split chain equals the explicit QR4C chain with the
    !      full w BITWISE at every owned cell; w_e = 0, w_i = w gives bitwise the pure
    !      do_wimpl TDMA result, and that TDMA equals the closed-form recursion of its rows
    !      (r = |w|dt/h; donor-side (S + r*S_donor)/(1 + r), the open surface/bottom rows as
    !      diff_ver_part_impl_ale builds them) to 1e-13 relative, both signs
    use mpi
    use mod_precision,    only: WP, MP
    use mod_mesh,         only: t_mesh
    use mod_dyn,          only: t_dyn
    use mod_tracer,       only: t_tracer
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex, set_partition
    use mod_mesh_read,    only: read_mesh
    use mod_mesh_areas,   only: compute_geometry
    use mod_part_bounds,  only: owned_bounds
    use oce_wsplit,       only: wsplit_implicit_fraction, wsplit_check_params
    use oce_ale,          only: compute_CFLz, compute_Wvel_split
    use oce_adv_tra_ver,  only: adv_tra_ver_qr4c
    use oce_adv_tra_flux, only: oce_tra_adv_flux2dtracer
    use oce_ale_tracer,   only: diff_ver_part_impl_ale
    implicit none

    character(len=512) :: mesh_dir
    type(t_partit) :: partit
    type(t_mesh)   :: mesh
    type(t_dyn)    :: dyn
    type(t_tracer) :: tracers
    integer :: nfail, iset, nsw
    integer :: nNodO, nNodL, nEdgeO, nElemO, nl
    ! ---- part W ----
    integer,       parameter :: nset  = 4
    integer,       parameter :: ngrid = 10000
    real(kind=WP), parameter :: cu_top = 20.0_WP
    real(kind=WP), parameter :: pmin(nset) = [0.5_WP, 0.9_WP, 0.0_WP, 0.25_WP]
    real(kind=WP), parameter :: pmax(nset) = [1.0_WP, 1.0_WP, 1.0_WP, 0.5_WP]
    real(kind=WP), parameter :: eps = epsilon(1.0_WP)
    ! ---- part S ----
    integer,       parameter :: nsset = 3
    real(kind=WP), parameter :: smin(nsset) = [0.5_WP, 0.9_WP, 0.25_WP]
    real(kind=WP), parameter :: smax(nsset) = [1.0_WP, 1.0_WP, 0.5_WP]
    ! ---- part X ----
    real(kind=WP), parameter :: dt   = 1800.0_WP     ! s
    real(kind=WP), parameter :: h0   = 10.0_WP       ! uniform layer thickness [m]
    real(kind=WP), parameter :: cq   = 1.0e-3_WP     ! quadratic coefficient of T(z)
    real(kind=WP), parameter :: t0   = 10.0_WP
    integer,       parameter :: np   = 5             ! probe cell
    integer,       parameter :: nlev_min = 10        ! columns probed: nlevels_nod2D >= this
    integer,       parameter :: nscan = 601
    real(kind=WP), parameter :: x_maxcfl = 1.0_WP, x_mincfl = 0.5_WP
    real(kind=WP), parameter :: cu_top_x = 3.0_WP*x_maxcfl
    real(kind=WP), parameter :: dcu = cu_top_x/real(nscan-1, WP)
    ! T'(z_p) = bq + 2*cq*z_p = 4*cq*h0 at z_p = -(np-0.5)*h0
    real(kind=WP), parameter :: bq = cq*h0*(2.0_WP*real(np, WP) + 3.0_WP)
    real(kind=WP), parameter :: jkink = cq*h0*h0*(1.0_WP + x_maxcfl)/2.0_WP   ! X0
    real(kind=WP), allocatable :: tref(:,:), tsplit(:,:), tother(:,:)
    real(kind=WP), allocatable :: flux_v(:,:), flux_h(:,:), dttf_h(:,:), dttf_v(:,:)
    real(kind=WP), allocatable :: heat_flux(:), water_flux(:), virtual_salt(:), relax_salt(:), real_salt_flux(:)
    real(kind=WP), allocatable :: delta(:,:)

    nfail = 0
    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call par_init(partit)

    !=========================================================================
    ! Part W - the limiting function, closed forms
    !=========================================================================
    do iset = 1, nset
        call part_w(pmin(iset), pmax(iset))
    end do

    !=========================================================================
    ! W8 - parameter validation
    !=========================================================================
    call check_true('W8 accepts (0,1)',    wsplit_check_params(0.0_WP, 1.0_WP))
    call check_true('W8 accepts (0.5,1)',  wsplit_check_params(0.5_WP, 1.0_WP))
    call check_true('W8 accepts (1,1)',    wsplit_check_params(1.0_WP, 1.0_WP))
    call check_true('W8 rejects (-0.1,1)', .not. wsplit_check_params(-0.1_WP, 1.0_WP))
    call check_true('W8 rejects (1.1,1)',  .not. wsplit_check_params(1.1_WP, 1.0_WP))
    call check_true('W8 rejects (0.5,0)',  .not. wsplit_check_params(0.5_WP, 0.0_WP))

    !=========================================================================
    ! pi mesh (the scaffold of test_vinv / test_ivertvisc)
    !=========================================================================
    if (partit%npes > 1) call set_partition(partit, trim(mesh_dir))
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.true.)
    call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
    nl = mesh%nl
    allocate(dyn%w(nl, nNodL), dyn%w_e(nl, nNodL), dyn%w_i(nl, nNodL), dyn%cfl_z(nl, nNodL))
    dyn%w = 0.0_WP; dyn%w_e = 0.0_WP; dyn%w_i = 0.0_WP; dyn%cfl_z = 0.0_WP

    !=========================================================================
    ! Part S - compute_Wvel_split identities on owned+halo faces
    !=========================================================================
    do iset = 1, nsset
        call part_s(smin(iset), smax(iset))
    end do

    !=========================================================================
    ! Part X - transition smoothness on the non-FCT tracer path
    !=========================================================================
    call setup_x()
    call part_x()

    if (partit%mype == 0) then
        if (nfail == 0) then
            write(*,'(a)') 'test_wsplit: OK'
        else
            write(*,'(a,i0,a)') 'test_wsplit: ', nfail, ' FAILURE(S)'
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

    ! --- MPI reductions (identity at np 1) ------------------------------------------
    function gmax(x) result(m)
        real(kind=WP), intent(in) :: x
        real(kind=WP) :: m
        integer :: ierr
        m = x
        if (partit%npes > 1) call MPI_Allreduce(x, m, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                                                partit%MPI_COMM_FESOM, ierr)
    end function gmax

    function gmin(x) result(m)
        real(kind=WP), intent(in) :: x
        real(kind=WP) :: m
        integer :: ierr
        m = x
        if (partit%npes > 1) call MPI_Allreduce(x, m, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                                                partit%MPI_COMM_FESOM, ierr)
    end function gmin

    function gsum(i) result(s)
        integer, intent(in) :: i
        integer :: s, ierr
        s = i
        if (partit%npes > 1) call MPI_Allreduce(i, s, 1, MPI_INTEGER, MPI_SUM, &
                                                partit%MPI_COMM_FESOM, ierr)
    end function gsum

    function gall(l) result(a)
        logical, intent(in) :: l
        logical :: a
        integer :: ierr
        a = l
        if (partit%npes > 1) call MPI_Allreduce(l, a, 1, MPI_LOGICAL, MPI_LAND, &
                                                partit%MPI_COMM_FESOM, ierr)
    end function gall

    ! --- the branch formulas of the specification (NOT the module's code) -----------
    pure function f_mid(cu, cmin, cmax) result(f)
        real(kind=WP), intent(in) :: cu, cmin, cmax
        real(kind=WP) :: f, x, ff
        x  = cu - cmin
        ff = 4.0_WP*cmax*(cmax - cmin)
        f  = x*x/(ff + x*x)
    end function f_mid

    pure function f_top(cu, cmax) result(f)
        real(kind=WP), intent(in) :: cu, cmax
        real(kind=WP) :: f
        f = (cu - cmax)/cu
    end function f_top

    subroutine part_w(cmin, cmax)
        real(kind=WP), intent(in) :: cmin, cmax
        character(len=24) :: tag
        real(kind=WP) :: d, ff, ccut, slope, h, f, cue
        real(kind=WP) :: cap_err, ql, qr, q_bound, sl_errl, sl_errr, dec, dec_bound, fbig
        real(kind=WP) :: cug(ngrid), fg(ngrid), fl(ngrid), cueg(ngrid), fexp(ngrid)
        real(kind=WP) :: cu_cap(5)
        integer :: k
        logical :: ok

        d    = cmax - cmin
        ff   = 4.0_WP*cmax*d
        ccut = 2.0_WP*cmax - cmin
        write(tag,'(a,f5.2,a,f5.2,a)') '(', cmin, ',', cmax, ')'
        if (partit%mype == 0) write(*,'(a)') '  --- parameters '//trim(tag)//'  ---'

        ! ---- W1: fully explicit below and at Cu_min, exact zero ----------------------
        ok = wsplit_implicit_fraction(0.0_WP, cmin, cmax) == 0.0_WP
        ok = ok .and. wsplit_implicit_fraction(0.5_WP*cmin, cmin, cmax) == 0.0_WP
        ok = ok .and. wsplit_implicit_fraction(cmin, cmin, cmax) == 0.0_WP
        call check_true('W1 f == 0 exactly for Cu <= Cu_min '//tag, ok)

        ! ---- W2: the explicit Courant number is capped at Cu_max on the top branch -----
        cu_cap = [ccut, 2.0_WP, 5.0_WP, 50.0_WP, 1.0e4_WP]
        cap_err = 0.0_WP
        do k = 1, 5
            f   = wsplit_implicit_fraction(cu_cap(k), cmin, cmax)
            cue = cu_cap(k)*(1.0_WP - f)
            cap_err = max(cap_err, abs(cue - cmax)/cu_cap(k))
        end do
        if (partit%mype == 0) write(*,'(a,es10.2)') '  W2 max |Cu_e - Cu_max|/Cu on the cap = ', cap_err
        call check_true('W2 cap |Cu*(1-f) - Cu_max| <= 1e-14*Cu '//tag, cap_err <= 1.0e-14_WP)

        ! ---- W3: limit and range -------------------------------------------------------
        fbig = wsplit_implicit_fraction(1.0e6_WP*cmax, cmin, cmax)
        if (partit%mype == 0) write(*,'(a,es10.2)') '  W3 1 - f(1e6*Cu_max)                   = ', 1.0_WP - fbig
        call check_true('W3 f(1e6*Cu_max) > 1 - 2e-6 '//tag, fbig > 1.0_WP - 2.0e-6_WP)
        do k = 1, ngrid
            cug(k) = cu_top*real(k-1, WP)/real(ngrid-1, WP)
        end do
        fg = wsplit_implicit_fraction(cug, cmin, cmax)         ! array (elemental) call
        call check_true('W3 0 <= f <= 1 on the grid '//tag, all(fg >= 0.0_WP) .and. all(fg <= 1.0_WP))

        ! ---- W4: C^0 joints -- branch formulas agree where the branches meet ----------
        ! At Cu_min the lower branch is 0 and the middle branch has x = 0 -> 0.
        ok = abs(f_mid(cmin, cmin, cmax) - 0.0_WP) <= 1.0e-14_WP
        call check_true('W4 branches agree at Cu_min '//tag, ok)
        ! At Cu_cut both branches give D/Cu_cut.
        ok = abs(f_mid(ccut, cmin, cmax) - f_top(ccut, cmax)) <= 1.0e-14_WP
        call check_true('W4 branches agree at Cu_cut '//tag, ok)
        ok = abs(f_mid(ccut, cmin, cmax) - d/ccut) <= 1.0e-14_WP .and. &
             abs(f_top(ccut, cmax)       - d/ccut) <= 1.0e-14_WP
        call check_true('W4 both branches equal D/Cu_cut at Cu_cut '//tag, ok)
        ! The module's value at the joints is the joint value.
        f = wsplit_implicit_fraction(ccut, cmin, cmax)
        ok = abs(f - d/ccut) <= 1.0e-14_WP
        call check_true('W4 module f(Cu_cut) == D/Cu_cut '//tag, ok)
        ! And the module reproduces each branch formula on its own branch.
        do k = 1, ngrid
            if (cug(k) <= cmin) then
                fexp(k) = 0.0_WP
            else if (cug(k) < ccut) then
                fexp(k) = f_mid(cug(k), cmin, cmax)
            else
                fexp(k) = f_top(cug(k), cmax)
            end if
        end do
        cap_err = maxval(abs(fg - fexp))
        if (partit%mype == 0) write(*,'(a,es10.2)') '  W4 max |f - branch formula| on the grid = ', cap_err
        call check_true('W4 module matches the branch formulas (4 eps) '//tag, cap_err <= 4.0_WP*eps)

        ! ---- W5: C^1 joints, one-sided difference quotients --------------------------
        h = 1.0e-6_WP
        ql = (wsplit_implicit_fraction(cmin, cmin, cmax) - wsplit_implicit_fraction(cmin - h, cmin, cmax))/h
        qr = (wsplit_implicit_fraction(cmin + h, cmin, cmax) - wsplit_implicit_fraction(cmin, cmin, cmax))/h
        q_bound = 2.0_WP*h/ff
        if (partit%mype == 0) write(*,'(a,2es10.2,a,es10.2)') '  W5 at Cu_min: q_left, q_right = ', ql, qr, &
                                                             '  bound 2h/F = ', q_bound
        call check_true('W5 slope 0 at Cu_min (both sides) '//tag, abs(ql) <= q_bound .and. abs(qr) <= q_bound)
        slope = cmax/(ccut*ccut)
        ql = (wsplit_implicit_fraction(ccut, cmin, cmax) - wsplit_implicit_fraction(ccut - h, cmin, cmax))/h
        qr = (wsplit_implicit_fraction(ccut + h, cmin, cmax) - wsplit_implicit_fraction(ccut, cmin, cmax))/h
        sl_errl = abs(ql - slope)/slope
        sl_errr = abs(qr - slope)/slope
        if (partit%mype == 0) write(*,'(a,es10.3,a,2es10.2)') '  W5 at Cu_cut: slope Cu_max/Cu_cut^2 = ', slope, &
                                                              '  rel err left/right = ', sl_errl, sl_errr
        call check_true('W5 slope Cu_max/Cu_cut^2 at Cu_cut (both sides, 1e-5) '//tag, &
                        sl_errl <= 1.0e-5_WP .and. sl_errr <= 1.0e-5_WP)

        ! ---- W6: monotone f and monotone explicit Courant number -----------------------
        ok = .true.
        do k = 2, ngrid
            ok = ok .and. (fg(k) >= fg(k-1))
        end do
        call check_true('W6 f non-decreasing on the grid '//tag, ok)
        cueg = cug*(1.0_WP - fg)
        ! Cu_e is exactly Cu_max on the saturated branch up to the cancellation noise of
        ! 1-f: f = (Cu-Cu_max)/Cu has an ABSOLUTE error <= eps (two roundings), 1-f is
        ! then exact (Sterbenz), and Cu*(1-f) inherits eps*Cu plus eps/2*Cu_max from the
        ! product. A decrease between neighbours of up to 2*eps*Cu + eps*Cu_max <= 4*eps*Cu
        ! is therefore round-off, not a monotonicity defect. Measured: 3.9e-15 at Cu ~ 17.5
        ! for every parameter set -- it scales with Cu, NOT with Cu_max (17.5 ulp of
        ! Cu_max = 1 but 35.5 ulp of Cu_max = 0.5), so the bound is in units of eps*Cu.
        dec = 0.0_WP
        ok  = .true.
        do k = 2, ngrid
            dec_bound = 4.0_WP*eps*max(cug(k), cmax)
            if (cueg(k-1) - cueg(k) > dec_bound) ok = .false.
            dec = max(dec, (cueg(k-1) - cueg(k))/(eps*max(cug(k), cmax)))
        end do
        if (partit%mype == 0) write(*,'(a,f8.2,a)') '  W6 max decrease of Cu_e between grid neighbours = ', dec, &
                                                   ' eps*Cu  (bound 4)'
        call check_true('W6 Cu_e non-decreasing up to saturated-branch noise (4 eps*Cu) '//tag, ok)

        ! ---- W7: elemental -- the array call IS the scalar loop --------------------------
        do k = 1, ngrid
            fl(k) = wsplit_implicit_fraction(cug(k), cmin, cmax)
        end do
        call check_true('W7 elemental array call == scalar loop bitwise '//tag, all(fg == fl))
    end subroutine part_w

    !=========================================================================
    ! Part S
    !=========================================================================
    subroutine part_s(cmin, cmax)
        real(kind=WP), intent(in) :: cmin, cmax
        character(len=24) :: tag
        integer :: n, nz, nzmin, nzmax, n_lo, n_mid, n_hi, sgn
        real(kind=WP) :: ccut, f, err1, err4, cfl, w
        logical :: ok1, ok2, ok3, ok4

        ccut = 2.0_WP*cmax - cmin
        write(tag,'(a,f5.2,a,f5.2,a)') '(', cmin, ',', cmax, ')'
        if (partit%mype == 0) write(*,'(a)') '  --- part S, parameters '//trim(tag)//'  ---'

        ! w: nonzero at EVERY level (the bottom face included -- the routine must not care),
        ! both signs; cfl_z: 0, 0.05, ..., 5 (contains 0, Cu_min, Cu_max, Cu_cut exactly)
        do n = 1, nNodL
            do nz = 1, nl
                sgn = merge(1, -1, mod(nz + n, 2) == 0)
                dyn%w(nz, n)     = real(sgn, WP)*1.0e-3_WP*(1.0_WP + 0.03_WP*real(nz, WP)) &
                                   *(1.0_WP + 0.1_WP*real(mod(n, 7), WP))
                dyn%cfl_z(nz, n) = 5.0_WP*real(mod(7*nz + 13*n, 101), WP)/100.0_WP
            end do
        end do
        dyn%use_wsplit    = .true.
        dyn%wsplit_mincfl = cmin
        dyn%wsplit_maxcfl = cmax
        dyn%w_e = -999.0_WP; dyn%w_i = -999.0_WP          ! must be overwritten at every face
        call compute_Wvel_split(dyn, mesh, partit)

        ok1 = .true.; ok2 = .true.; ok4 = .true.
        err1 = 0.0_WP; err4 = 0.0_WP; n_lo = 0; n_mid = 0; n_hi = 0
        do n = 1, nNodL
            if (mesh%nlevels_nod2D(n) <= 0) cycle
            nzmin = mesh%ulevels_nod2D(n); nzmax = mesh%nlevels_nod2D(n)
            do nz = nzmin, nzmax
                w   = dyn%w(nz, n)
                cfl = dyn%cfl_z(nz, n)
                ! S1
                err1 = max(err1, abs(dyn%w_e(nz, n) + dyn%w_i(nz, n) - w)/spacing(w))
                if (abs(dyn%w_e(nz, n) + dyn%w_i(nz, n) - w) > spacing(w)) ok1 = .false.
                ! S2
                f = wsplit_implicit_fraction(cfl, cmin, cmax)
                if (dyn%w_i(nz, n) /= f*w) ok2 = .false.
                ! S4
                if (cfl <= cmin) then
                    n_lo = n_lo + 1
                    if (dyn%w_i(nz, n) /= 0.0_WP) ok4 = .false.
                else if (cfl >= ccut) then
                    n_hi = n_hi + 1
                    err4 = max(err4, abs(abs(dyn%w_e(nz, n))*cfl/abs(w) - cmax)/cfl)
                else
                    n_mid = n_mid + 1
                end if
            end do
        end do
        err1 = gmax(err1); err4 = gmax(err4)
        n_lo = gsum(n_lo); n_mid = gsum(n_mid); n_hi = gsum(n_hi)
        if (partit%mype == 0) then
            write(*,'(a,f6.2,a)')          '  S1 max |w_e + w_i - w|                 = ', err1, ' ulp(w)'
            write(*,'(a,es10.2,a)')        '  S4 max |Cu_e - Cu_max|/cfl_z on the cap = ', err4, '  (bound 1e-14)'
            write(*,'(a,i0,a,i0,a,i0)')    '  S  faces below Cu_min / bend / capped   = ', n_lo, ' / ', n_mid, ' / ', n_hi
        end if
        call check_true('S1 |w_e + w_i - w| <= 1 ulp(w) at every owned+halo face '//tag, gall(ok1))
        call check_true('S2 w_i == f(cfl_z)*w bitwise '//tag, gall(ok2))
        call check_true('S4 w_i == 0 for cfl_z <= Cu_min '//tag, gall(ok4))
        call check_true('S4 |w_e|*cfl_z/|w| == Cu_max to 1e-14*cfl_z for cfl_z >= Cu_cut '//tag, err4 <= 1.0e-14_WP)
        call check_true('S  every cfl_z class populated '//tag, n_lo > 0 .and. n_mid > 0 .and. n_hi > 0)

        ! ---- S3: the off path ----------------------------------------------------------
        dyn%use_wsplit = .false.
        dyn%w_e = -999.0_WP; dyn%w_i = -999.0_WP
        call compute_Wvel_split(dyn, mesh, partit)
        ok3 = .true.
        do n = 1, nNodL
            if (mesh%nlevels_nod2D(n) <= 0) cycle
            nzmin = mesh%ulevels_nod2D(n); nzmax = mesh%nlevels_nod2D(n)
            do nz = nzmin, nzmax
                if (dyn%w_e(nz, n) /= dyn%w(nz, n)) ok3 = .false.
                if (dyn%w_i(nz, n) /= 0.0_WP .or. sign(1.0_WP, dyn%w_i(nz, n)) < 0.0_WP) ok3 = .false.
            end do
        end do
        call check_true('S3 use_wsplit=.false.: w_e == w, w_i == +0 bitwise '//tag, gall(ok3))
    end subroutine part_s

    !=========================================================================
    ! Part X
    !=========================================================================
    subroutine setup_x()
        ! uniform layers h0 (hnode = hnode_new), the depth arrays QR4C reads, a one-tracer
        ! state on the non-FCT path, Kv = 0, zero surface fluxes, zero horizontal fluxes.
        integer :: n, nz, nzmin, nzmax
        real(kind=WP) :: z
        allocate(mesh%hnode(nl-1, nNodL), mesh%hnode_new(nl-1, nNodL))
        allocate(mesh%zbar_3d_n(nl, nNodL), mesh%Z_3d_n(nl-1, nNodL))
        mesh%hnode = 0.0_MP; mesh%zbar_3d_n = 0.0_MP; mesh%Z_3d_n = 0.0_MP
        do n = 1, nNodL
            if (mesh%nlevels_nod2D(n) <= 0) cycle
            nzmin = mesh%ulevels_nod2D(n); nzmax = mesh%nlevels_nod2D(n)
            do nz = nzmin, nzmax - 1
                mesh%hnode(nz, n)  = real(h0, MP)
                mesh%Z_3d_n(nz, n) = real(-(real(nz, WP) - 0.5_WP)*h0, MP)
            end do
            do nz = nzmin, nzmax
                mesh%zbar_3d_n(nz, n) = real(-real(nz - 1, WP)*h0, MP)
            end do
        end do
        mesh%hnode_new = mesh%hnode

        tracers%num_tracers = 1
        allocate(tracers%data(1))
        allocate(tracers%data(1)%values(nl-1, nNodL))
        tracers%data(1)%values      = 0.0_WP
        tracers%data(1)%ID          = 1          ! temperature: bc_surface = -dt*(0 + T*0) = 0
        tracers%data(1)%tra_adv_ver = 'QR4C'
        tracers%data(1)%tra_adv_lim = 'NONE'     ! /= 'FCT' -> do_wimpl when use_wsplit
        tracers%data(1)%tra_adv_pv  = 1.0_WP     ! num_ord = 1: pure 4th-order centred
        tracers%data(1)%i_vert_diff = .true.
        allocate(dyn%work%Kv(nl, nNodL)); dyn%work%Kv = 0.0_WP
        allocate(heat_flux(nNodL), water_flux(nNodL), virtual_salt(nNodL), relax_salt(nNodL), real_salt_flux(nNodL))
        heat_flux = 0.0_WP; water_flux = 0.0_WP; virtual_salt = 0.0_WP; relax_salt = 0.0_WP; real_salt_flux = 0.0_WP
        allocate(flux_v(nl, nNodL), flux_h(nl-1, nEdgeO), dttf_h(nl-1, nNodL), dttf_v(nl-1, nNodL))
        flux_v = 0.0_WP; flux_h = 0.0_WP; dttf_h = 0.0_WP; dttf_v = 0.0_WP

        ! the quadratic profile at the cell mid-depths
        allocate(tref(nl-1, nNodL), tsplit(nl-1, nNodL), tother(nl-1, nNodL))
        tref = 0.0_WP
        do n = 1, nNodL
            if (mesh%nlevels_nod2D(n) <= 0) cycle
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n) - 1
                z = real(mesh%Z_3d_n(nz, n), WP)
                tref(nz, n) = t0 + bq*z + cq*z*z
            end do
        end do
        allocate(delta(nscan, nNodO))
    end subroutine setup_x

    subroutine set_uniform_w(w0)
        ! uniform interior w (faces nzmin+1..nzmax-1); surface and bottom faces 0
        real(kind=WP), intent(in) :: w0
        integer :: n, nz
        dyn%w = 0.0_WP
        do n = 1, nNodL
            if (mesh%nlevels_nod2D(n) <= 0) cycle
            do nz = mesh%ulevels_nod2D(n) + 1, mesh%nlevels_nod2D(n) - 1
                dyn%w(nz, n) = w0
            end do
        end do
    end subroutine set_uniform_w

    subroutine apply_chain(we, wi, with_impl, tin, tout)
        ! The non-FCT vertical sequence of one tracer step: QR4C on we -> flux divergence ->
        ! the ALE reconstruct of diff_tracers_ale (verbatim; the (hnode-hnode_new) term is 0
        ! here) -> (with_impl) the do_wimpl TDMA of diff_ver_part_impl_ale on wi.
        real(kind=WP), intent(in)  :: we(nl, nNodL), wi(nl, nNodL), tin(nl-1, nNodL)
        logical,       intent(in)  :: with_impl
        real(kind=WP), intent(out) :: tout(nl-1, nNodL)
        integer :: n, nzmin, nzmax
        tracers%data(1)%values = tin
        call adv_tra_ver_qr4c(we, tracers%data(1)%values, mesh, tracers%data(1)%tra_adv_pv, &
                              flux_v, o_init_zero=.true., partit=partit)
        dttf_h = 0.0_WP; dttf_v = 0.0_WP
        call oce_tra_adv_flux2dtracer(dt, dttf_h, dttf_v, flux_h, flux_v, mesh, partit=partit)
        do n = 1, nNodO
            nzmax = mesh%nlevels_nod2D(n) - 1
            nzmin = mesh%ulevels_nod2D(n)
            dttf_v(nzmin:nzmax, n) = dttf_v(nzmin:nzmax, n) + tracers%data(1)%values(nzmin:nzmax, n)* &
                                     (mesh%hnode(nzmin:nzmax, n) - mesh%hnode_new(nzmin:nzmax, n))
            tracers%data(1)%values(nzmin:nzmax, n) = tracers%data(1)%values(nzmin:nzmax, n) + &
                                     dttf_v(nzmin:nzmax, n)/mesh%hnode_new(nzmin:nzmax, n)
        end do
        if (with_impl) then
            dyn%w_i = wi
            dyn%use_wsplit = .true.                    ! do_wimpl
            call diff_ver_part_impl_ale(1, dt, dyn, tracers, mesh, heat_flux, water_flux, &
                                        virtual_salt, relax_salt, real_salt_flux, 0.0_WP, partit)
        end if
        tout = tracers%data(1)%values
    end subroutine apply_chain

    subroutine run_scan(cmin, cmax, wsign, ok_expl, n_expl)
        ! Delta(Cu) at the probe cell of every owned column, Cu on the scan grid; the
        ! velocity is derived from Cu and cfl_z is rebuilt by compute_CFLz (so the split
        ! sees exactly what the model would). X3a is checked on the fly: for Cu <= Cu_min the
        ! split chain must equal the explicit chain on the full w bitwise.
        real(kind=WP), intent(in)  :: cmin, cmax, wsign
        logical,       intent(out) :: ok_expl
        integer,       intent(out) :: n_expl
        integer :: k, n
        real(kind=WP) :: cu, w0
        dyn%use_wsplit    = .true.
        dyn%wsplit_mincfl = cmin
        dyn%wsplit_maxcfl = cmax
        ok_expl = .true.; n_expl = 0
        do k = 1, nscan
            cu = cu_top_x*real(k-1, WP)/real(nscan-1, WP)
            w0 = wsign*cu*h0/(2.0_WP*dt)
            call set_uniform_w(w0)
            call compute_CFLz(dyn, mesh, dt, partit)
            call compute_Wvel_split(dyn, mesh, partit)
            call apply_chain(dyn%w_e, dyn%w_i, .true., tref, tsplit)
            do n = 1, nNodO
                delta(k, n) = tsplit(np, n) - tref(np, n)
            end do
            if (cu <= cmin) then
                call apply_chain(dyn%w, dyn%w_i, .false., tref, tother)
                n_expl = n_expl + 1
                if (any(tsplit(:, 1:nNodO) /= tother(:, 1:nNodO))) ok_expl = .false.
            end if
        end do
    end subroutine run_scan

    subroutine analyse_scan(kk, ratio, d2k_min, d2k_max, nprobe)
        ! per probed column: D1(k) = Delta(k+1) - Delta(k), D2(k) = D1(k) - D1(k-1);
        ! ratio = max over columns of max|D2|/(dCu*max|D1|); D2 at the kink index kk.
        integer,       intent(in)  :: kk
        real(kind=WP), intent(out) :: ratio, d2k_min, d2k_max
        integer,       intent(out) :: nprobe
        integer :: n, k
        real(kind=WP) :: d1(nscan-1), d2max, d1max, d2k
        ratio = 0.0_WP; d2k_min = huge(1.0_WP); d2k_max = -huge(1.0_WP); nprobe = 0
        do n = 1, nNodO
            if (mesh%nlevels_nod2D(n) < nlev_min) cycle
            nprobe = nprobe + 1
            do k = 1, nscan - 1
                d1(k) = delta(k+1, n) - delta(k, n)
            end do
            d1max = maxval(abs(d1))
            d2max = 0.0_WP
            do k = 2, nscan - 1
                d2max = max(d2max, abs(d1(k) - d1(k-1)))
            end do
            ratio = max(ratio, d2max/(dcu*d1max))
            d2k = d1(kk) - d1(kk-1)
            d2k_min = min(d2k_min, d2k)
            d2k_max = max(d2k_max, d2k)
        end do
        ratio = gmax(ratio); d2k_min = gmin(d2k_min); d2k_max = gmax(d2k_max); nprobe = gsum(nprobe)
    end subroutine analyse_scan

    subroutine part_x()
        integer :: kk, isg, n, nz, nzmin, nzmax, nprobe, n_expl
        real(kind=WP) :: wsign, w0, r, ratio_s, ratio_h, d2k_min, d2k_max, d2k_err, err_cf, cu_kk
        real(kind=WP) :: s(nl-1)
        logical :: ok_expl, ok_bit
        character(len=8) :: stag

        kk = 1 + nint(x_maxcfl/dcu)                      ! the kink of the hard switch
        cu_kk = cu_top_x*real(kk-1, WP)/real(nscan-1, WP)
        call check_true('X0 the scan grid has a point exactly at Cu_max', cu_kk == x_maxcfl)
        if (partit%mype == 0) write(*,'(a,es10.3,a,i0,a)') '  X0 closed-form kink J = cq*h^2*(1+Cu_max)/2 = ', jkink, &
                                                          '   (kink at scan point ', kk, ')'

        do isg = 1, 2
            wsign = merge(1.0_WP, -1.0_WP, isg == 1)
            stag  = merge('(w > 0)', '(w < 0)', isg == 1)
            ! ---- X1: the smooth function ----------------------------------------------------
            call run_scan(x_mincfl, x_maxcfl, wsign, ok_expl, n_expl)
            call analyse_scan(kk, ratio_s, d2k_min, d2k_max, nprobe)
            if (partit%mype == 0) then
                write(*,'(a,a,a,f8.3,a,i0,a)') '  X1 smooth  (0.5,1) ', stag, ': max|D2|/(dCu*max|D1|) = ', ratio_s, &
                                              '  (bound 10; ', nprobe, ' columns)'
                write(*,'(a,a,a,f8.4,a)')      '  X1 smooth  (0.5,1) ', stag, ': D2(Cu_max)/(J*dCu)      = ', &
                                              d2k_max/(jkink*dcu), '  (no kink)'
            end if
            call check_true('X1 second differences <= 10*dCu*max|D1| (smooth) '//stag, ratio_s <= 10.0_WP .and. nprobe > 0)
            call check_true('X3a Cu <= Cu_min: split chain == explicit QR4C chain on the full w bitwise '//stag, &
                            gall(ok_expl) .and. n_expl > 0)

            ! ---- X2: the degenerate hard switch -- the positive control -------------------
            call run_scan(x_maxcfl, x_maxcfl, wsign, ok_expl, n_expl)
            call analyse_scan(kk, ratio_h, d2k_min, d2k_max, nprobe)
            d2k_err = max(abs(d2k_min - jkink*dcu), abs(d2k_max - jkink*dcu))/(jkink*dcu)
            if (partit%mype == 0) then
                write(*,'(a,a,a,f8.3,a)')      '  X2 hard    (1,1)   ', stag, ': max|D2|/(dCu*max|D1|) = ', ratio_h, &
                                              '  (must exceed 10)'
                write(*,'(a,a,a,f8.4,a,es9.2,a)') '  X2 hard    (1,1)   ', stag, ': D2(Cu_max)/(J*dCu)      = ', &
                                              d2k_min/(jkink*dcu), '  (min over columns; max rel err ', d2k_err, ')'
            end if
            call check_true('X2 hard switch violates the X1 criterion (X1 has teeth) '//stag, ratio_h > 10.0_WP)
            call check_true('X2 hard switch: D2 at the kink >= 0.5*|J|*dCu '//stag, d2k_min >= 0.5_WP*jkink*dcu)
            call check_true('X2 hard switch: D2 at the kink == J*dCu to 5% (pins J) '//stag, d2k_err <= 0.05_WP)

            ! ---- X3b: the fully implicit limit, w_e = 0, w_i = w at Cu = 2 -------------------
            w0 = wsign*2.0_WP*h0/(2.0_WP*dt)
            call set_uniform_w(w0)
            dyn%w_e = 0.0_WP
            call apply_chain(dyn%w_e, dyn%w, .true., tref, tsplit)     ! the chain
            tracers%data(1)%values = tref                               ! the TDMA alone
            dyn%w_i = dyn%w
            dyn%use_wsplit = .true.
            call diff_ver_part_impl_ale(1, dt, dyn, tracers, mesh, heat_flux, water_flux, &
                                        virtual_salt, relax_salt, real_salt_flux, 0.0_WP, partit)
            tother = tracers%data(1)%values
            ok_bit = all(tsplit(:, 1:nNodO) == tother(:, 1:nNodO))
            call check_true('X3b w_e = 0: chain == pure do_wimpl TDMA bitwise '//stag, gall(ok_bit))
            ! closed-form recursion of the upwind rows (r = |w|dt/h): the donor-side row
            ! (h + |w|dt)*S - |w|dt*S_donor = h*S_old; the open end rows h*S - |w|dt*S_donor
            ! = h*S_old (surface for w > 0, bottom for w < 0) and (h + |w|dt)*S = h*S_old
            ! (bottom for w > 0, surface for w < 0)
            r = abs(w0)*dt/h0
            err_cf = 0.0_WP
            do n = 1, nNodO
                if (mesh%nlevels_nod2D(n) <= 0) cycle
                nzmin = mesh%ulevels_nod2D(n); nzmax = mesh%nlevels_nod2D(n)
                s = 0.0_WP
                if (w0 > 0.0_WP) then
                    s(nzmax-1) = tref(nzmax-1, n)/(1.0_WP + r)
                    do nz = nzmax - 2, nzmin + 1, -1
                        s(nz) = (tref(nz, n) + r*s(nz+1))/(1.0_WP + r)
                    end do
                    s(nzmin) = tref(nzmin, n) + r*s(nzmin+1)
                else
                    s(nzmin) = tref(nzmin, n)/(1.0_WP + r)
                    do nz = nzmin + 1, nzmax - 2
                        s(nz) = (tref(nz, n) + r*s(nz-1))/(1.0_WP + r)
                    end do
                    s(nzmax-1) = tref(nzmax-1, n) + r*s(nzmax-2)
                end if
                err_cf = max(err_cf, maxval(abs(tother(nzmin:nzmax-1, n) - s(nzmin:nzmax-1))) &
                                     /maxval(abs(tref(nzmin:nzmax-1, n))))
            end do
            err_cf = gmax(err_cf)
            if (partit%mype == 0) write(*,'(a,a,a,es10.2)') '  X3b implicit limit ', stag, &
                                                          ': max rel err vs the row recursion = ', err_cf
            call check_true('X3b pure TDMA == closed-form upwind recursion (1e-13) '//stag, err_cf <= 1.0e-13_WP)
        end do
    end subroutine part_x
end program test_wsplit
