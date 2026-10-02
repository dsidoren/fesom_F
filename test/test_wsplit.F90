program test_wsplit
    ! Smooth Courant-number-dependent explicit/implicit vertical-velocity split
    ! (oce_wsplit, docs/plans/2026-10-02-wsplit-smooth.md).
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
    ! comments.
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
    ! WHAT IS ASSERTED, for the parameter sets (0.5,1), (0.9,1), (0,1), (0.25,0.5):
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
    ! Part W needs no mesh and no communication; it runs identically on every rank.
    use mpi
    use mod_precision,    only: WP
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex
    use oce_wsplit,       only: wsplit_implicit_fraction, wsplit_check_params
    implicit none

    type(t_partit) :: partit
    integer :: nfail, iset
    integer,       parameter :: nset  = 4
    integer,       parameter :: ngrid = 10000
    real(kind=WP), parameter :: cu_top = 20.0_WP
    real(kind=WP), parameter :: pmin(nset) = [0.5_WP, 0.9_WP, 0.0_WP, 0.25_WP]
    real(kind=WP), parameter :: pmax(nset) = [1.0_WP, 1.0_WP, 1.0_WP, 0.5_WP]
    real(kind=WP), parameter :: eps = epsilon(1.0_WP)

    nfail = 0
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
end program test_wsplit
