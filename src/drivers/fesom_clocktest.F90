program fesom_clocktest
    ! M8a-1 unit-level GATE for mod_clock: advance `clock` from a cold-start 0 1 1948 and
    ! verify timenew/daynew/yearnew/month/day_in_month against a HAND-COMPUTED table (literals
    ! below, NOT recomputed by the same code), then check clock_nsteps for s/d/m/y units.
    ! Proves the calendar arithmetic (day/month/year rollovers, the >86400/>ndpyr triggers,
    ! the num_day_in_month month scan, get_run_steps) before it drives forcing in M8b/c.
    !
    ! Anchor: step_per_day=48 (dt=1800), include_fleapyear=.false. (CORE2 noleap). Hand table
    ! (each milestone = the clock state AFTER k `clock` calls; day = exactly 48 steps so the
    ! >86400 strict test keeps step 48 on day 1 and rolls on step 49):
    !   k=48     86400.0  day  1  1948  month  1  dim  1      (last step of day 1)
    !   k=49      1800.0  day  2  1948  month  1  dim  2      (day rollover)
    !   k=1488   86400.0  day 31  1948  month  1  dim 31      (last step of Jan)
    !   k=1489    1800.0  day 32  1948  month  2  dim  1      (MONTH rollover Jan->Feb)
    !   k=17520  86400.0  day365  1948  month 12  dim 31      (last step of the year)
    !   k=17521   1800.0  day  1  1949  month  1  dim  1      (YEAR rollover, day 366->1)
    !
    !   FESOM3_RESTART_IN  dir for the test1.clock cold-start file (default ./)
    use mpi
    use mod_precision,    only: WP
    use mod_config,       only: dt, step_per_day, include_fleapyear, runid, RestartInPath, &
                                run_length, run_length_unit
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex
    use mod_clock,        only: clock, clock_init, clock_nsteps, &
                                timenew, daynew, yearnew, month, day_in_month, ndpyr, fleapyear, r_restart
    implicit none

    type(t_partit) :: partit
    character(512) :: rin
    integer        :: k, file_unit, ierr, nfail

    call par_init(partit)

    ! --- anchor config (the CORE2 noleap timestep) ---
    step_per_day      = 48
    dt                = 86400.0_WP / real(step_per_day, WP)   ! = 1800 s; clock reads mod_config%dt
    include_fleapyear = .false.
    runid             = 'test1'
    call get_environment_variable('FESOM3_RESTART_IN', rin)
    if (len_trim(rin) == 0) rin = './'
    RestartInPath = trim(rin)

    nfail = 0

    ! --- cold-start clock file: both lines equal => r_restart=.false. ---
    if (partit%mype == 0) then
        open(newunit=file_unit, file=trim(RestartInPath)//trim(runid)//'.clock', &
             status='replace', action='write')
        write(file_unit,*) 0.0_WP, 1, 1948
        write(file_unit,*) 0.0_WP, 1, 1948
        close(file_unit)
    end if
    call MPI_Barrier(partit%MPI_COMM_FESOM, ierr)

    call clock_init(partit)

    ! --- init state: 0 1 1948, month 1, dim 1, ndpyr 365, cold start ---
    if (partit%mype == 0) then
        write(*,'(a)') '--- clock_init (cold start) ---'
        call check_r('timenew',      timenew,      0.0_WP)
        call check_i('daynew',       daynew,       1)
        call check_i('yearnew',      yearnew,      1948)
        call check_i('month',        month,        1)
        call check_i('day_in_month', day_in_month, 1)
        call check_i('ndpyr',        ndpyr,        365)
        call check_i('fleapyear',    fleapyear,    0)
        if (r_restart) then
            write(*,'(a)') '  FAIL r_restart: got T want F'; nfail = nfail + 1
        end if

        ! --- clock_nsteps (get_run_steps) BEFORE advancing (reads yearnew=1948, month=1) ---
        write(*,'(a)') '--- clock_nsteps ---'
        run_length_unit = 's'; run_length = 5;   call check_i('nsteps s,5',   clock_nsteps(partit), 5)
        run_length_unit = 'd'; run_length = 730; call check_i('nsteps d,730', clock_nsteps(partit), 35040)
        run_length_unit = 'm'; run_length = 1;   call check_i('nsteps m,1',   clock_nsteps(partit), 1488)
        run_length_unit = 'm'; run_length = 2;   call check_i('nsteps m,2',   clock_nsteps(partit), 2832)
        run_length_unit = 'y'; run_length = 1;   call check_i('nsteps y,1',   clock_nsteps(partit), 17520)
        run_length_unit = 'y'; run_length = 2;   call check_i('nsteps y,2',   clock_nsteps(partit), 35040)

        write(*,'(a)') '--- clock advance (milestones) ---'
    end if

    ! --- advance and check the hand table at each milestone ---
    do k = 1, 17521
        call clock
        if (partit%mype /= 0) cycle
        select case (k)
        case (48);    call milestone(k, 86400.0_WP, 1,   1948, 1,  1)   ! last step of day 1
        case (49);    call milestone(k,  1800.0_WP, 2,   1948, 1,  2)
        case (1488);  call milestone(k, 86400.0_WP, 31,  1948, 1,  31)
        case (1489);  call milestone(k,  1800.0_WP, 32,  1948, 2,  1)
        case (17520); call milestone(k, 86400.0_WP, 365, 1948, 12, 31)
        case (17521); call milestone(k,  1800.0_WP, 1,   1949, 1,  1)
        end select
    end do

    if (partit%mype == 0) then
        write(*,*)
        if (nfail == 0) then
            write(*,'(a)') 'fesom_clocktest: PASS (all milestones + clock_nsteps max|Δ|=0)'
        else
            write(*,'(a,i0,a)') 'fesom_clocktest: FAIL (', nfail, ' mismatches)'
        end if
    end if

    call MPI_Bcast(nfail, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierr)
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    if (nfail /= 0) error stop 1

contains

    subroutine check_i(label, got, want)
        character(len=*), intent(in) :: label
        integer,          intent(in) :: got, want
        if (got == want) then
            write(*,'(a,a,a,i0)') '  ok   ', label, ' = ', got
        else
            write(*,'(a,a,a,i0,a,i0)') '  FAIL ', label, ' got ', got, ' want ', want
            nfail = nfail + 1
        end if
    end subroutine check_i

    subroutine check_r(label, got, want)
        character(len=*), intent(in) :: label
        real(kind=WP),    intent(in) :: got, want
        if (got == want) then
            write(*,'(a,a,a,es14.6)') '  ok   ', label, ' = ', got
        else
            write(*,'(a,a,a,es14.6,a,es14.6)') '  FAIL ', label, ' got ', got, ' want ', want
            nfail = nfail + 1
        end if
    end subroutine check_r

    subroutine milestone(k, t, d, y, mo, dim)
        integer,       intent(in) :: k, d, y, mo, dim
        real(kind=WP), intent(in) :: t
        character(len=32) :: tag
        write(tag,'(a,i0,a)') 'k=', k, ' '
        call check_r(trim(tag)//'timenew',      timenew,      t)
        call check_i(trim(tag)//'daynew',       daynew,       d)
        call check_i(trim(tag)//'yearnew',      yearnew,      y)
        call check_i(trim(tag)//'month',        month,        mo)
        call check_i(trim(tag)//'day_in_month', day_in_month, dim)
    end subroutine milestone

end program fesom_clocktest
