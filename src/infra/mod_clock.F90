module mod_clock
    ! Model calendar clock — state + per-step advance + cold-start init + run-length
    ! -> total-steps mapping. Transcribed VERBATIM from FESOM2 v2.7.3
    ! src/gen_modules_clock.F90 (module g_clock) and the get_run_steps routine in
    ! src/gen_model_setup.F90:334-382. The arithmetic (the >86400 / >ndpyr rollovers,
    ! the num_day_in_month month scan, the cold-start old==new test) is byte-faithful so a
    ! long run stays max|Δ|=0 vs FESOM2 by induction.
    !
    ! Deviations from g_clock, all intentional and recorded in the M8 plan:
    !   * use mod_config (not g_config) for dt / include_fleapyear / runid / RestartInPath /
    !     run_length / run_length_unit / step_per_day.
    !   * r_restart lives HERE (g_clock read it from g_config); the output-file creation that
    !     consumes it is M9, so mod_config stays edit-free.
    !   * clock_finish / clock_newyear (the .clock WRITE on restart) and the use_transit lines
    !     are OMITTED — both are M9 / out of scope.
    use mod_precision, only: WP
    use mod_config,    only: dt, include_fleapyear, runid, RestartInPath, &
                             run_length, run_length_unit, step_per_day
    use mod_partit,    only: t_partit
    use, intrinsic :: iso_fortran_env, only: error_unit
    use mpi
    implicit none
    public
    save

    real(kind=WP)  :: timeold, timenew     ! time in a day, unit: sec
    integer        :: dayold, daynew       ! day in a year
    integer        :: yearold, yearnew     ! year before and after time step
    integer        :: yearstart            ! year when simulation started
    integer        :: month, day_in_month  ! month and day in a month
    integer        :: fleapyear            ! 1 fleapyear, 0 not
    integer        :: ndpyr                ! number of days in yearnew
    integer        :: num_day_in_month(0:1,12)
    character(4)   :: cyearold, cyearnew   ! year as character string
    character(2)   :: cmonth               ! month as character string
    logical        :: r_restart = .false.  ! NEW home (g_clock read this from g_config)
    data num_day_in_month(0,:) /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
    data num_day_in_month(1,:) /31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

contains
    !
    !--------------------------------------------------------------------------------
    !
    subroutine clock
        ! Advance the clock one timestep: time -> day -> year -> month. Called at the TOP
        ! of each step (mirror FESOM2 fesom_module.F90:673). g_clock::clock VERBATIM.
        implicit none
        integer          :: i
        real(kind=WP)    :: aux1, aux2
        !
        timeold=timenew
        dayold=daynew
        yearold=yearnew

        ! update time
        timenew=timenew+dt

        ! update day
        if (timenew>86400._WP) then  !assumed that time step is less than one day!
           daynew=daynew+1
           timenew=timenew-86400._WP
        endif

        ! update year
        if (daynew>ndpyr) then
           daynew=1
           yearnew=yearnew+1
           call check_fleapyr(yearnew, fleapyear)
           ndpyr=365+fleapyear
           write(cyearold,'(i4)') yearold
           write(cyearnew,'(i4)') yearnew
        endif

        ! find month and dayinmonth at new time step
        aux1=0
        do i=1,12
           aux2=aux1+num_day_in_month(fleapyear,i)
           if(daynew>aux1 .and. daynew<=aux2) then
              month=i
              write(cmonth, '(I2.2)') month
              day_in_month=daynew-aux1
              exit
           end if
           aux1=aux2
        end do

    end subroutine clock
    !
    !--------------------------------------------------------------------------------
    !
    subroutine clock_init(partit)
        ! Initialise the clock for this run from RestartInPath//runid//'.clock' (2 lines:
        ! old, new). Cold start = both lines equal (r_restart=.false., yearold=yearnew-1);
        ! restart = lines differ (r_restart=.true.). g_clock::clock_init VERBATIM minus the
        ! use_transit lines.
        implicit none
        type(t_partit), intent(in), target    :: partit
        integer                               :: i, daystart
        real(kind=WP)                         :: aux1, aux2, timestart
        integer                               :: ierr
        integer                               :: file_unit
        character(512)                        :: errmsg

        ! init clock for this run - read clock file FIRST
        open(newunit=file_unit, file=trim(RestartInPath)//trim(runid)//'.clock', action='read', &
            status='old', iostat=ierr, iomsg=errmsg)
        if (ierr /= 0) then
          write (unit=error_unit, fmt='(3A)') &
            '### error: can not open file ', trim(RestartInPath)//trim(runid)//'.clock', &
            ', error: ' // trim(errmsg)
          call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        read(unit=file_unit, fmt=*) timeold, dayold, yearold
        read(unit=file_unit, fmt=*) timenew, daynew, yearnew
        close(unit=file_unit)
        if(daynew==0) daynew=1

        ! the model initialized at - set AFTER reading clock file
        timestart=timenew
        daystart=daynew
        yearstart=yearnew

        ! check if this is a restart or not
        ! For initial run: clock file has same values on both lines (old=new)
        ! For restart: clock file has different old/new values
        if(yearnew==yearold .and. daynew==dayold .and. timenew==timeold) then
           r_restart=.false.
           yearold=yearnew-1 !required for checking if create new output files
        else
           r_restart=.true.
        end if

        ! year as character string
        write(cyearold,'(i4)') yearold
        write(cyearnew,'(i4)') yearnew

        ! if restart model at beginning of a day, set timenew to be zero
        if (timenew==86400._WP) then
           timenew=0.0_WP
           daynew=daynew+1
        endif

        ! check fleap year
        call check_fleapyr(yearnew, fleapyear)
        ndpyr=365+fleapyear

        ! find month and dayinmonth at the new time step
        aux1=0
        do i=1,12
           aux2=aux1+num_day_in_month(fleapyear,i)
           if(daynew>aux1 .and. daynew<=aux2) then
              month=i
              day_in_month=daynew-aux1
              exit
           end if
           aux1=aux2
        end do

        if(partit%mype==0) then
            if(r_restart) then
                write(*,*)
                print *, achar(27)//'[31m'    //'____________________________________________________________'//achar(27)//'[0m'
                print *, achar(27)//'[5;7;31m'//' --> THIS IS A RESTART RUN !!!                              '//achar(27)//'[0m'
                write(*,"(A, F8.2, I4, I5)") '     > clock restarted at time:', timenew, daynew, yearnew
                write(*,"(A, I5)") '     > yearstart for annual_event:', yearstart
                write(*,*)
            else
                write(*,*)
                print *, achar(27)//'[32m'  //'____________________________________________________________'//achar(27)//'[0m'
                print *, achar(27)//'[7;32m'//' --> THIS IS A INITIALISATION RUN !!!                       '//achar(27)//'[0m'
                write(*,"(A, F8.2, I4, I5)")'     > clock initialized at time:', timenew, daynew, yearnew
                write(*,"(A, I5)") '     > yearstart for annual_event:', yearstart
                write(*,*)
            end if
        end if

    end subroutine clock_init
    !
    !-------------------------------------------------------------------------------
    !
    integer function clock_nsteps(partit) result(nsteps)
        ! Map run_length + run_length_unit -> total step count. Transcribed VERBATIM from
        ! get_run_steps (gen_model_setup.F90:334-382): 's'->run_length; 'd'->step_per_day*
        ! run_length; 'm'->per-month sum of step_per_day*num_day_in_month over the spanned
        ! months; 'y'->per-year sum of step_per_day*(365+fleapyear). The 'm'/'y' loops
        ! re-evaluate check_fleapyr each period so the M9 1948-2009 leap-correct roadmap holds
        ! (for the noleap M8 target 'y' reduces to run_length*365*step_per_day). clock must be
        ! initialised (clock_init) before this is called.
        implicit none
        type(t_partit), intent(in) :: partit
        integer :: i, temp_year, temp_mon, temp_fleapyear, ierr

        if(run_length_unit=='s') then
           nsteps=run_length
        elseif(run_length_unit=='d') then
           nsteps=step_per_day*run_length
        elseif(run_length_unit=='m') then
           nsteps=0
           temp_mon=month-1
           temp_year=yearnew
           temp_fleapyear=fleapyear
           do i=1,run_length
              temp_mon=temp_mon+1
              if(temp_mon>12) then
                 temp_year=temp_year+1
                 temp_mon=1
                 call check_fleapyr(temp_year, temp_fleapyear)
              end if
              nsteps=nsteps+step_per_day*num_day_in_month(temp_fleapyear,temp_mon)
           end do
        elseif(run_length_unit=='y') then
           nsteps=0
           do i=1,run_length
              temp_year=yearnew+i-1
              call check_fleapyr(temp_year, temp_fleapyear)
              nsteps=nsteps+step_per_day*(365+temp_fleapyear)
           end do
        else
           write(*,*) 'Run length unit ', run_length_unit, ' is not defined.'
           write(*,*) 'Please check and update the code.'
           call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        if(partit%mype==0) write(*,*) nsteps, ' steps to run for ', trim(runid), ' job submission'
    end function clock_nsteps
    !
    !----------------------------------------------------------------------------
    !
    subroutine check_fleapyr(year, flag)
        ! flag=1 if `year` is a leap year AND leap years are enabled, else 0. With
        ! include_fleapyear=.false. (the CORE2 noleap target) this always returns 0.
        implicit none
        integer, intent(in) :: year
        integer, intent(out):: flag

        flag=0

        if(.not.include_fleapyear) return
        call is_fleapyr(year, flag)
    end subroutine check_fleapyr

    subroutine is_fleapyr(year, flag)
        ! The Gregorian leap-year rule (dead for the CORE2 noleap target).
        implicit none
        integer, intent(in) :: year
        integer, intent(out):: flag
        flag=0
        if ((mod(year,4)==0.and.mod(year,100)/=0) .or. mod(year,400)==0) then
           flag=1
        endif
    end subroutine is_fleapyr
    !
    !----------------------------------------------------------------------------
    !
end module mod_clock
