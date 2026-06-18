program test_params
    ! M0.2 gate: constants are bit-identical to the FESOM2 v2.7.3 values, and the
    ! write-once namelist readers accept a valid namelist and reject a malformed one.
    use mpi
    use mod_precision,  only: WP, MP, MPI_WP, MAX_PATH
    use mod_constants,  only: pi, rad, r_earth, omega, density_0, density_0_r, &
                              g, vcpw, inv_vcpw, small
    use mod_config,     only: read_config, runid, step_per_day, dt, which_ALE
    use mod_param_phys, only: read_param_phys, mix_scheme, Fer_GM, Redi, A_ver
    implicit none

    integer :: ierr, nfail
    nfail = 0
    call MPI_Init(ierr)

    ! ---- precision scaffolding (WP value is build-dependent) ----
#if defined(USE_SINGLE_PRECISION)
    call check_int('WP==4 (single)', WP, 4)
    call check_true('MPI_WP==MPI_REAL', MPI_WP == MPI_REAL)
#else
    call check_int('WP==8 (anchor)', WP, 8)
    call check_true('MPI_WP==MPI_DOUBLE_PRECISION', MPI_WP == MPI_DOUBLE_PRECISION)
#endif
    call check_int('MP==max(WP,4)',  MP, max(WP,4))
    call check_int('MAX_PATH',       MAX_PATH, 4096)

    ! ---- constants: bit-identical to FESOM2 v2.7.3 oce_modules.F90 ----
    call check_real('pi',          real(pi,WP),       3.14159265358979_WP)
    call check_real('rad',         real(rad,WP),      real(3.14159265358979_MP/180.0_MP, WP))
    call check_real('r_earth',     real(r_earth,WP),  6367500.0_WP)
    call check_real('omega',       real(omega,WP),    real(2*3.14159265358979_MP/(3600.0_MP*24.0_MP), WP))
    call check_real('density_0',   density_0,         1030.0_WP)
    call check_real('density_0_r', density_0_r,       1.0_WP/1030.0_WP)
    call check_real('g',           g,                 9.81_WP)
    call check_real('vcpw',        vcpw,              4.2e6_WP)
    call check_real('inv_vcpw',    inv_vcpw,          1.0_WP/4.2e6_WP)
    call check_real('small',       small,             1.0e-8_WP)

    ! ---- namelist accept ----
    call write_good_config()
    call read_config('test_params.config', ierr)
    call check_true('read_config accepts valid file', ierr == 0)
    call check_true('runid read',        trim(runid) == 'pitest')
    call check_int ('step_per_day read', step_per_day, 32)
    call check_real('dt derived',        dt, 86400.0_WP/32.0_WP)
    call check_true('which_ALE read',    trim(which_ALE) == 'linfs')

    call write_good_oce()
    call read_param_phys('test_params.oce', ierr)
    call check_true('read_param_phys accepts valid file', ierr == 0)
    call check_true('mix_scheme read', trim(mix_scheme) == 'PP')
    call check_true('Fer_GM read',     .not. Fer_GM)
    call check_true('Redi read',       .not. Redi)
    call check_real('A_ver read',      A_ver, 1.0e-4_WP)

    ! ---- namelist reject ----
    call write_bad_config()
    call read_config('test_params.bad', ierr)
    call check_true('read_config rejects malformed file', ierr > 0)

    call read_config('this_file_does_not_exist.config', ierr)
    call check_true('read_config rejects missing file', ierr > 0)

    ! ---- summary ----
    if (nfail == 0) then
        write(*,'(a)') 'test_params: ALL PASS'
    else
        write(*,'(a,i0,a)') 'test_params: ', nfail, ' FAILURE(S)'
    end if
    call MPI_Finalize(ierr)
    if (nfail /= 0) error stop 1

contains

    subroutine check_true(name, cond)
        character(len=*), intent(in) :: name
        logical,          intent(in) :: cond
        if (.not. cond) then
            nfail = nfail + 1
            write(*,'(a)') '  FAIL: '//name
        end if
    end subroutine

    subroutine check_int(name, got, want)
        character(len=*), intent(in) :: name
        integer,          intent(in) :: got, want
        if (got /= want) then
            nfail = nfail + 1
            write(*,'(a,i0,a,i0)') '  FAIL: '//name//' got=', got, ' want=', want
        end if
    end subroutine

    subroutine check_real(name, got, want)
        character(len=*), intent(in) :: name
        real(kind=WP),    intent(in) :: got, want
        if (got /= want) then   ! bit-exact equality intended
            nfail = nfail + 1
            write(*,'(a,es24.16,a,es24.16)') '  FAIL: '//name//' got=', got, ' want=', want
        end if
    end subroutine

    subroutine write_good_config()
        integer :: u
        open(newunit=u, file='test_params.config', status='replace', action='write')
        write(u,'(a)') '&modelname'
        write(u,'(a)') "  runid = 'pitest'"
        write(u,'(a)') '/'
        write(u,'(a)') '&timestep'
        write(u,'(a)') '  step_per_day = 32'
        write(u,'(a)') '/'
        write(u,'(a)') '&ale_def'
        write(u,'(a)') "  which_ALE = 'linfs'"
        write(u,'(a)') '/'
        close(u)
    end subroutine

    subroutine write_good_oce()
        integer :: u
        open(newunit=u, file='test_params.oce', status='replace', action='write')
        write(u,'(a)') '&oce_dyn'
        write(u,'(a)') "  mix_scheme = 'PP'"
        write(u,'(a)') '  Fer_GM = .false.'
        write(u,'(a)') '  Redi   = .false.'
        write(u,'(a)') '  A_ver  = 1.e-4'
        write(u,'(a)') '/'
        close(u)
    end subroutine

    subroutine write_bad_config()
        ! Unknown variable in &timestep -> namelist runtime parse error.
        integer :: u
        open(newunit=u, file='test_params.bad', status='replace', action='write')
        write(u,'(a)') '&timestep'
        write(u,'(a)') '  step_per_day = 32'
        write(u,'(a)') '  this_variable_does_not_exist = 5'
        write(u,'(a)') '/'
        close(u)
    end subroutine

end program test_params
