program test_io_means
    ! M9 Task 2.6 GATE for the namelist.io parser (mod_io_means::means_read_namelist). The smoke driver
    ! (fesom_outputsmoke) exercises the writer KNOBS via FESOM3_* env, but NOT the namelist read itself —
    ! this test does. It writes a known namelist.io, parses it, and asserts the &nml_general global knobs
    ! and the &nml_list rows round-trip. Crucially it pins the FESOM2 derived-type-array namelist trick
    ! (one t_io_entry filled per 5 flat values 'id',freq,'unit',prec,'op') on BOTH Intel and GNU — the
    ! one portability risk of the parser. Pure parse: no mesh / MPI collectives needed.
    use mpi
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex
    use mod_io_means,     only: t_io_config, t_io_entry, means_read_namelist, MEANS_MAXF
    implicit none

    type(t_partit)    :: partit
    type(t_io_config) :: cfg
    type(t_io_entry)  :: list(MEANS_MAXF)
    integer           :: nlist, u, nfail
    logical           :: ok
    character(len=*), parameter :: path = 'test_namelist_io.tmp'

    call par_init(partit)
    nfail = 0

    if (partit%mype == 0) then
        ! ---- write a known namelist.io (non-default knobs + a 4-row io_list) -----------------------
        open(newunit=u, file=path, status='replace', action='write', form='formatted')
        write(u,'(a)') '&nml_general'
        write(u,'(a)') '  n_writers      = 3'
        write(u,'(a)') '  chunk_time     = 4'
        write(u,'(a)') '  chunk_vert     = 12'
        write(u,'(a)') '  chunk_horiz    = 250000'
        write(u,'(a)') "  compressor     = 'lz4'"
        write(u,'(a)') "  filesplit_freq = 'm'"
        write(u,'(a)') "  vec_frame      = 'native'"
        write(u,'(a)') '/'
        write(u,'(a)') '&nml_list'
        write(u,'(a)') "  io_list = 'ssh   ', 1, 'm', 4, 'snap',"
        write(u,'(a)') "            'temp  ', 2, 'd', 8, 'mean',"
        write(u,'(a)') "            'unod  ', 1, 'y', 4, 'mean',"
        write(u,'(a)') "            'vnod  ', 1, 'y', 4, 'mean',"
        write(u,'(a)') '/'
        close(u)

        ! ---- parse + assert ------------------------------------------------------------------------
        call means_read_namelist(path, cfg, list, nlist, ok)
        call expect_l(ok, .true., 'parse ok', nfail)
        ! &nml_general
        call expect_i(cfg%n_writers,   3,      'n_writers',   nfail)
        call expect_i(cfg%chunk_time,  4,      'chunk_time',  nfail)
        call expect_i(cfg%chunk_vert,  12,     'chunk_vert',  nfail)
        call expect_i(cfg%chunk_horiz, 250000, 'chunk_horiz', nfail)
        call expect_s(trim(cfg%compressor),     'lz4',    'compressor',     nfail)
        call expect_s(trim(cfg%filesplit_freq), 'm',      'filesplit_freq', nfail)
        call expect_s(trim(cfg%vec_frame),      'native', 'vec_frame',      nfail)
        ! &nml_list (the derived-type-array parse — 5 values per entry)
        call expect_i(nlist, 4, 'nlist', nfail)
        call expect_row(list(1), 'ssh',  1, 'm', 4, 'snap', nfail)
        call expect_row(list(2), 'temp', 2, 'd', 8, 'mean', nfail)
        call expect_row(list(3), 'unod', 1, 'y', 4, 'mean', nfail)
        call expect_row(list(4), 'vnod', 1, 'y', 4, 'mean', nfail)

        ! ---- a missing file must report ok=.false. (env-only fallback) -----------------------------
        call means_read_namelist('test_namelist_io_absent.tmp', cfg, list, nlist, ok)
        call expect_l(ok, .false., 'absent file => ok=.false.', nfail)

        open(newunit=u, file=path, status='old'); close(u, status='delete')

        if (nfail == 0) then
            write(*,'(a)') 'test_io_means: PASS (namelist.io &nml_general + &nml_list parsed correctly)'
        else
            write(*,'(a,i0,a)') 'test_io_means: FAIL (', nfail, ' check(s))'
        end if
    end if

    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    if (nfail /= 0) error stop 1

contains

    subroutine expect_i(got, want, name, nf)
        integer,          intent(in)    :: got, want
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: nf
        if (got /= want) then
            write(*,'(a,i0,a,i0)') '  BAD '//name//': got ', got, ' want ', want
            nf = nf + 1
        end if
    end subroutine expect_i

    subroutine expect_s(got, want, name, nf)
        character(len=*), intent(in)    :: got, want, name
        integer,          intent(inout) :: nf
        if (trim(got) /= trim(want)) then
            write(*,'(a)') '  BAD '//name//": got '"//trim(got)//"' want '"//trim(want)//"'"
            nf = nf + 1
        end if
    end subroutine expect_s

    subroutine expect_l(got, want, name, nf)
        logical,          intent(in)    :: got, want
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: nf
        if (got .neqv. want) then
            write(*,'(a,l1,a,l1)') '  BAD '//name//': got ', got, ' want ', want
            nf = nf + 1
        end if
    end subroutine expect_l

    subroutine expect_row(e, id, freq, unit, prec, op, nf)
        type(t_io_entry), intent(in)    :: e
        character(len=*), intent(in)    :: id, unit, op
        integer,          intent(in)    :: freq, prec
        integer,          intent(inout) :: nf
        call expect_s(trim(e%id),   id,   'row '//id//' id',   nf)
        call expect_i(e%freq,       freq, 'row '//id//' freq', nf)
        call expect_s(trim(e%unit), unit, 'row '//id//' unit', nf)
        call expect_i(e%precision,  prec, 'row '//id//' prec', nf)
        call expect_s(trim(e%op),   op,   'row '//id//' op',   nf)
    end subroutine expect_row

end program test_io_means
