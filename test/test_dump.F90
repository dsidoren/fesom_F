program test_dump
    ! M0.6 gate (Fortran half): write node + element dumps via mod_dump, read the
    ! raw binary back, verify record structure + values round-trip exactly. The
    ! dump_diff.py self-check is a separate ctest (tools/dump_diff.py --selftest).
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use mpi
    use mod_precision, only: WP
    use mod_dump
    implicit none

    integer, parameter :: NN = 3140, NE = 5839, NL = 5
    real(kind=WP), allocatable :: density(:,:), eta(:), Av(:,:), pgf_x(:)
    integer,       allocatable :: myList_nod(:), myList_elem(:), nlev_nod(:), nlev_elem(:)
    integer :: ierr, i, j, nfail

    nfail = 0
    call MPI_Init(ierr)

    allocate(myList_nod(NN), myList_elem(NE), nlev_nod(NN), nlev_elem(NE))
    myList_nod  = [(i, i=1, NN)]    ! identity: probe gids resolve to local==gid
    myList_elem = [(i, i=1, NE)]
    nlev_nod = NL; nlev_elem = NL

    allocate(density(NL, NN), eta(NN), Av(NL, NE), pgf_x(NE))
    do j = 1, NN
        do i = 1, NL
            density(i, j) = real(j*100 + i, WP)
        end do
        eta(j) = real(j*100, WP)
    end do
    do j = 1, NE
        do i = 1, NL
            Av(i, j) = real(j*10 + i, WP)
        end do
        pgf_x(j) = real(j*10, WP)
    end do

    call dump_init(0, NN, myList_nod, NE, myList_elem, &
                   node_prefix='test_dump_nod', elem_prefix='test_dump_elem')
    call check(dump_is_active(), 'dump active with explicit prefixes')

    call dump_node   (DUMP_SUBSTEP_PRESSURE_BV, 1, 'density', density, nlev_nod)
    call dump_node_2d(DUMP_SUBSTEP_ETA_N,       1, 'eta',     eta)
    call dump_elem   (DUMP_SUBSTEP_PGF,         1, 'Av',      Av, nlev_elem)
    call dump_elem_2d(DUMP_SUBSTEP_PGF,         1, 'pgf_x',   pgf_x)
    call dump_finalize()

    call verify_node()
    call verify_elem()

    if (nfail == 0) then
        write(*,'(a)') 'test_dump: ALL PASS'
    else
        write(*,'(a,i0,a)') 'test_dump: ', nfail, ' FAILURE(S)'
    end if
    call MPI_Finalize(ierr)
    if (nfail /= 0) error stop 1

contains

    subroutine check(cond, name)
        logical, intent(in) :: cond
        character(len=*), intent(in) :: name
        if (.not. cond) then
            nfail = nfail + 1
            write(*,'(a)') '  FAIL: '//name
        end if
    end subroutine

    subroutine verify_node()
        integer :: u, ios, nrec, density_recs, eta_recs
        integer(int32) :: step, substep, gid, nlev
        character(len=24) :: name24
        real(real64) :: vals(64)
        open(newunit=u, file='test_dump_nod.00000', form='unformatted', &
             access='stream', status='old', action='read', iostat=ios)
        call check(ios == 0, 'node dump file exists')
        if (ios /= 0) return
        nrec = 0; density_recs = 0; eta_recs = 0
        do
            read(u, iostat=ios) step, substep, gid, nlev
            if (ios /= 0) exit
            read(u, iostat=ios) name24
            read(u, iostat=ios) vals(1:nlev)
            nrec = nrec + 1
            call check(step == 1, 'node rec step==1')
            select case (trim(name24))
            case ('density')
                density_recs = density_recs + 1
                call check(substep == DUMP_SUBSTEP_PRESSURE_BV, 'density substep')
                call check(nlev == NL, 'density nlev==NL')
                call check(vals(1) == real(gid*100 + 1, real64), 'density value(1)')
                call check(vals(nlev) == real(gid*100 + NL, real64), 'density value(NL)')
            case ('eta')
                eta_recs = eta_recs + 1
                call check(substep == DUMP_SUBSTEP_ETA_N, 'eta substep')
                call check(nlev == 1, 'eta nlev==1')
                call check(vals(1) == real(gid*100, real64), 'eta value')
            case default
                call check(.false., 'unexpected node field '//trim(name24))
            end select
        end do
        close(u)
        call check(density_recs == DUMP_NPROBES_NOD, 'density: 1 record per probe')
        call check(eta_recs == DUMP_NPROBES_NOD, 'eta: 1 record per probe')
    end subroutine verify_node

    subroutine verify_elem()
        integer :: u, ios, av_recs, pgf_recs
        integer(int32) :: step, substep, gid, nlev
        character(len=24) :: name24
        real(real64) :: vals(64)
        open(newunit=u, file='test_dump_elem.00000', form='unformatted', &
             access='stream', status='old', action='read', iostat=ios)
        call check(ios == 0, 'elem dump file exists')
        if (ios /= 0) return
        av_recs = 0; pgf_recs = 0
        do
            read(u, iostat=ios) step, substep, gid, nlev
            if (ios /= 0) exit
            read(u, iostat=ios) name24
            read(u, iostat=ios) vals(1:nlev)
            select case (trim(name24))
            case ('Av')
                av_recs = av_recs + 1
                call check(nlev == NL, 'Av nlev==NL')
                call check(vals(1) == real(gid*10 + 1, real64), 'Av value(1)')
            case ('pgf_x')
                pgf_recs = pgf_recs + 1
                call check(nlev == 1, 'pgf_x nlev==1')
                call check(vals(1) == real(gid*10, real64), 'pgf_x value')
            case default
                call check(.false., 'unexpected elem field '//trim(name24))
            end select
        end do
        close(u)
        call check(av_recs == DUMP_NPROBES_ELEM, 'Av: 1 record per probe')
        call check(pgf_recs == DUMP_NPROBES_ELEM, 'pgf_x: 1 record per probe')
    end subroutine verify_elem

end program test_dump
