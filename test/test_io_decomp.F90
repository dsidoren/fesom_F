program test_io_decomp
    ! M9 Task 1.1 GATE for mod_io_decomp: redistribute a LABELLED field (value = canonical id) from a
    ! synthetic compute partition to the canonical chunked writer layout, and assert every writer slot
    ! holds its canonical id (pad slots = fill). Synthetic partition = ROUND-ROBIN ownership (rank r
    ! owns global ids g with mod(g-1,npes)==r) so the Alltoallv genuinely reorders — no mesh needed.
    !
    ! Covers the plan's "1-rank buf[g]==g" (np=1: round-robin collapses to identity, Alltoallv = self
    ! copy) AND the multi-rank redistribution (np=2/8) in isolation, for real / int / 3D fields, over
    ! several (N, C, n_writers) combos incl. partial last chunks and writer SUBSETS (n_writers<npes).
    use mpi
    use mod_precision,    only: WP
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex
    use mod_io_decomp
    implicit none

    type(t_partit) :: partit
    integer :: npes, mype, comm, ierr, nfail

    call par_init(partit)
    comm = partit%MPI_COMM_FESOM; mype = partit%mype; npes = partit%npes
    nfail = 0

    call run_case(100, 7,   npes)            ! partial chunks (15), all writers
    call run_case(100, 10,  2)               ! writer SUBSET (2 of npes)
    call run_case(100, 1,   npes)            ! C=1 (100 chunks)
    call run_case(100, 100, 1)               ! single chunk, single writer
    call run_case(37,  5,   max(1, npes-1))  ! odd N, fewer writers
    call run_case(64,  8,   npes)            ! evenly divides

    call MPI_Allreduce(MPI_IN_PLACE, nfail, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
    if (mype == 0) then
        if (nfail == 0) then
            write(*,'(a)') 'test_io_decomp: PASS (canonical redistribution real+int+3D, max|Δ|=0)'
        else
            write(*,'(a,i0,a)') 'test_io_decomp: FAIL (', nfail, ' mismatches)'
        end if
    end if
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    if (nfail /= 0) error stop 1

contains

    subroutine run_case(N, C, nw)
        integer, intent(in) :: N, C, nw
        type(t_io_decomp) :: D
        integer, allocatable :: myList(:), fi(:), bufi(:)
        real(WP), allocatable :: f2(:), buf2(:), f3(:,:), buf3(:,:)
        integer :: myDim, i, g, k, L
        integer, parameter :: NLEV = 4

        ! round-robin ownership
        myDim = 0
        do g = 1, N
            if (mod(g-1, npes) == mype) myDim = myDim + 1
        end do
        allocate(myList(max(1,myDim)))
        k = 0
        do g = 1, N
            if (mod(g-1, npes) == mype) then
                k = k + 1; myList(k) = g
            end if
        end do

        call decomp_init(D, C, nw, N, myList, myDim, comm, mype, npes)

        ! 2D real: value = gid
        allocate(f2(max(1,myDim)), buf2(max(1,D%w_nbuf)))
        do i = 1, myDim; f2(i) = real(myList(i), WP); end do
        call decomp_redistribute(D, f2, buf2, -1.0_WP)
        do k = 1, D%w_nbuf
            g = D%w_base_gid + k
            if (g <= N) then
                if (buf2(k) /= real(g, WP))   call bump(N, C, nw, 'r2', g)
            else
                if (buf2(k) /= -1.0_WP)       call bump(N, C, nw, 'r2pad', g)
            end if
        end do

        ! 2D int: value = gid*10
        allocate(fi(max(1,myDim)), bufi(max(1,D%w_nbuf)))
        do i = 1, myDim; fi(i) = myList(i)*10; end do
        call decomp_redistribute(D, fi, bufi, -999)
        do k = 1, D%w_nbuf
            g = D%w_base_gid + k
            if (g <= N) then
                if (bufi(k) /= g*10)          call bump(N, C, nw, 'i2', g)
            else
                if (bufi(k) /= -999)          call bump(N, C, nw, 'i2pad', g)
            end if
        end do

        ! 3D real: value(L) = gid*1000 + L
        allocate(f3(NLEV, max(1,myDim)), buf3(NLEV, max(1,D%w_nbuf)))
        do i = 1, myDim
            do L = 1, NLEV; f3(L,i) = real(myList(i)*1000 + L, WP); end do
        end do
        call decomp_redistribute(D, f3, buf3, -1.0_WP)
        do k = 1, D%w_nbuf
            g = D%w_base_gid + k
            do L = 1, NLEV
                if (g <= N) then
                    if (buf3(L,k) /= real(g*1000 + L, WP)) call bump(N, C, nw, 'r3', g)
                else
                    if (buf3(L,k) /= -1.0_WP)              call bump(N, C, nw, 'r3pad', g)
                end if
            end do
        end do

        deallocate(myList, f2, buf2, fi, bufi, f3, buf3)
    end subroutine run_case

    subroutine bump(N, C, nw, tag, g)
        integer,          intent(in) :: N, C, nw, g
        character(len=*), intent(in) :: tag
        nfail = nfail + 1
        if (nfail <= 6) write(*,'(a,i0,a,4(a,i0))') '  [rank ', mype, '] FAIL ', &
            trim(tag)//' N=', N, ' C=', C, ' nw=', nw, ' gid=', g
    end subroutine bump

end program test_io_decomp
