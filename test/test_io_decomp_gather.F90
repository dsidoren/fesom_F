program test_io_decomp_gather
    ! Task 1.1 GATE for mod_io_decomp: gather ∘ redistribute == identity on owned entities.
    !
    ! Build a synthetic t_io_decomp via decomp_init (the module is unit-testable without a mesh —
    ! pass a round-robin myList of canonical ids), fill a DETERMINISTIC field_owned (no RNG, so it is
    ! reproducible), push it to the canonical writer layout (decomp_redistribute), pull it back
    ! (decomp_gather), and assert field_owned2 == field_owned BITWISE (max|Δ|=0) for every owned i.
    !
    ! Round-robin ownership (rank r owns global ids g with mod(g-1,npes)==r) makes each writer receive
    ! from BOTH ranks at np=2 (real cross-rank exchange); at np=1 the Alltoallv is a self-copy. The
    ! store holds only canonical owned values, so the partial last chunk's pad slots are never read by
    ! gather (recv_target indexes real entities only). Covers a 2D and a 3D field over several
    ! (N, C, n_writers) combos incl. partial last chunks and writer SUBSETS (n_writers<npes).
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

    call run_case(100, 7,   npes)            ! partial last chunk (15 chunks), all writers
    call run_case(100, 10,  2)               ! 2 writers: at np=2 both own chunks => real cross-rank exchange
    call run_case(64,  8,   npes)            ! evenly divides
    call run_case(37,  5,   max(1, npes-1))  ! odd N, fewer writers
    call run_case(100, 1,   npes)            ! C=1 (100 chunks)
    call run_case(100, 100, 1)               ! single chunk, single writer

    call MPI_Allreduce(MPI_IN_PLACE, nfail, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
    if (mype == 0) then
        if (nfail == 0) then
            write(*,'(a)') 'test_io_decomp_gather: PASS (gather∘redistribute == identity, max|Δ|=0, 2D+3D)'
        else
            write(*,'(a,i0,a)') 'test_io_decomp_gather: FAIL (', nfail, ' mismatches)'
        end if
    end if
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    if (nfail /= 0) error stop 1

contains

    subroutine run_case(N, C, nw)
        integer, intent(in) :: N, C, nw
        type(t_io_decomp) :: D
        integer, allocatable :: myList(:)
        real(WP), allocatable :: f2(:), buf2(:), g2(:)
        real(WP), allocatable :: f3(:,:), buf3(:,:), g3(:,:)
        integer :: myDim, i, g, k, L
        integer, parameter :: NLEV = 4

        ! round-robin ownership: rank r owns global ids g with mod(g-1,npes)==r
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

        ! --- 2D: deterministic owned field -> redistribute -> gather -> assert bitwise identity ---
        allocate(f2(max(1,myDim)), buf2(max(1,D%w_nbuf)), g2(max(1,myDim)))
        do i = 1, myDim
            f2(i) = real(myList(i), WP) * 1.5_WP - 0.25_WP
        end do
        g2 = -huge(1.0_WP)                          ! sentinel: any unfilled owned slot fails the compare
        call decomp_redistribute(D, f2, buf2, -1.0_WP)
        call decomp_gather(D, buf2, g2)
        do i = 1, myDim
            if (g2(i) /= f2(i)) call bump(N, C, nw, 'r2', myList(i))
        end do

        ! --- 3D: per-level offset -> redistribute -> gather -> assert bitwise identity ---
        allocate(f3(NLEV, max(1,myDim)), buf3(NLEV, max(1,D%w_nbuf)), g3(NLEV, max(1,myDim)))
        do i = 1, myDim
            do L = 1, NLEV
                f3(L,i) = real(myList(i), WP) * 1.5_WP - 0.25_WP + real(L, WP) * 1000.0_WP
            end do
        end do
        g3 = -huge(1.0_WP)
        call decomp_redistribute(D, f3, buf3, -1.0_WP)
        call decomp_gather(D, buf3, g3)
        do i = 1, myDim
            do L = 1, NLEV
                if (g3(L,i) /= f3(L,i)) call bump(N, C, nw, 'r3', myList(i))
            end do
        end do

        deallocate(myList, f2, buf2, g2, f3, buf3, g3)
    end subroutine run_case

    subroutine bump(N, C, nw, tag, g)
        integer,          intent(in) :: N, C, nw, g
        character(len=*), intent(in) :: tag
        nfail = nfail + 1
        if (nfail <= 6) write(*,'(a,i0,a,4(a,i0))') '  [rank ', mype, '] FAIL ', &
            trim(tag)//' N=', N, ' C=', C, ' nw=', nw, ' gid=', g
    end subroutine bump

end program test_io_decomp_gather
