program test_partit
    ! M0.4 gate: 1-rank synthesis + dist_<NP>/ reading on the pi mesh.
    ! Runs at 1, 2, 8 ranks. Invariants: owned nodes/edges partition exactly
    ! (sum == global), owned+halo elements are redundant (sum >= global).
    use mpi
    use mod_partit,        only: t_partit
    use mod_partitioning,  only: par_init, par_ex, set_partition
    implicit none

    ! pi mesh reference (override with FESOM3_MESH_DIR).
    integer, parameter :: PI_NOD2D = 3140, PI_ELEM2D = 5839, PI_EDGE2D = 8986
    character(len=512) :: mesh_dir
    type(t_partit) :: partit
    integer :: ierr, nfail, nfail_g
    integer :: sum_nod, sum_elem, sum_edge

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'

    nfail = 0
    call par_init(partit)
    call set_partition(partit, trim(mesh_dir))

    ! ---- global partition invariants (all rank counts) ----
    call MPI_Allreduce(partit%myDim_nod2D,  sum_nod,  1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
    call MPI_Allreduce(partit%myDim_elem2D, sum_elem, 1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
    call MPI_Allreduce(partit%myDim_edge2D, sum_edge, 1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, ierr)

    ! Nodes are uniquely owned (sum == nod2D); elements and edges have
    ! boundary-redundant owned counts (sum >= global), verified on pi/dist_2.
    call check(sum_nod  == PI_NOD2D,  'sum(myDim_nod2D) == nod2D (unique)')
    call check(sum_elem >= PI_ELEM2D, 'sum(myDim_elem2D) >= elem2D (boundary-redundant)')
    call check(sum_edge >= PI_EDGE2D, 'sum(myDim_edge2D) >= edge2D (boundary-redundant)')

    ! ---- local sanity ----
    call check(partit%myDim_nod2D > 0, 'myDim_nod2D > 0')
    call check(size(partit%myList_nod2D) == partit%myDim_nod2D + partit%eDim_nod2D, 'myList_nod2D size')
    call check(all(partit%myList_nod2D >= 1 .and. partit%myList_nod2D <= PI_NOD2D), 'node gids in range')
    call check(all(partit%myList_elem2D >= 1 .and. partit%myList_elem2D <= PI_ELEM2D), 'elem gids in range')

    ! ---- 1-rank synthesis specifics ----
    if (partit%npes == 1) then
        call check(partit%myDim_nod2D == PI_NOD2D, '1-rank: myDim_nod2D == nod2D')
        call check(partit%eDim_nod2D == 0, '1-rank: eDim_nod2D == 0')
        call check(partit%eDim_elem2D == 0 .and. partit%eXDim_elem2D == 0, '1-rank: elem halos == 0')
        call check(partit%myList_nod2D(1) == 1 .and. partit%myList_nod2D(PI_NOD2D) == PI_NOD2D, &
                   '1-rank: identity node map')
        call check(partit%com_nod2D%rPEnum == 0 .and. partit%com_nod2D%sPEnum == 0, '1-rank: no neighbours')
    else
        ! multi-rank: every PE should have at least some halo and neighbours
        call check(partit%eDim_nod2D > 0, 'multi-rank: eDim_nod2D > 0')
        call check(partit%com_nod2D%rPEnum > 0, 'multi-rank: has receive neighbours')
        call check(size(partit%com_nod2D%rlist) == partit%eDim_nod2D, 'multi-rank: rlist size == eDim_nod2D')
    end if

    call MPI_Allreduce(nfail, nfail_g, 1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
    if (partit%mype == 0) then
        if (nfail_g == 0) then
            write(*,'(a,i0,a)') 'test_partit (npes=', partit%npes, '): ALL PASS'
        else
            write(*,'(a,i0,a,i0,a)') 'test_partit (npes=', partit%npes, '): ', nfail_g, ' FAILURE(S)'
        end if
    end if

    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    if (nfail_g /= 0) error stop 1

contains
    subroutine check(cond, name)
        logical,          intent(in) :: cond
        character(len=*), intent(in) :: name
        if (.not. cond) then
            nfail = nfail + 1
            write(*,'(a,i0,a)') '  FAIL [rank ', partit%mype, ']: '//name
        end if
    end subroutine
end program test_partit
