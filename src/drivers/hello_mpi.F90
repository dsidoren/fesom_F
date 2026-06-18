program hello_mpi
    ! M0.1 build-system smoke test: links the (near-empty) fesom3 library,
    ! initialises MPI, reports rank/size, and finalises cleanly.
    use mpi
    use mod_fesom_version, only: fesom3_version_string
    implicit none

    integer :: ierr, rank, nprocs

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank,    ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

    if (rank == 0) then
        write(*,'(a)')        trim(fesom3_version_string())//" — hello-MPI"
        write(*,'(a,i0,a)')   "Running on ", nprocs, " rank(s)."
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    write(*,'(a,i0,a,i0)') "  rank ", rank, " of ", nprocs

    call MPI_Finalize(ierr)
end program hello_mpi
