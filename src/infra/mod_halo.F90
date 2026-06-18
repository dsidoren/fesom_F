module mod_halo
    ! Generic MPI halo exchange (decision D7). The generic exchange_nod/exchange_elem
    ! interfaces + com_struct ARE the language-agnostic portability contract; MPI is
    ! the single implementation (manual pack -> Isend/Irecv -> unpack). Halo is
    ! BROADCAST-ONLY: each halo entry is OVERWRITTEN with its owner's value (never
    ! additive). rlist/slist hold LOCAL indices; the field array is indexed locally
    ! (1..myDim+eDim). Structure follows the dwarf gen_halo_exchange.F90, written
    ! clean (GPU/g2g/OpenACC bloat omitted; luse_g2g kept as an inert foundation arg).
    !
    ! Deferred (optimizations, not needed for byte-identity): begin/end overlap split,
    ! multi-field batched exchange, com_elem2D_full default path.
    use mpi
    use mod_precision, only: WP, MPI_WP
    use mod_partit,    only: t_partit, com_struct
    implicit none
    private
    public :: exchange_nod, exchange_elem, stale_halo_max_nod, stale_halo_max_elem

    interface exchange_nod
        module procedure exchange_nod_2D_r, exchange_nod_3D_r, exchange_nod_2D_i
    end interface

    interface exchange_elem
        module procedure exchange_elem_2D_r, exchange_elem_3D_r
    end interface

    integer, parameter :: HALO_TAG = 1

contains

    ! ===================== node wrappers =====================
    subroutine exchange_nod_2D_r(arr, partit, luse_g2g)
        real(kind=WP),  intent(inout) :: arr(:)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g   ! GPU-direct hook (inert in v1)
        call core_2D_r(arr, partit%com_nod2D, partit%MPI_COMM_FESOM)
    end subroutine

    subroutine exchange_nod_3D_r(arr, partit, luse_g2g)
        real(kind=WP),  intent(inout) :: arr(:,:)   ! (nl, nod_size)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call core_3D_r(arr, partit%com_nod2D, partit%MPI_COMM_FESOM)
    end subroutine

    subroutine exchange_nod_2D_i(arr, partit, luse_g2g)
        integer,        intent(inout) :: arr(:)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call core_2D_i(arr, partit%com_nod2D, partit%MPI_COMM_FESOM)
    end subroutine

    ! ===================== element wrappers ==================
    subroutine exchange_elem_2D_r(arr, partit, luse_g2g)
        real(kind=WP),  intent(inout) :: arr(:)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call core_2D_r(arr, partit%com_elem2D, partit%MPI_COMM_FESOM)
    end subroutine

    subroutine exchange_elem_3D_r(arr, partit, luse_g2g)
        real(kind=WP),  intent(inout) :: arr(:,:)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call core_3D_r(arr, partit%com_elem2D, partit%MPI_COMM_FESOM)
    end subroutine

    ! ===================== cores =============================
    subroutine core_2D_r(arr, com, comm)
        real(kind=WP),    intent(inout) :: arr(:)
        type(com_struct), intent(in)    :: com
        integer,          intent(in)    :: comm
        real(kind=WP), allocatable :: sbuf(:), rbuf(:)
        integer, allocatable :: req(:)
        integer :: i, off, n, nreq, ierr
        if (com%rPEnum == 0 .and. com%sPEnum == 0) return
        allocate(req(com%rPEnum + com%sPEnum))
        allocate(rbuf(max(com%rptr(com%rPEnum+1)-1, 0)))
        allocate(sbuf(max(com%sptr(com%sPEnum+1)-1, 0)))
        nreq = 0
        do i = 1, com%rPEnum
            off = com%rptr(i); n = com%rptr(i+1) - com%rptr(i)
            nreq = nreq + 1
            call MPI_Irecv(rbuf(off), n, MPI_WP, com%rPE(i), HALO_TAG, comm, req(nreq), ierr)
        end do
        do i = 1, com%sPEnum
            off = com%sptr(i); n = com%sptr(i+1) - com%sptr(i)
            sbuf(off:off+n-1) = arr(com%slist(off:off+n-1))
            nreq = nreq + 1
            call MPI_Isend(sbuf(off), n, MPI_WP, com%sPE(i), HALO_TAG, comm, req(nreq), ierr)
        end do
        call MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, ierr)
        do i = 1, com%rPEnum
            off = com%rptr(i); n = com%rptr(i+1) - com%rptr(i)
            arr(com%rlist(off:off+n-1)) = rbuf(off:off+n-1)   ! broadcast-only
        end do
    end subroutine core_2D_r

    subroutine core_2D_i(arr, com, comm)
        integer,          intent(inout) :: arr(:)
        type(com_struct), intent(in)    :: com
        integer,          intent(in)    :: comm
        integer, allocatable :: sbuf(:), rbuf(:), req(:)
        integer :: i, off, n, nreq, ierr
        if (com%rPEnum == 0 .and. com%sPEnum == 0) return
        allocate(req(com%rPEnum + com%sPEnum))
        allocate(rbuf(max(com%rptr(com%rPEnum+1)-1, 0)))
        allocate(sbuf(max(com%sptr(com%sPEnum+1)-1, 0)))
        nreq = 0
        do i = 1, com%rPEnum
            off = com%rptr(i); n = com%rptr(i+1) - com%rptr(i)
            nreq = nreq + 1
            call MPI_Irecv(rbuf(off), n, MPI_INTEGER, com%rPE(i), HALO_TAG, comm, req(nreq), ierr)
        end do
        do i = 1, com%sPEnum
            off = com%sptr(i); n = com%sptr(i+1) - com%sptr(i)
            sbuf(off:off+n-1) = arr(com%slist(off:off+n-1))
            nreq = nreq + 1
            call MPI_Isend(sbuf(off), n, MPI_INTEGER, com%sPE(i), HALO_TAG, comm, req(nreq), ierr)
        end do
        call MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, ierr)
        do i = 1, com%rPEnum
            off = com%rptr(i); n = com%rptr(i+1) - com%rptr(i)
            arr(com%rlist(off:off+n-1)) = rbuf(off:off+n-1)
        end do
    end subroutine core_2D_i

    subroutine core_3D_r(arr, com, comm)
        real(kind=WP),    intent(inout) :: arr(:,:)   ! (nl, n)
        type(com_struct), intent(in)    :: com
        integer,          intent(in)    :: comm
        real(kind=WP), allocatable :: sbuf(:,:), rbuf(:,:)
        integer, allocatable :: req(:)
        integer :: i, off, n, nreq, ierr, nl, k
        if (com%rPEnum == 0 .and. com%sPEnum == 0) return
        nl = size(arr, 1)
        allocate(req(com%rPEnum + com%sPEnum))
        allocate(rbuf(nl, max(com%rptr(com%rPEnum+1)-1, 0)))
        allocate(sbuf(nl, max(com%sptr(com%sPEnum+1)-1, 0)))
        nreq = 0
        do i = 1, com%rPEnum
            off = com%rptr(i); n = com%rptr(i+1) - com%rptr(i)
            nreq = nreq + 1
            call MPI_Irecv(rbuf(1,off), nl*n, MPI_WP, com%rPE(i), HALO_TAG, comm, req(nreq), ierr)
        end do
        do i = 1, com%sPEnum
            off = com%sptr(i); n = com%sptr(i+1) - com%sptr(i)
            do k = 0, n-1
                sbuf(:, off+k) = arr(:, com%slist(off+k))
            end do
            nreq = nreq + 1
            call MPI_Isend(sbuf(1,off), nl*n, MPI_WP, com%sPE(i), HALO_TAG, comm, req(nreq), ierr)
        end do
        call MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, ierr)
        do i = 1, com%rPEnum
            off = com%rptr(i); n = com%rptr(i+1) - com%rptr(i)
            do k = 0, n-1
                arr(:, com%rlist(off+k)) = rbuf(:, off+k)   ! broadcast-only
            end do
        end do
    end subroutine core_3D_r

    ! ============ stale-halo probe (exchange-and-compare) ============
    real(kind=WP) function stale_halo_max_nod(arr, partit) result(maxdiff)
        real(kind=WP),  intent(in) :: arr(:)
        type(t_partit), intent(in) :: partit
        maxdiff = probe(arr, partit%com_nod2D, partit%MPI_COMM_FESOM)
    end function

    real(kind=WP) function stale_halo_max_elem(arr, partit) result(maxdiff)
        real(kind=WP),  intent(in) :: arr(:)
        type(t_partit), intent(in) :: partit
        maxdiff = probe(arr, partit%com_elem2D, partit%MPI_COMM_FESOM)
    end function

    real(kind=WP) function probe(arr, com, comm) result(maxdiff)
        ! Exchange a copy and report max|arr(halo) - exchanged(halo)|. >0 => stale halo.
        real(kind=WP),    intent(in) :: arr(:)
        type(com_struct), intent(in) :: com
        integer,          intent(in) :: comm
        real(kind=WP), allocatable :: cpy(:)
        integer :: i, off, n
        maxdiff = 0.0_WP
        cpy = arr
        call core_2D_r(cpy, com, comm)
        do i = 1, com%rPEnum
            off = com%rptr(i); n = com%rptr(i+1) - com%rptr(i)
            maxdiff = max(maxdiff, maxval(abs(arr(com%rlist(off:off+n-1)) - cpy(com%rlist(off:off+n-1)))))
        end do
    end function probe

end module mod_halo
