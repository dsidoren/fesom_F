module mod_advhor_dump
    ! Write FESOM3's horizontal-advection fields in the SAME self-describing binary
    ! format as the FESOM2 oracle's src/fesom_advhor_dump.F90, so tools/advhor_diff.py
    ! compares them field-by-field for the M1.1 operator byte-gate. Single-rank only
    ! (local == global), matching the FESOM2 1-rank advection-gate run.
    !
    ! Format (stream, little-endian) — same record layout as mod_geom_dump, magic
    ! "FADVHDMP":
    !   char8  "FADVHDMP" | int32 nod2D,elem2D,edge2D,nl
    !   records: char24 name | int32 dtype(0=r64,1=i32) | int32 d1,d2 | data
    ! 3-D real arrays are written with d1 = size(dim1)*size(dim2), d2 = size(dim3)
    ! (column-major contiguous), the same memory order on both sides.
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use mod_precision, only: MP
    implicit none
    private
    public :: advhor_dump_open, advhor_dump_close, wr_r2, wr_r3, wr_i1, wr_i2

contains

    subroutine advhor_dump_open(u, path, nod2D, elem2D, edge2D, nl)
        integer,          intent(out) :: u
        character(len=*), intent(in)  :: path
        integer,          intent(in)  :: nod2D, elem2D, edge2D, nl
        integer :: ios
        open(newunit=u, file=trim(path), status='replace', form='unformatted', &
             access='stream', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,'(a)') 'advhor_dump_open: cannot open '//trim(path); error stop 1
        end if
        write(u) 'FADVHDMP'
        write(u) int(nod2D,int32), int(elem2D,int32), int(edge2D,int32), int(nl,int32)
    end subroutine advhor_dump_open

    subroutine advhor_dump_close(u)
        integer, intent(in) :: u
        close(u)
    end subroutine advhor_dump_close

    subroutine wr_r2(u, name, a)
        integer, intent(in) :: u
        character(len=*), intent(in) :: name
        real(kind=MP), intent(in) :: a(:,:)
        character(len=24) :: nm
        nm = name
        write(u) nm, int(0,int32), int(size(a,1),int32), int(size(a,2),int32)
        write(u) real(a, real64)
    end subroutine wr_r2

    subroutine wr_r3(u, name, a)
        integer, intent(in) :: u
        character(len=*), intent(in) :: name
        real(kind=MP), intent(in) :: a(:,:,:)
        character(len=24) :: nm
        nm = name
        write(u) nm, int(0,int32), int(size(a,1)*size(a,2),int32), int(size(a,3),int32)
        write(u) real(a, real64)
    end subroutine wr_r3

    subroutine wr_i1(u, name, a)
        integer, intent(in) :: u
        character(len=*), intent(in) :: name
        integer, intent(in) :: a(:)
        character(len=24) :: nm
        nm = name
        write(u) nm, int(1,int32), int(size(a),int32), int(1,int32)
        write(u) int(a, int32)
    end subroutine wr_i1

    subroutine wr_i2(u, name, a)
        integer, intent(in) :: u
        character(len=*), intent(in) :: name
        integer, intent(in) :: a(:,:)
        character(len=24) :: nm
        nm = name
        write(u) nm, int(1,int32), int(size(a,1),int32), int(size(a,2),int32)
        write(u) int(a, int32)
    end subroutine wr_i2

end module mod_advhor_dump
