program fesom_zarrsmoke
    ! Stage-0 smoke test for mod_io_zarr (M9 Task 0.2 + 0.3). Writes TWO Zarr v2 stores with three
    ! toy arrays each, identical values, differing only in codec + consolidation:
    !   zarrsmoke.zarr      codec none, consolidated
    !   zarrsmoke_lz4.zarr  codec lz4,  consolidated
    ! Arrays (chunks deliberately do NOT divide the shape => partial-chunk padding; distinct per-cell
    ! values => a C-order transpose bug is caught by tools/zarr_diff.py):
    !   arr1d_f8   shape [7]   chunks [3]   value(i)     = i*1.5
    !   arr2d_f4   shape [5,3] chunks [2,2] value(i0,i1) = i0*100 + i1
    !   arr2d_i4   shape [5,3] chunks [2,2] value(i0,i1) = i0*100 + i1
    !
    !   FESOM3_ZARR_OUT   output directory for the stores (default: ./)
    use, intrinsic :: iso_fortran_env, only: real64
    use mod_io_zarr
    implicit none

    character(len=4096) :: outdir

    call get_environment_variable('FESOM3_ZARR_OUT', outdir)
    if (len_trim(outdir) == 0) outdir = '.'

    call write_smoke(trim(outdir)//'/zarrsmoke.zarr',     'none')
    call write_smoke(trim(outdir)//'/zarrsmoke_lz4.zarr', 'lz4')

    print '(A)', 'ZARRSMOKE OK'

contains

    subroutine write_smoke(store_path, codec)
        character(len=*), intent(in) :: store_path, codec
        type(t_zarr_store) :: store
        type(t_zarr_array) :: a1, a2, a3
        type(t_zarr_attrs) :: at
        real(real64) :: d1(7)
        real(real64) :: d2(5,3)
        integer      :: d3(5,3)
        integer      :: i0, i1

        ! store + group attrs
        call zarr_attrs_init(at)
        call zattr_str(at, 'Conventions', 'FESOM3-zarrsmoke')
        call zattr_int(at, 'smoke_version', 1)
        call zarr_create_store(store, store_path, at)

        ! 1-D f8
        do i0 = 1, 7
            d1(i0) = real(i0, real64) * 1.5_real64
        end do
        call zarr_array_init(a1, 'arr1d_f8', dims=[7], chunks=[3], dtype='<f8', codec=codec)
        call zarr_attrs_init(at)
        call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=2) :: 'x'])
        call zattr_str(at, 'units', 'meter')
        call zarr_define_array(store, a1, at)
        call zarr_write_whole(store, a1, d1)

        ! 2-D f4
        do i1 = 1, 3
            do i0 = 1, 5
                d2(i0,i1) = real(i0*100 + i1, real64)
            end do
        end do
        call zarr_array_init(a2, 'arr2d_f4', dims=[5,3], chunks=[2,2], dtype='<f4', codec=codec)
        call zarr_attrs_init(at)
        call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=2) :: 'd0','d1'])
        call zarr_define_array(store, a2, at)
        call zarr_write_whole(store, a2, d2)

        ! 2-D i4
        do i1 = 1, 3
            do i0 = 1, 5
                d3(i0,i1) = i0*100 + i1
            end do
        end do
        call zarr_array_init(a3, 'arr2d_i4', dims=[5,3], chunks=[2,2], dtype='<i4', codec=codec)
        call zarr_attrs_init(at)
        call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=2) :: 'd0','d1'])
        call zarr_define_array(store, a3, at)
        call zarr_write_whole(store, a3, d3)

        call zarr_consolidate(store)
    end subroutine write_smoke

end program fesom_zarrsmoke
