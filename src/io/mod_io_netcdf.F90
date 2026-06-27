module mod_io_netcdf
    ! Thin netCDF read wrapper for FESOM3 (M2.10 forcing read onward).
    !
    ! FESOM2 reads forcing through two layers — gen_surface_forcing (nf90_*, the
    ! F90 `netcdf` module) for the lon/lat/time axes, and forcing_provider_netcdf
    ! (nf_*, the F77 netcdf.inc) for the data slices via nf_get_vara_real into a
    ! real(4) buffer. The BYTES read are identical either way (IEEE deserialization,
    ! no arithmetic), so FESOM3 uses the F90 `netcdf` module uniformly. The data
    ! slices are read as real(4) (the on-disk type) and promoted to WP by the
    ! caller, exactly mirroring FESOM2 (real(4) sbcdata -> real(wp) in the bilinear).
    !
    ! Name lookups accept the same case/alias fallbacks FESOM2 uses (LON/lon/
    ! longitude, LAT/lat/latitude, TIME/time), so the same CORE2/JRA files resolve.
    use, intrinsic :: iso_fortran_env, only: real32, real64, error_unit
    use mod_precision, only: WP
    use netcdf
    implicit none
    private
    public :: nc_open_read, nc_close, nc_dimlen, nc_get_axis_dp, &
              nc_get_slice_r4, nc_get_slice_dp, nc_get_var3d_dp, &
              nc_get_att_text, nc_get_att_dp, nc_varid

contains

    subroutine nc_check(st, ctx)
        integer, intent(in) :: st
        character(len=*), intent(in) :: ctx
        if (st /= nf90_noerr) then
            write(error_unit,*) '[mod_io_netcdf] netCDF error in ', trim(ctx), ': ', &
                                trim(nf90_strerror(st))
            error stop 1
        end if
    end subroutine nc_check

    integer function nc_open_read(path) result(ncid)
        character(len=*), intent(in) :: path
        call nc_check(nf90_open(trim(path), nf90_nowrite, ncid), 'open '//trim(path))
    end function nc_open_read

    subroutine nc_close(ncid)
        integer, intent(in) :: ncid
        call nc_check(nf90_close(ncid), 'close')
    end subroutine nc_close

    ! First dimension whose name matches any of `names`; returns its length.
    integer function nc_dimlen(ncid, names) result(n)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: names(:)
        integer :: i, did, st
        n = -1
        do i = 1, size(names)
            st = nf90_inq_dimid(ncid, trim(names(i)), did)
            if (st == nf90_noerr) then
                call nc_check(nf90_inquire_dimension(ncid, did, len=n), 'dimlen '//trim(names(i)))
                return
            end if
        end do
        write(error_unit,*) '[mod_io_netcdf] no dimension matched: ', (trim(names(i))//' ', i=1,size(names))
        error stop 1
    end function nc_dimlen

    ! Variable id for the first matching name.
    integer function nc_varid(ncid, names) result(vid)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: names(:)
        integer :: i, st
        do i = 1, size(names)
            st = nf90_inq_varid(ncid, trim(names(i)), vid)
            if (st == nf90_noerr) return
        end do
        write(error_unit,*) '[mod_io_netcdf] no variable matched: ', (trim(names(i))//' ', i=1,size(names))
        error stop 1
    end function nc_varid

    ! Read a 1-D coordinate/axis variable into a real64 array (whatever its on-disk
    ! type — netCDF converts). Used for lon/lat/time axes.
    subroutine nc_get_axis_dp(ncid, names, arr)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: names(:)
        real(real64), intent(out) :: arr(:)
        integer :: vid
        vid = nc_varid(ncid, names)
        call nc_check(nf90_get_var(ncid, vid, arr), 'get_axis '//trim(names(1)))
    end subroutine nc_get_axis_dp

    ! Read a real(4) 2-D lon x lat slice at time index `t_indx` from a (lon,lat,time)
    ! variable. buf must be shape (nLon, nLat). Mirrors FESOM2 nf_get_vara_real with
    ! start=(1,1,t), count=(nLon,nLat,1).
    subroutine nc_get_slice_r4(ncid, names, t_indx, buf)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: names(:)
        integer, intent(in) :: t_indx
        real(real32), intent(out) :: buf(:,:)
        integer :: vid
        vid = nc_varid(ncid, names)
        call nc_check(nf90_get_var(ncid, vid, buf, &
             start=[1, 1, t_indx], count=[size(buf,1), size(buf,2), 1]), &
             'get_slice '//trim(names(1)))
    end subroutine nc_get_slice_r4

    ! Read a real(8) 2-D lon x lat slice at time index `t_indx` from a (lon,lat,time)
    ! variable. buf must be shape (nLon, nLat). The real64 analog of nc_get_slice_r4 —
    ! mirrors FESOM2 read_other_NetCDF (gen_modules_read_NetCDF.F90:102) which reads the
    ! 2D slice straight into a real64 ncdata via nf90_get_var with start=(1,1,itime),
    ! count=(lonlen,latlen,1). On-disk float -> real64 is an exact netCDF conversion.
    subroutine nc_get_slice_dp(ncid, names, t_indx, buf)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: names(:)
        integer, intent(in) :: t_indx
        real(real64), intent(out) :: buf(:,:)
        integer :: vid
        vid = nc_varid(ncid, names)
        call nc_check(nf90_get_var(ncid, vid, buf, &
             start=[1, 1, t_indx], count=[size(buf,1), size(buf,2), 1]), &
             'get_slice_dp '//trim(names(1)))
    end subroutine nc_get_slice_dp

    ! Read a scalar real attribute (e.g. 'missing_value') into real64. Mirrors FESOM2
    ! read_other_NetCDF (gen_modules_read_NetCDF.F90:105) nf90_get_att into a real64 miss;
    ! an on-disk float attribute (1.e30f / -99.f) converts float -> real64 exactly, so the
    ! ==miss equality against the (also float->real64) data values matches the oracle.
    real(real64) function nc_get_att_dp(ncid, varnames, attname, stat) result(val)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: varnames(:)
        character(len=*), intent(in) :: attname
        integer, intent(out), optional :: stat   ! present => return status, don't error-stop
        integer :: vid, st
        val = 0.0_real64
        vid = nc_varid(ncid, varnames)
        st  = nf90_get_att(ncid, vid, trim(attname), val)
        if (present(stat)) then
            ! caller tolerates an absent attribute (e.g. the Sweeney chl file has no missing_value)
            stat = st
        else
            call nc_check(st, 'get_att_dp '//trim(attname))
        end if
    end function nc_get_att_dp

    ! Read a full real(8) (lon,lat,depth) cube from a 3-D variable into buf, whose
    ! shape (nLon, nLat, nDepth) must match the on-disk extents. Mirrors the FESOM2
    ! 3D-IC read (gen_ic3d.F90 getcoeffld:369, nf_get_vara_double, start=1, edges=full)
    ! into a real(WP) ncdata buffer — the bytes are identical (netCDF converts the
    ! on-disk type, here double, to real64 deterministically; no auto fill/scale in
    ! the base library, so NaN land sentinels pass through verbatim for the caller's
    ! ieee_is_nan mask).
    subroutine nc_get_var3d_dp(ncid, names, buf)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: names(:)
        real(real64), intent(out) :: buf(:,:,:)
        integer :: vid
        vid = nc_varid(ncid, names)
        call nc_check(nf90_get_var(ncid, vid, buf, &
             start=[1, 1, 1], count=[size(buf,1), size(buf,2), size(buf,3)]), &
             'get_var3d '//trim(names(1)))
    end subroutine nc_get_var3d_dp

    ! Read a text attribute (e.g. 'calendar') from a variable; returns '' if absent.
    function nc_get_att_text(ncid, varnames, attname) result(txt)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: varnames(:)
        character(len=*), intent(in) :: attname
        character(len=64) :: txt
        integer :: vid, st, alen, k
        txt = ''
        vid = nc_varid(ncid, varnames)
        st = nf90_inquire_attribute(ncid, vid, trim(attname), len=alen)
        if (st /= nf90_noerr) return
        st = nf90_get_att(ncid, vid, trim(attname), txt)
        if (st /= nf90_noerr) then; txt = ''; return; end if
        ! netCDF C-string attributes can carry a trailing NUL (and other control bytes); blank
        ! them so trim()/string compares work. This was the JRA55 trap: calendar='gregorian'+char(0)
        ! made forcing_julday's `trim(calendar)=='gregorian'` test fail -> noleap branch (CORE2's
        ! 'NOLEAP' was immune since it falls to the same else branch regardless).
        do k = 1, len(txt)
            if (iachar(txt(k:k)) < 32) txt(k:k) = ' '
        end do
    end function nc_get_att_text

end module mod_io_netcdf
