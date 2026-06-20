program fesom_ncsmoke
    ! Smoke test: validate FESOM3 can read the CORE2 forcing stub via mod_io_netcdf,
    ! and that the bytes match ncdump. Expected (u_10.1948.nc):
    !   LON(1)=0  LON(2)=1.875   LAT(1)=-88.542   TIME=0,6,12,18,24
    !   nLon=192 nLat=94 nTime=5
    use, intrinsic :: iso_fortran_env, only: real32, real64
    use mod_io_netcdf
    implicit none
    character(len=*), parameter :: F = &
        '/home/a/a270088/port2/fesom2/test/input/global/u_10.1948.nc'
    integer :: ncid, nlon, nlat, ntime
    real(real64), allocatable :: lon(:), lat(:), tim(:)
    real(real32), allocatable :: slice(:,:)
    character(len=64) :: cal

    ncid  = nc_open_read(F)
    nlon  = nc_dimlen(ncid, ['LON','lon'])
    nlat  = nc_dimlen(ncid, ['LAT','lat'])
    ntime = nc_dimlen(ncid, ['TIME','time'])
    print '(A,3I6)', 'dims (nLon,nLat,nTime): ', nlon, nlat, ntime

    allocate(lon(nlon), lat(nlat), tim(ntime), slice(nlon,nlat))
    call nc_get_axis_dp(ncid, ['LON','lon'], lon)
    call nc_get_axis_dp(ncid, ['LAT','lat'], lat)
    call nc_get_axis_dp(ncid, ['TIME','time'], tim)
    cal = nc_get_att_text(ncid, ['TIME','time'], 'calendar')
    call nc_get_slice_r4(ncid, ['U_10_MOD','uas     '], 1, slice)

    print '(A,2F12.5)', 'LON(1:2)   = ', lon(1), lon(2)
    print '(A,2F12.5)', 'LAT(1),LAT(end) = ', lat(1), lat(nlat)
    print '(A,5F6.1)',  'TIME       = ', tim(1:ntime)
    print '(A,A)',      'calendar   = ', trim(cal)
    print '(A,2F12.6)', 'U_10(1,1),U_10(2,1) = ', slice(1,1), slice(2,1)
    call nc_close(ncid)
    print '(A)', 'NCSMOKE OK'
end program fesom_ncsmoke
