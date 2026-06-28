module mod_io_means
    ! FESOM3 field output as per-variable-per-year Zarr v2 stores (M9 Stage 2, Tasks 2.1 + 2.2).
    !
    ! Each registered field gets its OWN store, split by year: <out_dir>/<name>.fesom.<YYYY>.zarr,
    ! holding the data var (time, nod2) + a CF `time` coord + embedded `lon`/`lat` node coords. This
    ! mirrors FESOM2's per-variable-per-year mean files (<name>.fesom.<year>.nc), but as Zarr.
    !
    ! Reuses the proven Stage-1 stack verbatim:
    !   - mod_io_decomp : canonical-order MPI_Alltoallv redistribution -> writer subset (no rank-0
    !                     gather), so the store is partition-INDEPENDENT (dist_2 == dist_8).
    !   - mod_io_zarr   : the Zarr v2 writer; store-create ordering = rank 0 creates store + defines
    !                     all arrays -> barrier -> writers write their chunks.
    !
    ! NEW vs Stage 1 (mesh.diag is write-once): a GROWING time dimension. With chunk_time pinned to 1
    ! (v1), one output record == one fresh chunk file per writer, so "append" = write the new data
    ! chunks [t,c] + the time-coord chunk [t] + bump each .zarray shape[0] to t+1 (zarr_rewrite_zarray)
    ! — NO read-modify-write of a partial time-chunk (chunk_time>1 is deferred to Task 2.6).
    !
    ! MEAN vs SNAPSHOT (Task 2.2) — transcribed from io_meandata.F90:update_means/compute_means:
    !   - every step `means_accumulate` does, IN THE OUTPUT PRECISION, mean: acc += value, count++ ;
    !     snapshot: acc = value, count = 1  (io_meandata.F90:2107/2127/2142).
    !   - at output `means_write` writes acc / count (divide in the output precision: real64 for <f8,
    !     real32 for <f4 — io_meandata.F90:2335/2353), then zeroes acc + count for the next interval.
    !   So float32 means accumulate AND divide in float32, byte-matching FESOM2's r4 mean stream.
    !
    ! Field types: node SCALARS 2-D (ssh/sst/...) + 3-D (T/S/w), and node VECTOR pairs (unod/vnod, the
    ! velocity at nodes — FESOM2's default velocity output, dynamics%uvnode) which are r2g-rotated to
    ! geographic at write time (Task 2.5, means_define_vector3d). Element-based u/v (dynamics%uv) is the
    ! alternative FESOM2 offers but never enables by default — left as a future elem-decomp addition.
    !
    ! Caller protocol (all ranks, same order):
    !   means_init(io, dir, mesh[, partit])              once
    !   means_define_node2d(io, name, long, units, ...)  per field, once (mean=.true. for a mean stream)
    !   ... each step ...
    !   means_accumulate(io, 'ssh', eta_n(1:myDim)) ...  EVERY step (mean sums; snapshot overwrites)
    !   ... when an output interval is due ...
    !   means_begin(io, year, time_sec)                  opens per-year stores on year change
    !   means_write(io)                                  divide+write all fields; reset; t++
    use mpi
    use, intrinsic :: iso_fortran_env, only: real64, real32
    use mod_precision,   only: WP
    use mod_constants,   only: rad
    use mod_mesh,        only: t_mesh
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: is_multirank, local_dims
    use mod_mesh_rotate, only: vector_r2g
    use mod_io_zarr
    use mod_io_decomp
    implicit none
    private
    public :: t_io_means, means_init, means_define_node2d, means_define_node3d, means_define_vector3d
    public :: means_accumulate, means_begin, means_write, means_finalize

    interface means_accumulate
        module procedure means_accumulate_2d, means_accumulate_3d
    end interface means_accumulate

    integer, parameter :: MEANS_MAXF = 128
    ! netCDF standard fill values (== FESOM2 io_meandata.F90 NC_FILL_*); below-bottom levels get this,
    ! xarray masks _FillValue -> NaN (CF). Plan Task 2.4: nlevels-based masking (cleaner than FESOM2's
    ! value-based abs(acc)<1e-30 quirk; valid-level VALUES still byte-match FESOM2).
    real(real64), parameter :: NC_FILL_DOUBLE = 9.9692099683868690e+36_real64
    real(real32), parameter :: NC_FILL_FLOAT  = 9.9692099683868690e+36_real32

    ! one registered output field (its per-year store + array handles + accumulator)
    type :: t_mean_field
        character(len=64)  :: name = '', long_name = '', units = '', std = ''
        character(len=8)   :: dtype = '<f4'        ! float32 default (Task 2.6 exposes the knob)
        logical            :: is_mean = .false.     ! mean stream (Task 2.2) vs snapshot (Task 2.1)
        integer            :: ndim = 2              ! 2 = node scalar (time,nod2); 3 = (time,nz,nod2)
        integer            :: nlev = 1              ! vertical size (3-D): nl-1 (layers) or nl (levels)
        integer            :: voff = 0              ! valid levels per node = nlevels_nod2D - voff
        character(len=8)   :: vdim = ''             ! 'nz1' (layers) / 'nz' (levels)
        type(t_zarr_store) :: store                 ! current-year store
        type(t_zarr_array) :: a_data, a_time, a_lon, a_lat, a_vert
        logical            :: open  = .false.       ! store created for the current year
        ! owned accumulator, in the OUTPUT precision (FESOM2-faithful) — one set per field
        real(real64), allocatable :: acc8(:),    acc8_3d(:,:)
        real(real32), allocatable :: acc4(:),    acc4_3d(:,:)
        integer                   :: count = 0
        ! vector pairing (Task 2.5): the (x,y) components of a velocity/wind pair are r2g-rotated
        ! TOGETHER at write time (io_meandata.F90:io_r2g). is_vec_x marks the x-component (drives the
        ! paired write); is_vec_y the y-component (its data is written by its partner, skipped in the
        ! per-field data loop). vec_partner = the other component's field index.
        logical :: is_vec_x = .false., is_vec_y = .false.
        integer :: vec_partner = 0
    end type t_mean_field

    type :: t_io_means
        character(len=512) :: out_dir  = '.'
        character(len=32)  :: calendar = 'standard'
        integer            :: nf = 0
        type(t_mean_field) :: f(MEANS_MAXF)
        type(t_io_decomp)  :: Dn                     ! node decomp (node scalars)
        integer            :: nNodO = 0
        integer            :: mype = 0, comm = MPI_COMM_SELF, npes = 1
        logical            :: mr = .false.
        integer            :: cur_year = -2000000000 ! sentinel "no year opened yet"
        integer            :: t = 0                  ! record index within the current year (0-based)
        real(real64)       :: cur_time = 0.0_real64  ! time-coord value for the in-flight record
        integer            :: nl = 0                 ! mesh%nl (level count)
        real(WP), allocatable :: lon_owned(:), lat_owned(:)  ! cached owned node coords (deg, geographic)
        real(WP), allocatable :: rlon_owned(:), rlat_owned(:)! cached owned ROTATED node coords (rad) for r2g
        integer,  allocatable :: nlev_owned(:)               ! cached owned nlevels_nod2D (3-D mask)
        real(WP), allocatable :: depth_nz(:), depth_nz1(:)   ! cached vertical coords (-zbar / -Z)
        ! vector output frame (Task 2.5): .true. => r2g-rotate vector pairs to GEOGRAPHIC before write
        ! (= FESOM2 vec_autorotate=.true., the production default); .false. => write NATIVE rotated-mesh
        ! components. FESOM3 default = geographic (the scientifically-useful frame).
        logical               :: autorotate = .true.
    end type t_io_means

contains

    ! ----------------------------------------------------------------- setup / registration

    subroutine means_init(io, out_dir, mesh, partit, chunk_horiz, n_writers, calendar, vec_frame)
        type(t_io_means), intent(out)          :: io
        character(len=*), intent(in)           :: out_dir
        type(t_mesh),     intent(in)           :: mesh
        type(t_partit),   intent(in), optional :: partit
        integer,          intent(in), optional :: chunk_horiz, n_writers
        character(len=*), intent(in), optional :: calendar
        character(len=*), intent(in), optional :: vec_frame    ! 'geographic' (default) | 'native'
        integer :: C, nw, i, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF
        character(len=16) :: frame
        io%out_dir  = out_dir
        io%nf       = 0
        io%cur_year = -2000000000
        io%t        = 0
        io%calendar = 'standard'
        if (present(calendar)) io%calendar = calendar
        ! vector output frame: default geographic (autorotate); arg overrides default; env overrides arg.
        io%autorotate = .true.
        if (present(vec_frame)) io%autorotate = (trim(vec_frame) /= 'native')
        frame = ''
        call get_environment_variable('FESOM3_VEC_FRAME', frame)
        if (len_trim(frame) > 0) io%autorotate = (trim(frame) /= 'native')
        io%mr   = is_multirank(partit)
        io%mype = 0; io%comm = MPI_COMM_SELF; io%npes = 1
        if (present(partit)) then
            io%mype = partit%mype; io%comm = partit%MPI_COMM_FESOM; io%npes = partit%npes
        end if
        C = 500000; if (present(chunk_horiz)) C = chunk_horiz
        nw = 0;     if (present(n_writers))   nw = n_writers
        call read_env_int('FESOM3_CHUNK_HORIZ', C)
        call read_env_int('FESOM3_N_WRITERS',  nw)
        call local_dims(mesh, partit, io%nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF)
        call decomp_init_entity(io%Dn, C, nw, DECOMP_NODE, mesh, partit)
        io%nl = mesh%nl
        ! cache owned node lon/lat (deg) + nlevels (3-D mask) — embedded into each per-year store.
        ! Also cache the ROTATED node coords (radians) that vector_r2g needs at write time (flag_coord=0,
        ! exactly as io_meandata.F90:io_r2g passes mesh%coord_nod2D for node-based vectors).
        allocate(io%lon_owned(max(1,io%nNodO)),  io%lat_owned(max(1,io%nNodO)), &
                 io%rlon_owned(max(1,io%nNodO)), io%rlat_owned(max(1,io%nNodO)), &
                 io%nlev_owned(max(1,io%nNodO)))
        do i = 1, io%nNodO
            io%lon_owned(i)  = real(mesh%geo_coord_nod2D(1,i)/rad, WP)
            io%lat_owned(i)  = real(mesh%geo_coord_nod2D(2,i)/rad, WP)
            io%rlon_owned(i) = real(mesh%coord_nod2D(1,i), WP)
            io%rlat_owned(i) = real(mesh%coord_nod2D(2,i), WP)
            io%nlev_owned(i) = mesh%nlevels_nod2D(i)
        end do
        ! vertical coords (CF positive-down depths): nz = -zbar(1:nl), nz1 = -Z(1:nl-1)
        allocate(io%depth_nz(io%nl), io%depth_nz1(io%nl-1))
        do i = 1, io%nl;   io%depth_nz(i)  = real(-mesh%zbar(i), WP); end do
        do i = 1, io%nl-1; io%depth_nz1(i) = real(-mesh%Z(i),    WP); end do
    end subroutine means_init

    subroutine means_define_node2d(io, name, long_name, units, std, precision, mean)
        type(t_io_means), intent(inout)        :: io
        character(len=*), intent(in)           :: name, long_name, units
        character(len=*), intent(in), optional :: std, precision
        logical,          intent(in), optional :: mean
        integer :: k, n
        call zarr_check(io%nf < MEANS_MAXF, 'means_define_node2d: too many fields')
        io%nf = io%nf + 1
        k = io%nf
        io%f(k)%name      = name
        io%f(k)%long_name = long_name
        io%f(k)%units     = units
        io%f(k)%std       = ''
        if (present(std)) io%f(k)%std = std
        io%f(k)%dtype = '<f4'
        if (present(precision)) then
            if (trim(precision) == 'double' .or. trim(precision) == '<f8') io%f(k)%dtype = '<f8'
        end if
        io%f(k)%is_mean = .false.
        if (present(mean)) io%f(k)%is_mean = mean
        io%f(k)%open  = .false.
        io%f(k)%count = 0
        io%f(k)%ndim = 2
        n = max(1, io%nNodO)
        if (io%f(k)%dtype == '<f8') then
            allocate(io%f(k)%acc8(n)); io%f(k)%acc8 = 0.0_real64
        else
            allocate(io%f(k)%acc4(n)); io%f(k)%acc4 = 0.0_real32
        end if
    end subroutine means_define_node2d

    ! Register a 3-D node field. on_full_levels=.true. => nl levels (w; vdim 'nz', valid 1..nlevels);
    ! .false. => nl-1 layers (T/S; vdim 'nz1', valid 1..nlevels-1). Below-bottom => NC_FILL (CF mask).
    subroutine means_define_node3d(io, name, long_name, units, on_full_levels, std, precision, mean)
        type(t_io_means), intent(inout)        :: io
        character(len=*), intent(in)           :: name, long_name, units
        logical,          intent(in)           :: on_full_levels
        character(len=*), intent(in), optional :: std, precision
        logical,          intent(in), optional :: mean
        integer :: k, n
        call zarr_check(io%nf < MEANS_MAXF, 'means_define_node3d: too many fields')
        io%nf = io%nf + 1
        k = io%nf
        io%f(k)%name = name; io%f(k)%long_name = long_name; io%f(k)%units = units
        io%f(k)%std = ''; if (present(std)) io%f(k)%std = std
        io%f(k)%dtype = '<f4'
        if (present(precision)) then
            if (trim(precision) == 'double' .or. trim(precision) == '<f8') io%f(k)%dtype = '<f8'
        end if
        io%f(k)%is_mean = .false.; if (present(mean)) io%f(k)%is_mean = mean
        io%f(k)%ndim = 3
        io%f(k)%open = .false.; io%f(k)%count = 0
        if (on_full_levels) then
            io%f(k)%nlev = io%nl;   io%f(k)%voff = 0; io%f(k)%vdim = 'nz'
        else
            io%f(k)%nlev = io%nl-1; io%f(k)%voff = 1; io%f(k)%vdim = 'nz1'
        end if
        n = max(1, io%nNodO)
        if (io%f(k)%dtype == '<f8') then
            allocate(io%f(k)%acc8_3d(io%f(k)%nlev, n)); io%f(k)%acc8_3d = 0.0_real64
        else
            allocate(io%f(k)%acc4_3d(io%f(k)%nlev, n)); io%f(k)%acc4_3d = 0.0_real32
        end if
    end subroutine means_define_node3d

    ! Register a 3-D node VECTOR pair (Task 2.5): two node3d fields (x=zonal, y=meridional) linked so
    ! means_write r2g-rotates them TOGETHER to geographic (when autorotate) — the io_meandata io_r2g
    ! analog for unod/vnod. Each component is accumulated independently (means_accumulate as usual); the
    ! pairing matters only at write time. on_full_levels follows means_define_node3d (false => nl-1/nz1).
    subroutine means_define_vector3d(io, name_x, name_y, long_x, long_y, units, on_full_levels, &
                                     std_x, std_y, precision, mean)
        type(t_io_means), intent(inout)        :: io
        character(len=*), intent(in)           :: name_x, name_y, long_x, long_y, units
        logical,          intent(in)           :: on_full_levels
        character(len=*), intent(in), optional :: std_x, std_y, precision
        logical,          intent(in), optional :: mean
        integer :: kx, ky
        call means_define_node3d(io, name_x, long_x, units, on_full_levels, &
                                 std=std_x, precision=precision, mean=mean)
        kx = io%nf
        call means_define_node3d(io, name_y, long_y, units, on_full_levels, &
                                 std=std_y, precision=precision, mean=mean)
        ky = io%nf
        io%f(kx)%is_vec_x = .true.; io%f(kx)%vec_partner = ky
        io%f(ky)%is_vec_y = .true.; io%f(ky)%vec_partner = kx
    end subroutine means_define_vector3d

    ! ----------------------------------------------------------------- accumulate / write

    ! Accumulate a 2-D node field's owned value for THIS step (in the output precision). mean =>
    ! running sum + count; snapshot => overwrite (count=1). Call every step (FESOM2 update_means).
    subroutine means_accumulate_2d(io, name, owned)
        type(t_io_means), intent(inout) :: io
        character(len=*), intent(in)    :: name
        real(WP),         intent(in)    :: owned(:)
        integer :: k, i
        k = find_field(io, name)
        call zarr_check(k > 0, 'means_accumulate: unknown field '//trim(name))
        call zarr_check(io%f(k)%ndim == 2, 'means_accumulate: 2-D call on 3-D field '//trim(name))
        if (io%f(k)%is_mean) then
            if (io%f(k)%dtype == '<f8') then
                do i = 1, io%nNodO; io%f(k)%acc8(i) = io%f(k)%acc8(i) + real(owned(i), real64); end do
            else
                do i = 1, io%nNodO; io%f(k)%acc4(i) = io%f(k)%acc4(i) + real(owned(i), real32); end do
            end if
            io%f(k)%count = io%f(k)%count + 1
        else
            if (io%f(k)%dtype == '<f8') then
                do i = 1, io%nNodO; io%f(k)%acc8(i) = real(owned(i), real64); end do
            else
                do i = 1, io%nNodO; io%f(k)%acc4(i) = real(owned(i), real32); end do
            end if
            io%f(k)%count = 1
        end if
    end subroutine means_accumulate_2d

    ! Accumulate a 3-D node field owned(1:nlev, 1:nNodO) for THIS step (output precision).
    subroutine means_accumulate_3d(io, name, owned)
        type(t_io_means), intent(inout) :: io
        character(len=*), intent(in)    :: name
        real(WP),         intent(in)    :: owned(:,:)
        integer :: k, i, L, nlev
        k = find_field(io, name)
        call zarr_check(k > 0, 'means_accumulate: unknown field '//trim(name))
        call zarr_check(io%f(k)%ndim == 3, 'means_accumulate: 3-D call on 2-D field '//trim(name))
        nlev = io%f(k)%nlev
        if (io%f(k)%is_mean) then
            if (io%f(k)%dtype == '<f8') then
                do i = 1, io%nNodO; do L = 1, nlev
                    io%f(k)%acc8_3d(L,i) = io%f(k)%acc8_3d(L,i) + real(owned(L,i), real64); end do; end do
            else
                do i = 1, io%nNodO; do L = 1, nlev
                    io%f(k)%acc4_3d(L,i) = io%f(k)%acc4_3d(L,i) + real(owned(L,i), real32); end do; end do
            end if
            io%f(k)%count = io%f(k)%count + 1
        else
            if (io%f(k)%dtype == '<f8') then
                do i = 1, io%nNodO; do L = 1, nlev
                    io%f(k)%acc8_3d(L,i) = real(owned(L,i), real64); end do; end do
            else
                do i = 1, io%nNodO; do L = 1, nlev
                    io%f(k)%acc4_3d(L,i) = real(owned(L,i), real32); end do; end do
            end if
            io%f(k)%count = 1
        end if
    end subroutine means_accumulate_3d

    ! Start a record. On a year change, (re)create every field's per-year store + embed lon/lat.
    subroutine means_begin(io, year, time_sec)
        type(t_io_means), intent(inout) :: io
        integer,          intent(in)    :: year
        real(real64),     intent(in)    :: time_sec
        if (year /= io%cur_year) then
            call open_year(io, year)
            io%cur_year = year
            io%t = 0
        end if
        io%cur_time = time_sec
    end subroutine means_begin

    ! Write the in-flight record for ALL fields: divide each accumulator by its count (in the output
    ! precision), redistribute to canonical order, writers write chunks [t,c]; rank 0 appends the
    ! time-coord chunk [t] and bumps each data/time .zarray shape[0]; then reset accumulators, t++.
    subroutine means_write(io)
        type(t_io_means), intent(inout) :: io
        integer :: k, ierr
        ! --- data (collective redistribute per field; writers write their chunks) ---
        do k = 1, io%nf
            if (io%f(k)%is_vec_y) cycle                 ! its data is written by the vec_x partner (paired r2g)
            if (io%f(k)%is_vec_x) then
                call write_vector_3d(k, io%f(k)%vec_partner)
            else if (io%f(k)%ndim == 2) then
                call write_data_2d(k)
            else
                call write_data_3d(k)
            end if
        end do
        if (io%mr) call MPI_Barrier(io%comm, ierr)   ! all writers flushed record io%t
        ! --- metadata (rank 0): time coord chunk [t] + shape bumps ---
        if (io%mype == 0) then
            do k = 1, io%nf
                call zarr_write_chunk(io%f(k)%store, io%f(k)%a_time, [io%t], [io%cur_time])
                io%f(k)%a_data%dims(1) = io%t + 1
                io%f(k)%a_time%dims(1) = io%t + 1
                call zarr_rewrite_zarray(io%f(k)%store, io%f(k)%a_data)
                call zarr_rewrite_zarray(io%f(k)%store, io%f(k)%a_time)
            end do
        end if
        ! --- reset accumulators for the next interval (io_meandata clean_meanarrays) ---
        do k = 1, io%nf
            if (allocated(io%f(k)%acc8))    io%f(k)%acc8    = 0.0_real64
            if (allocated(io%f(k)%acc4))    io%f(k)%acc4    = 0.0_real32
            if (allocated(io%f(k)%acc8_3d)) io%f(k)%acc8_3d = 0.0_real64
            if (allocated(io%f(k)%acc4_3d)) io%f(k)%acc4_3d = 0.0_real32
            io%f(k)%count = 0
        end do
        io%t = io%t + 1
        if (io%mr) call MPI_Barrier(io%comm, ierr)

    contains

        subroutine write_data_2d(k)
            integer, intent(in) :: k
            real(WP), allocatable :: owned_mean(:), buf(:), chunk2(:,:)
            integer :: i, c, lo, cnt
            allocate(owned_mean(max(1,io%nNodO)), buf(max(1, io%Dn%w_nbuf)), chunk2(1, io%Dn%C))
            cnt = max(1, io%f(k)%count)
            if (io%f(k)%dtype == '<f8') then
                do i = 1, io%nNodO; owned_mean(i) = io%f(k)%acc8(i) / real(cnt, real64); end do
            else
                do i = 1, io%nNodO; owned_mean(i) = real(io%f(k)%acc4(i) / real(cnt, real32), WP); end do
            end if
            call decomp_redistribute(io%Dn, owned_mean(1:io%nNodO), buf, 0.0_WP)
            do c = io%Dn%w_first_chunk, io%Dn%w_last_chunk
                lo = (c - io%Dn%w_first_chunk)*io%Dn%C + 1
                chunk2(1, 1:io%Dn%C) = buf(lo:lo + io%Dn%C - 1)
                call zarr_write_chunk(io%f(k)%store, io%f(k)%a_data, [io%t, c], chunk2)
            end do
        end subroutine write_data_2d

        ! 3-D: divide (output precision), mask below-bottom (L > nlevels_nod2D - voff) to NC_FILL,
        ! redistribute level-by-level, writers write chunk [t, 0, c] = (1, nlev, C).
        subroutine write_data_3d(k)
            integer, intent(in) :: k
            real(WP), allocatable :: mean3(:,:), buf3(:,:), chunk3(:,:,:)
            integer :: i, L, c, lo, cnt, nlev, nvalid
            nlev = io%f(k)%nlev
            cnt  = max(1, io%f(k)%count)
            allocate(mean3(nlev, max(1,io%nNodO)), buf3(nlev, max(1, io%Dn%w_nbuf)), &
                     chunk3(1, nlev, io%Dn%C))
            do i = 1, io%nNodO
                nvalid = io%nlev_owned(i) - io%f(k)%voff
                do L = 1, nlev
                    if (L <= nvalid) then
                        if (io%f(k)%dtype == '<f8') then
                            mean3(L,i) = io%f(k)%acc8_3d(L,i) / real(cnt, real64)
                        else
                            mean3(L,i) = real(io%f(k)%acc4_3d(L,i) / real(cnt, real32), WP)
                        end if
                    else
                        mean3(L,i) = real(NC_FILL_DOUBLE, WP)
                    end if
                end do
            end do
            call decomp_redistribute(io%Dn, mean3(1:nlev,1:io%nNodO), buf3, real(NC_FILL_DOUBLE, WP))
            do c = io%Dn%w_first_chunk, io%Dn%w_last_chunk
                lo = (c - io%Dn%w_first_chunk)*io%Dn%C + 1
                chunk3(1, 1:nlev, 1:io%Dn%C) = buf3(1:nlev, lo:lo + io%Dn%C - 1)
                call zarr_write_chunk(io%f(k)%store, io%f(k)%a_data, [io%t, 0, c], chunk3)
            end do
        end subroutine write_data_3d

        ! 3-D vector pair (kx=zonal, ky=meridional). FESOM2 ORDER (io_meandata.F90): io_r2g rotates the
        ! ACCUMULATED SUM in place (L2265) BEFORE compute_means divides by addcounter (L2335) — and
        ! rotate-then-divide /= divide-then-rotate in floating point, so we match it: rotate the SUM,
        ! THEN divide. <f8: rotate the r8 sum, divide in r8. <f4: promote the r4 sum to WP(=r8), rotate,
        ! demote to r4 (io_r2g r4 branch L3028-3033), divide in r4. Below-bottom (nlevels mask) -> NC_FILL,
        ! not rotated. native frame (autorotate=.false.) skips rotation -> identical to a plain field.
        subroutine write_vector_3d(kx, ky)
            integer, intent(in) :: kx, ky
            real(WP), allocatable :: mx(:,:), my(:,:), bx(:,:), by(:,:), chunk3(:,:,:)
            integer :: i, L, c, lo, cnt, nlev, nvalid
            real(real64) :: sx8, sy8
            real(real32) :: sx4, sy4
            real(WP)     :: tx, ty
            nlev = io%f(kx)%nlev
            cnt  = max(1, io%f(kx)%count)
            allocate(mx(nlev, max(1,io%nNodO)), my(nlev, max(1,io%nNodO)), &
                     bx(nlev, max(1, io%Dn%w_nbuf)), by(nlev, max(1, io%Dn%w_nbuf)), &
                     chunk3(1, nlev, io%Dn%C))
            do i = 1, io%nNodO
                nvalid = io%nlev_owned(i) - io%f(kx)%voff
                do L = 1, nlev
                    if (L > nvalid) then
                        mx(L,i) = real(NC_FILL_DOUBLE, WP); my(L,i) = real(NC_FILL_DOUBLE, WP)
                        cycle
                    end if
                    if (io%f(kx)%dtype == '<f8') then
                        sx8 = io%f(kx)%acc8_3d(L,i); sy8 = io%f(ky)%acc8_3d(L,i)
                        if (io%autorotate) then
                            tx = real(sx8, WP); ty = real(sy8, WP)
                            call vector_r2g(tx, ty, io%rlon_owned(i), io%rlat_owned(i), 0)
                            sx8 = real(tx, real64); sy8 = real(ty, real64)
                        end if
                        mx(L,i) = real(sx8 / real(cnt, real64), WP)
                        my(L,i) = real(sy8 / real(cnt, real64), WP)
                    else
                        sx4 = io%f(kx)%acc4_3d(L,i); sy4 = io%f(ky)%acc4_3d(L,i)
                        if (io%autorotate) then
                            tx = real(sx4, WP); ty = real(sy4, WP)
                            call vector_r2g(tx, ty, io%rlon_owned(i), io%rlat_owned(i), 0)
                            sx4 = real(tx, real32); sy4 = real(ty, real32)
                        end if
                        mx(L,i) = real(sx4 / real(cnt, real32), WP)
                        my(L,i) = real(sy4 / real(cnt, real32), WP)
                    end if
                end do
            end do
            call decomp_redistribute(io%Dn, mx(1:nlev,1:io%nNodO), bx, real(NC_FILL_DOUBLE, WP))
            call decomp_redistribute(io%Dn, my(1:nlev,1:io%nNodO), by, real(NC_FILL_DOUBLE, WP))
            do c = io%Dn%w_first_chunk, io%Dn%w_last_chunk
                lo = (c - io%Dn%w_first_chunk)*io%Dn%C + 1
                chunk3(1, 1:nlev, 1:io%Dn%C) = bx(1:nlev, lo:lo + io%Dn%C - 1)
                call zarr_write_chunk(io%f(kx)%store, io%f(kx)%a_data, [io%t, 0, c], chunk3)
                chunk3(1, 1:nlev, 1:io%Dn%C) = by(1:nlev, lo:lo + io%Dn%C - 1)
                call zarr_write_chunk(io%f(ky)%store, io%f(ky)%a_data, [io%t, 0, c], chunk3)
            end do
        end subroutine write_vector_3d

    end subroutine means_write

    ! End of run. v1 keeps each <var>/.zarray current every record, so stores open with
    ! consolidated=False already; consolidated .zmetadata for field stores is a Task 2.6 polish.
    subroutine means_finalize(io)
        type(t_io_means), intent(inout) :: io
        integer :: ierr
        if (io%mr) call MPI_Barrier(io%comm, ierr)
    end subroutine means_finalize

    ! ----------------------------------------------------------------- internals

    ! Create + define every field's store for `year`, then writers embed the static lon/lat.
    subroutine open_year(io, year)
        type(t_io_means), intent(inout) :: io
        integer,          intent(in)    :: year
        integer :: k, ierr
        character(len=4)              :: cy
        character(len=:), allocatable :: tunits
        write(cy, '(i4.4)') year
        tunits = 'seconds since '//cy//'-01-01 00:00:00'
        do k = 1, io%nf
            io%f(k)%store%path = trim(io%out_dir)//'/'//trim(io%f(k)%name)//'.fesom.'//cy//'.zarr'
            call def_field_store(io%f(k), io%Dn%N, io%Dn%C, io%mype, tunits, trim(io%calendar))
        end do
        if (io%mr) call MPI_Barrier(io%comm, ierr)
        do k = 1, io%nf
            call put_static(io, io%f(k)%store, io%f(k)%a_lon, io%lon_owned)
            call put_static(io, io%f(k)%store, io%f(k)%a_lat, io%lat_owned)
            ! vertical coord is global/replicated -> rank 0 writes it whole (like mesh.diag nz/nz1)
            if (io%mype == 0 .and. io%f(k)%ndim == 3) then
                if (trim(io%f(k)%vdim) == 'nz') then
                    call zarr_write_whole(io%f(k)%store, io%f(k)%a_vert, io%depth_nz(1:io%nl))
                else
                    call zarr_write_whole(io%f(k)%store, io%f(k)%a_vert, io%depth_nz1(1:io%nl-1))
                end if
            end if
            io%f(k)%open = .true.
        end do
        if (io%mr) call MPI_Barrier(io%comm, ierr)
    end subroutine open_year

    ! Define one field store: data (time,nod2) grows from shape[0]=0; time (time,); lon/lat (nod2,).
    ! Array handles init on all ranks; dirs/.zarray/.zattrs written on rank 0 (store-create ordering).
    subroutine def_field_store(f, N, C, mype, tunits, calendar)
        type(t_mean_field), intent(inout) :: f
        integer,            intent(in)    :: N, C, mype
        character(len=*),   intent(in)    :: tunits, calendar
        type(t_zarr_attrs) :: gat, at
        ! handles (all ranks). 2-D: fill_value:null (a valid 0 isn't NaN-masked). 3-D: NC_FILL so
        ! xarray masks below-bottom -> NaN (CF). time/lon/lat/vert coords: null fill.
        if (f%ndim == 2) then
            call zarr_array_init(f%a_data, trim(f%name), [0, N], [1, C], f%dtype, has_fill=.false.)
        else
            call zarr_array_init(f%a_data, trim(f%name), [0, f%nlev, N], [1, f%nlev, C], f%dtype, &
                                 fill=NC_FILL_DOUBLE, has_fill=.true.)
            call zarr_array_init(f%a_vert, trim(f%vdim), [f%nlev], [f%nlev], '<f8', has_fill=.false.)
        end if
        call zarr_array_init(f%a_time, 'time', [0], [1], '<f8', has_fill=.false.)
        call zarr_array_init(f%a_lon,  'lon',  [N], [C], '<f8', has_fill=.false.)
        call zarr_array_init(f%a_lat,  'lat',  [N], [C], '<f8', has_fill=.false.)
        if (mype /= 0) return
        call zarr_attrs_init(gat)
        call zattr_str(gat, 'Conventions', 'CF-1.8')
        call zattr_str(gat, 'description', 'FESOM3 field output (Zarr v2, M9)')
        call zarr_create_store(f%store, trim(f%store%path), gat)
        ! data variable
        call zarr_attrs_init(at)
        if (f%ndim == 2) then
            call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: 'time', 'nod2'])
        else
            call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: 'time', f%vdim, 'nod2'])
        end if
        call zattr_str(at, 'long_name', trim(f%long_name))
        if (len_trim(f%units) > 0) call zattr_str(at, 'units', trim(f%units))
        if (len_trim(f%std)   > 0) call zattr_str(at, 'standard_name', trim(f%std))
        if (f%is_mean) then
            call zattr_str(at, 'cell_methods', 'time: mean')
        else
            call zattr_str(at, 'cell_methods', 'time: point')
        end if
        call zattr_str(at, 'coordinates', 'lon lat')
        call zarr_define_array(f%store, f%a_data, at)
        ! vertical coordinate (3-D only)
        if (f%ndim == 3) then
            call zarr_attrs_init(at)
            call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: f%vdim])
            call zattr_str(at, 'long_name', 'depth'); call zattr_str(at, 'units', 'meters')
            call zattr_str(at, 'positive', 'down')
            call zarr_define_array(f%store, f%a_vert, at)
        end if
        ! time coordinate (CF)
        call zarr_attrs_init(at)
        call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: 'time'])
        call zattr_str(at, 'long_name', 'time')
        call zattr_str(at, 'units', tunits)
        call zattr_str(at, 'calendar', calendar)
        call zarr_define_array(f%store, f%a_time, at)
        ! lon / lat node coords
        call zarr_attrs_init(at)
        call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: 'nod2'])
        call zattr_str(at, 'long_name', 'longitude'); call zattr_str(at, 'units', 'degrees_east')
        call zattr_str(at, 'standard_name', 'longitude')
        call zarr_define_array(f%store, f%a_lon, at)
        call zarr_attrs_init(at)
        call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: 'nod2'])
        call zattr_str(at, 'long_name', 'latitude'); call zattr_str(at, 'units', 'degrees_north')
        call zattr_str(at, 'standard_name', 'latitude')
        call zarr_define_array(f%store, f%a_lat, at)
    end subroutine def_field_store

    ! Redistribute a static 1-D node field and let each writer write its chunks [c] (lon/lat embed).
    subroutine put_static(io, store, arr, owned)
        type(t_io_means),   intent(in) :: io
        type(t_zarr_store), intent(in) :: store
        type(t_zarr_array), intent(in) :: arr
        real(WP),           intent(in) :: owned(:)
        real(WP), allocatable :: buf(:)
        integer :: c, lo
        allocate(buf(max(1, io%Dn%w_nbuf)))
        call decomp_redistribute(io%Dn, owned, buf, 0.0_WP)
        do c = io%Dn%w_first_chunk, io%Dn%w_last_chunk
            lo = (c - io%Dn%w_first_chunk)*io%Dn%C + 1
            call zarr_write_chunk(store, arr, [c], buf(lo:lo + io%Dn%C - 1))
        end do
    end subroutine put_static

    integer function find_field(io, name) result(k)
        type(t_io_means), intent(in) :: io
        character(len=*), intent(in) :: name
        integer :: i
        k = -1
        do i = 1, io%nf
            if (trim(io%f(i)%name) == trim(name)) then; k = i; return; end if
        end do
    end function find_field

    subroutine read_env_int(name, val)
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: val
        character(len=64) :: buf
        integer :: ios, tmp
        call get_environment_variable(name, buf, status=ios)
        if (ios == 0 .and. len_trim(buf) > 0) then
            read(buf, *, iostat=ios) tmp
            if (ios == 0) val = tmp
        end if
    end subroutine read_env_int

end module mod_io_means
