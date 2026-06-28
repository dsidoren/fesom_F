module mod_io_restart
    ! FESOM3 restart / checkpoint WRITE path (restart Stage 3, Task 3.2).
    !
    ! A checkpoint is an immutable folder  <RestartOutPath>/fesom.<YYYY>.<DDD>.<SSSSS>/  holding ONE
    ! single-variable, single-entity SNAPSHOT Zarr store per registered prognostic field, plus a tiny
    ! checkpoint.json (rank-0 provenance). Each <name>.zarr is an M9 field store in every respect
    ! (same lon/lat + _ARRAY_DIMENSIONS + UGRID embedding via mod_io_coords) EXCEPT it has NO time
    ! dimension and NO mean accumulation: the data var is (entity) for a 2-D field or (nlev, entity)
    ! for a 3-D field. So ushow / xarray render a checkpoint store verbatim, and history is just
    ! separate per-checkpoint folders (immutable, atomic, prunable) instead of records in a file.
    !
    ! Reuses the proven M9 stack verbatim — NOT mod_io_means (whose growing time-dim / per-period /
    ! cadence / accumulation are the wrong tool for an immutable snapshot):
    !   - mod_io_decomp : canonical-order MPI_Alltoallv redistribution -> writer subset (no rank-0
    !                     gather), so the store is partition-INDEPENDENT (np=2 write reads at np=8).
    !   - mod_io_coords : the shared lon/lat + _ARRAY_DIMENSIONS + CF/UGRID-attr embedding (Task 2.1).
    !   - mod_io_zarr   : the Zarr v2 chunk writer; store-create ordering = rank 0 creates store +
    !                     defines all arrays -> barrier -> writers write their chunks.
    !
    ! REGISTRATION, not binary dumps: a field descriptor holds a Fortran POINTER to the LIVE model
    ! array (eta_n, uv components, tracers, ice arrays, sigma), exactly like M9's register_output_var.
    ! The reader (Task 4.1) reads back INTO the same pointers. Pointer dummies (not copies) so a strided
    ! model section such as dyn%uv(1,:,:) associates without a copy and stays live for read-back.
    !
    ! Caller protocol (all ranks, same order):
    !   restart_init(R, mesh, partit, ...)                          once (decomps + cached coords)
    !   restart_register_field(R, name, units, entity, p2d=.. )     per field, once (associates a pointer)
    !   ... at a checkpoint instant ...
    !   restart_write(R, dir, year, day, time_sec, globalstep)      writes the folder + per-field stores
    !   restart_finalize(R)                                         end-of-run barrier
    !
    ! SCOPE (Task 3.2): the general per-field writer + the folder + checkpoint.json. Atomic tmp/rename
    ! + restart.latest + keep-N prune is Task 3.3; registering the full oce+ice field set is Task 3.4;
    ! the READ path is Stage 4. This file deliberately implements only the WRITE mechanism + one field.
    use mpi
    use, intrinsic :: iso_fortran_env, only: real64
    use mod_precision,   only: WP
    use mod_mesh,        only: t_mesh
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: is_multirank, local_dims
    use mod_io_zarr
    use mod_io_decomp
    use mod_io_coords    ! shared lon/lat + _ARRAY_DIMENSIONS + UGRID-attr embedding (Task 2.1)
    implicit none
    private

    ! Re-export the entity selectors so a caller needs only `use mod_io_restart`.
    public :: DECOMP_NODE, DECOMP_ELEM
    public :: t_restart, t_restart_field, RESTART_MAXF
    public :: restart_init, restart_register_field, restart_write_field, restart_write, restart_finalize

    integer, parameter :: RESTART_MAXF = 64
    integer, parameter :: RESTART_FORMAT_VERSION = 1

    ! One registered prognostic field: its identity + a live POINTER to the model array it serializes.
    ! ndim=2 uses p2d (entity,); ndim=3 uses p3d (nlev, entity). on_full_levels picks nl levels ('nz')
    ! vs nl-1 layers ('nz1'), the M9 distinction — a wrong level count breaks the shape/gate.
    type :: t_restart_field
        character(len=64) :: name   = ''
        character(len=32) :: units  = ''
        integer           :: entity = DECOMP_NODE      ! DECOMP_NODE | DECOMP_ELEM
        integer           :: ndim   = 2                ! 2 or 3
        logical           :: on_full_levels = .false.  ! 3-D: nl levels (true) vs nl-1 layers (false)
        character(len=8)  :: dtype  = '<f8'            ! restart default = full precision
        integer           :: nlev   = 1                ! vertical size (3-D): nl or nl-1
        character(len=8)  :: hdim   = 'nod2'           ! 'nod2' (node) / 'elem' (element)
        character(len=8)  :: vdim   = ''               ! 'nz' (levels) / 'nz1' (layers)
        real(WP), pointer :: p2d(:)   => null()        ! live array (2-D field)
        real(WP), pointer :: p3d(:,:) => null()        ! live array (3-D field)
    end type t_restart_field

    type :: t_restart
        integer               :: nf = 0
        type(t_restart_field) :: f(RESTART_MAXF)
        type(t_io_decomp)     :: Dn                    ! node decomp (node fields)
        type(t_io_decomp)     :: De                    ! element decomp (element fields)
        integer               :: nNodO = 0, nElemO = 0
        integer               :: nl = 0
        integer               :: mype = 0, comm = MPI_COMM_SELF, npes = 1
        logical               :: mr = .false.
        ! global writer knobs (mirror mod_io_means)
        integer               :: chunk_vert = 0        ! vertical chunk (0 => full depth single chunk)
        character(len=16)     :: compressor = 'none'   ! data-array codec: 'none' | 'lz4'
        ! cached owned coords for the embed (node + element-centroid), like means_init
        real(WP), allocatable :: lon_n(:), lat_n(:), rlon_n(:), rlat_n(:)
        integer,  allocatable :: nlev_n(:)
        real(WP), allocatable :: lon_e(:), lat_e(:), rlon_e(:), rlat_e(:)
        integer,  allocatable :: nlev_e(:)
        real(WP), allocatable :: depth_nz(:), depth_nz1(:)  ! -zbar / -Z (CF positive-down depths)
    end type t_restart

contains

    ! ----------------------------------------------------------------- setup / registration

    subroutine restart_init(R, mesh, partit, chunk_horiz, n_writers, chunk_vert, compressor)
        type(t_restart),  intent(out)          :: R
        type(t_mesh),     intent(in)           :: mesh
        type(t_partit),   intent(in), optional :: partit
        integer,          intent(in), optional :: chunk_horiz, n_writers, chunk_vert
        character(len=*), intent(in), optional :: compressor
        integer :: C, nw, i, nNodL, nEdgeO, nEdgeL, nElemL, nElemF
        character(len=16) :: cbuf
        R%nf = 0
        R%mr   = is_multirank(partit)
        R%mype = 0; R%comm = MPI_COMM_SELF; R%npes = 1
        if (present(partit)) then
            R%mype = partit%mype; R%comm = partit%MPI_COMM_FESOM; R%npes = partit%npes
        end if
        ! writer knobs: arg overrides default; FESOM3_* env overrides arg (mirror means_init).
        C  = 500000; if (present(chunk_horiz)) C  = chunk_horiz
        nw = 0;      if (present(n_writers))   nw = n_writers
        R%chunk_vert = 0; if (present(chunk_vert)) R%chunk_vert = chunk_vert
        R%compressor = 'none'; if (present(compressor)) R%compressor = compressor
        call read_env_int('FESOM3_CHUNK_HORIZ', C)
        call read_env_int('FESOM3_N_WRITERS',  nw)
        call read_env_int('FESOM3_CHUNK_VERT', R%chunk_vert)
        cbuf = ''; call get_environment_variable('FESOM3_COMPRESSOR', cbuf)
        if (len_trim(cbuf) > 0) R%compressor = cbuf
        call local_dims(mesh, partit, R%nNodO, nNodL, nEdgeO, nEdgeL, R%nElemO, nElemL, nElemF)
        call decomp_init_entity(R%Dn, C, nw, DECOMP_NODE, mesh, partit)
        call decomp_init_entity(R%De, C, nw, DECOMP_ELEM, mesh, partit)
        R%nl = mesh%nl
        ! cached owned coords (geographic deg lon/lat for the embed + rotated rad + nlevels mask).
        allocate(R%lon_n(max(1,R%nNodO)),  R%lat_n(max(1,R%nNodO)), &
                 R%rlon_n(max(1,R%nNodO)), R%rlat_n(max(1,R%nNodO)), R%nlev_n(max(1,R%nNodO)))
        call io_coords_compute(DECOMP_NODE, mesh, R%nNodO, R%lon_n, R%lat_n, R%rlon_n, R%rlat_n, R%nlev_n)
        allocate(R%lon_e(max(1,R%nElemO)),  R%lat_e(max(1,R%nElemO)), &
                 R%rlon_e(max(1,R%nElemO)), R%rlat_e(max(1,R%nElemO)), R%nlev_e(max(1,R%nElemO)))
        call io_coords_compute(DECOMP_ELEM, mesh, R%nElemO, R%lon_e, R%lat_e, R%rlon_e, R%rlat_e, R%nlev_e)
        ! vertical coords (CF positive-down depths): nz = -zbar(1:nl), nz1 = -Z(1:nl-1)
        allocate(R%depth_nz(R%nl), R%depth_nz1(R%nl-1))
        do i = 1, R%nl;   R%depth_nz(i)  = real(-mesh%zbar(i), WP); end do
        do i = 1, R%nl-1; R%depth_nz1(i) = real(-mesh%Z(i),    WP); end do
    end subroutine restart_init

    ! Register one prognostic field by associating a LIVE pointer to its model array. Pass p2d for a
    ! 2-D field, p3d for a 3-D field (on_full_levels selects nl vs nl-1). entity = DECOMP_NODE/ELEM.
    ! The pointer is stored, not copied, so the writer reads the live values and the reader (Task 4.1)
    ! can write back into the same array. Caller passes a Fortran pointer (e.g. peta => dyn%eta_n).
    subroutine restart_register_field(R, name, units, entity, p2d, p3d, on_full_levels, precision)
        type(t_restart),  intent(inout)        :: R
        character(len=*), intent(in)           :: name, units
        integer,          intent(in)           :: entity
        real(WP), pointer, intent(in), optional :: p2d(:)
        real(WP), pointer, intent(in), optional :: p3d(:,:)
        logical,          intent(in), optional :: on_full_levels
        character(len=*), intent(in), optional :: precision
        integer :: k
        call zarr_check(R%nf < RESTART_MAXF, 'restart: too many fields ('//trim(name)//')')
        R%nf = R%nf + 1; k = R%nf
        R%f(k)%name   = name
        R%f(k)%units  = units
        R%f(k)%entity = entity
        R%f(k)%hdim   = 'nod2'; if (entity == DECOMP_ELEM) R%f(k)%hdim = 'elem'
        R%f(k)%dtype  = restart_dtype(precision)
        if (present(p3d)) then
            R%f(k)%ndim = 3
            R%f(k)%p3d => p3d
            R%f(k)%on_full_levels = .false.
            if (present(on_full_levels)) R%f(k)%on_full_levels = on_full_levels
            if (R%f(k)%on_full_levels) then
                R%f(k)%nlev = R%nl;   R%f(k)%vdim = 'nz'
            else
                R%f(k)%nlev = R%nl-1; R%f(k)%vdim = 'nz1'
            end if
        else
            R%f(k)%ndim = 2
            R%f(k)%nlev = 1
            if (present(p2d)) R%f(k)%p2d => p2d
        end if
    end subroutine restart_register_field

    ! ----------------------------------------------------------------- checkpoint write

    ! Write one full checkpoint: create the folder fesom.<YYYY>.<DDD>.<SSSSS>/ (zero-padded; SSSSS =
    ! int(time_sec) = sec-of-day / timenew), write each registered field's snapshot store into it, and
    ! write checkpoint.json (rank 0). time_sec is the new-clock sec-of-day; Stage 5 passes timenew.
    ! (Atomic tmp/rename + restart.latest is Task 3.3 — a plain mkdir + direct write is used here.)
    subroutine restart_write(R, checkpoint_dir, year, day, time_sec, globalstep)
        type(t_restart),  intent(in), target :: R
        character(len=*), intent(in)         :: checkpoint_dir
        integer,          intent(in)         :: year, day
        real(real64),     intent(in)         :: time_sec
        integer,          intent(in), optional :: globalstep
        integer :: k, ierr, gstep
        character(len=4) :: cy
        character(len=3) :: cd
        character(len=5) :: cs
        character(len=:), allocatable :: folder
        gstep = 0; if (present(globalstep)) gstep = globalstep
        write(cy, '(i4.4)') year
        write(cd, '(i3.3)') day
        write(cs, '(i5.5)') int(time_sec)
        folder = trim(checkpoint_dir)//'/fesom.'//cy//'.'//cd//'.'//cs
        if (R%mype == 0) call zarr_mkdir(folder)          ! recursive mkdir -p (also creates the parent)
        if (R%mr) call MPI_Barrier(R%comm, ierr)          ! folder must exist before any store write
        do k = 1, R%nf
            call restart_write_field(R, R%f(k), trim(folder)//'/'//trim(R%f(k)%name)//'.zarr')
        end do
        if (R%mype == 0) call write_checkpoint_json(R, trim(folder)//'/checkpoint.json', year, day, &
                                                    time_sec, gstep)
        if (R%mr) call MPI_Barrier(R%comm, ierr)
    end subroutine restart_write

    ! Write ONE field's snapshot store (general: node OR element, 2-D OR 3-D). Resolves the field's
    ! entity decomp + cached coords from R, creates a single-variable single-entity store (data var
    ! (entity) or (nlev, entity), <f8 by default, codec from the namelist), embeds lon/lat +
    ! _ARRAY_DIMENSIONS + UGRID attrs via the Stage 2 helper, then decomp_redistributes the live array
    ! to canonical writer buffers and writes the data chunks. Store-create ordering: rank 0 defines ->
    ! barrier -> every writer writes its chunks.
    subroutine restart_write_field(R, f, store_path)
        type(t_restart),       intent(in), target :: R
        type(t_restart_field), intent(in)         :: f
        character(len=*),      intent(in)         :: store_path
        type(t_io_decomp), pointer :: D
        real(WP),          pointer :: own_lon(:), own_lat(:)
        integer                    :: nO
        type(t_zarr_store) :: store
        type(t_zarr_array) :: a_data, a_vert, a_lon, a_lat
        type(t_zarr_attrs) :: gat, at
        real(WP), allocatable :: buf(:), buf3(:,:)
        integer :: c, lo, ierr, cv, nvc, vc, L0, cvn
        ! resolve the entity context (node vs element-centroid)
        if (f%entity == DECOMP_ELEM) then
            D => R%De; nO = R%nElemO; own_lon => R%lon_e; own_lat => R%lat_e
        else
            D => R%Dn; nO = R%nNodO; own_lon => R%lon_n; own_lat => R%lat_n
        end if
        cv = vchunk_eff(R%chunk_vert, max(1, f%nlev))
        ! ---- array handles (all ranks) ----
        if (f%ndim == 2) then
            call zarr_array_init(a_data, trim(f%name), [D%N], [D%C], trim(f%dtype), &
                                 has_fill=.false., codec=trim(R%compressor))
        else
            call zarr_array_init(a_data, trim(f%name), [f%nlev, D%N], [cv, D%C], trim(f%dtype), &
                                 has_fill=.false., codec=trim(R%compressor))
            call zarr_array_init(a_vert, trim(f%vdim), [f%nlev], [f%nlev], '<f8', has_fill=.false.)
        end if
        call io_coords_init_lonlat(a_lon, a_lat, D%N, D%C)
        ! the chunk writers (every writer) reference store%path -> set on ALL ranks
        store%path = trim(store_path)
        ! ---- rank 0 defines the store + arrays + attrs (store-create ordering) ----
        if (R%mype == 0) then
            call zarr_attrs_init(gat)
            call zattr_str(gat, 'Conventions', 'CF-1.8')
            call zattr_str(gat, 'description', 'FESOM3 restart checkpoint field (Zarr v2 snapshot)')
            call zarr_create_store(store, trim(store_path), gat)
            call zarr_attrs_init(at)
            if (f%ndim == 2) then
                call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: f%hdim])
            else
                call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: f%vdim, f%hdim])
            end if
            if (len_trim(f%units) > 0) call zattr_str(at, 'units', trim(f%units))
            call zattr_str(at, 'coordinates', 'lon lat')
            call zarr_define_array(store, a_data, at)
            if (f%ndim == 3) then
                call zarr_attrs_init(at)
                call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: f%vdim])
                call zattr_str(at, 'long_name', 'depth'); call zattr_str(at, 'units', 'meters')
                call zattr_str(at, 'positive', 'down')
                call zarr_define_array(store, a_vert, at)
            end if
            call io_coords_define_lonlat(store, a_lon, a_lat, f%hdim)
        end if
        if (R%mr) call MPI_Barrier(R%comm, ierr)
        ! ---- every writer redistributes the live field + writes its data chunks ----
        if (f%ndim == 2) then
            allocate(buf(max(1, D%w_nbuf)))
            call decomp_redistribute(D, f%p2d(1:nO), buf, 0.0_WP)            ! collective
            do c = D%w_first_chunk, D%w_last_chunk
                lo = (c - D%w_first_chunk)*D%C + 1
                call zarr_write_chunk(store, a_data, [c], buf(lo:lo + D%C - 1))
            end do
        else
            allocate(buf3(f%nlev, max(1, D%w_nbuf)))
            call decomp_redistribute(D, f%p3d(1:f%nlev,1:nO), buf3, 0.0_WP)  ! collective
            nvc = (f%nlev + cv - 1)/cv
            do c = D%w_first_chunk, D%w_last_chunk
                lo = (c - D%w_first_chunk)*D%C + 1
                do vc = 0, nvc - 1
                    L0 = vc*cv; cvn = min(cv, f%nlev - L0)
                    call zarr_write_chunk(store, a_data, [vc, c], buf3(L0+1:L0+cvn, lo:lo + D%C - 1))
                end do
            end do
        end if
        ! embedded lon/lat coord data (shared helper redistributes + writes per-chunk on each writer)
        call io_coords_put(D, store, a_lon, own_lon)
        call io_coords_put(D, store, a_lat, own_lat)
        ! vertical coord is global/replicated -> rank 0 writes it whole
        if (R%mype == 0 .and. f%ndim == 3) then
            if (trim(f%vdim) == 'nz') then
                call zarr_write_whole(store, a_vert, R%depth_nz(1:R%nl))
            else
                call zarr_write_whole(store, a_vert, R%depth_nz1(1:R%nl-1))
            end if
        end if
        if (R%mr) call MPI_Barrier(R%comm, ierr)
    end subroutine restart_write_field

    subroutine restart_finalize(R)
        type(t_restart), intent(in) :: R
        integer :: ierr
        if (R%mr) call MPI_Barrier(R%comm, ierr)
    end subroutine restart_finalize

    ! ----------------------------------------------------------------- checkpoint.json (rank 0)

    ! Plain Fortran string writes (no JSON library): the small provenance manifest the reader uses for
    ! the time-vs-clock safety check + npes provenance. format_version/year/day/globalstep/npes_wrote
    ! are integers; time_sec a real (the sec-of-day passed in); fesom_git from $FESOM3_GIT else 'unknown'.
    subroutine write_checkpoint_json(R, path, year, day, time_sec, gstep)
        type(t_restart),  intent(in) :: R
        character(len=*), intent(in) :: path
        integer,          intent(in) :: year, day, gstep
        real(real64),     intent(in) :: time_sec
        character(len=64) :: git
        character(len=32) :: rbuf
        character(len=:), allocatable :: doc, rstr
        write(rbuf, '(es24.16)') time_sec
        rstr = trim(adjustl(rbuf))
        git = ''; call get_environment_variable('FESOM3_GIT', git)
        if (len_trim(git) == 0) git = 'unknown'
        doc = '{'//char(10)// &
              '  "format_version": '//itoa(RESTART_FORMAT_VERSION)//','//char(10)// &
              '  "year": '//itoa(year)//','//char(10)// &
              '  "day": '//itoa(day)//','//char(10)// &
              '  "time_sec": '//rstr//','//char(10)// &
              '  "globalstep": '//itoa(gstep)//','//char(10)// &
              '  "fesom_git": "'//trim(git)//'",'//char(10)// &
              '  "npes_wrote": '//itoa(R%npes)//char(10)// &
              '}'
        call write_text_file(path, doc)
    end subroutine write_checkpoint_json

    ! ----------------------------------------------------------------- small helpers

    ! map a precision arg ('double'|'<f8'|'8' => <f8; '<f4'|'4'|'single' => <f4) to a Zarr dtype.
    ! restart DEFAULT = <f8 (full precision; lossy restart is not a v1 option).
    function restart_dtype(precision) result(dt)
        character(len=*), intent(in), optional :: precision
        character(len=8) :: dt
        dt = '<f8'
        if (present(precision)) then
            if (trim(precision) == '<f4' .or. trim(precision) == '4' .or. trim(precision) == 'single') &
                dt = '<f4'
        end if
    end function restart_dtype

    ! effective vertical chunk: full depth when chunk_vert is 0 or >= nlev, else chunk_vert.
    integer function vchunk_eff(cv, nlev)
        integer, intent(in) :: cv, nlev
        if (cv <= 0 .or. cv >= nlev) then
            vchunk_eff = max(1, nlev)
        else
            vchunk_eff = cv
        end if
    end function vchunk_eff

    function itoa(i) result(s)
        integer, intent(in)           :: i
        character(len=:), allocatable :: s
        character(len=32) :: tmp
        write(tmp, '(i0)') i
        s = trim(tmp)
    end function itoa

    ! Write a text blob verbatim + a trailing newline via stream (any length is safe). Mirrors
    ! mod_io_zarr's private write_text_file (kept local so mod_io_restart owns its manifest I/O).
    subroutine write_text_file(path, text)
        character(len=*), intent(in) :: path, text
        integer :: u, ios
        open(newunit=u, file=trim(path), status='replace', action='write', &
             form='unformatted', access='stream', iostat=ios)
        call zarr_check(ios == 0, 'restart: open(write) '//trim(path))
        write(u) text//char(10)
        close(u)
    end subroutine write_text_file

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

end module mod_io_restart
