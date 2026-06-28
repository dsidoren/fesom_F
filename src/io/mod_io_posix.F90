module mod_io_posix
    ! Thin Fortran wrappers over POSIX filesystem ops for the restart writer (Task 3.1).
    !
    ! These bind libc symbols DIRECTLY via iso_c_binding — exactly like mod_io_zarr's
    ! c_mkdir (mod_io_zarr.F90:34-44) — so we NEVER fork. The M9 lesson: a forking
    ! execute_command_line('mv'/'rm') segfaults AFTER MPI_Init on Levante's OpenMPI/vader
    ! stack (the project's MPI-fork gotcha), which is why c_mkdir exists; mirror it.
    !
    ! No C source file is compiled into the build (mod_io_zarr binds mkdir / liblz4 the
    ! same way, and the repo ships zero .c files); a new .F90 is auto-included by the
    ! src/io/*.F90 CONFIGURE_DEPENDS GLOB in the top-level CMakeLists.
    !
    ! The restart writer (Stage 3) needs: atomic finalize + pointer-file swap (rename),
    ! durability before rename (fsync of the containing directory), and tmp-cleanup /
    ! keep-N prune of a multi-file Zarr store directory tree (posix_rmtree).
    !
    ! Each public wrapper takes a Fortran character(len=*) path, appends c_null_char, and
    ! returns the C int status (0 = success, -1 = error with errno set). x86-64 Linux,
    ! little-endian LP64 — matches the rest of the build.
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_ptr, c_funptr, c_funloc, &
                                           c_null_char
    implicit none
    private

    public :: posix_rename, posix_unlink, posix_rmdir, posix_fsync_dir, posix_rmtree

    ! POSIX open(2) flag and nftw(3) flags we use (Linux x86-64 glibc values).
    integer(c_int), parameter :: O_RDONLY  = 0_c_int
    integer(c_int), parameter :: FTW_PHYS  = 1_c_int   ! physical walk; do not follow symlinks
    integer(c_int), parameter :: FTW_DEPTH = 8_c_int   ! post-order: a dir's contents before the dir

    interface
        ! int rename(const char *oldpath, const char *newpath);  — atomic within one filesystem.
        function c_rename(oldp, newp) bind(C, name="rename") result(r)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: oldp, newp
            integer(c_int) :: r
        end function c_rename

        ! int unlink(const char *pathname);  — remove a file (or symlink).
        function c_unlink(path) bind(C, name="unlink") result(r)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: path
            integer(c_int) :: r
        end function c_unlink

        ! int rmdir(const char *pathname);  — remove an EMPTY directory.
        function c_rmdir(path) bind(C, name="rmdir") result(r)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: path
            integer(c_int) :: r
        end function c_rmdir

        ! int open(const char *pathname, int flags);  — C's open is variadic
        ! (int open(const char*, int, ...)); the 2-arg form is the O_RDONLY case we use
        ! (no mode without O_CREAT). The x86-64 SysV ABI passes a normal (non-variadic)
        ! call cleanly here since glibc never reads the absent mode vararg for O_RDONLY.
        function c_open(path, flags) bind(C, name="open") result(fd)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: path
            integer(c_int), value :: flags
            integer(c_int) :: fd
        end function c_open

        ! int fsync(int fd);  — flush a descriptor (here a directory fd) to stable storage.
        function c_fsync(fd) bind(C, name="fsync") result(r)
            import :: c_int
            integer(c_int), value :: fd
            integer(c_int) :: r
        end function c_fsync

        ! int close(int fd);
        function c_close(fd) bind(C, name="close") result(r)
            import :: c_int
            integer(c_int), value :: fd
            integer(c_int) :: r
        end function c_close

        ! int remove(const char *pathname);  — unlink() for files, rmdir() for empty dirs.
        ! Called from the nftw callback with the C path pointer it already holds, so the
        ! argument is a by-value c_ptr (no Fortran copy, no re-null-termination needed).
        function c_remove_ptr(path) bind(C, name="remove") result(r)
            import :: c_ptr, c_int
            type(c_ptr), value :: path
            integer(c_int)     :: r
        end function c_remove_ptr

        ! int nftw(const char *dirpath, int (*fn)(const char*, const struct stat*, int,
        !          struct FTW*), int nopenfd, int flags);
        ! With FTW_DEPTH it visits a directory's contents BEFORE the directory itself, so
        ! the callback can remove() every entry bottom-up. We never dereference the struct
        ! stat pointer, so the nftw-vs-nftw64 stat layout is irrelevant on LP64 x86-64.
        function c_nftw(dirpath, fn, nopenfd, flags) bind(C, name="nftw") result(r)
            import :: c_char, c_funptr, c_int
            character(kind=c_char), dimension(*), intent(in) :: dirpath
            type(c_funptr), value :: fn
            integer(c_int), value :: nopenfd, flags
            integer(c_int)        :: r
        end function c_nftw
    end interface

    ! Failure accumulator for the nftw callback. The C-driven callback can only see module
    ! state (not posix_rmtree's locals), so we count remove() failures here. SAFE because
    ! posix_rmtree is called serially from a single rank (the restart writer's rank-0
    ! tmp-cleanup / keep-N prune, and the gate), never concurrently; reset before each walk.
    integer, save :: rmtree_nerr = 0

contains

    ! --- atomic rename: tmp -> final store, and the restart.latest pointer-file swap ----
    integer function posix_rename(oldpath, newpath) result(r)
        character(len=*), intent(in) :: oldpath, newpath
        r = int(c_rename(trim(oldpath)//c_null_char, trim(newpath)//c_null_char))
    end function posix_rename

    ! --- remove a single file -----------------------------------------------------------
    integer function posix_unlink(path) result(r)
        character(len=*), intent(in) :: path
        r = int(c_unlink(trim(path)//c_null_char))
    end function posix_unlink

    ! --- remove an EMPTY directory ------------------------------------------------------
    integer function posix_rmdir(path) result(r)
        character(len=*), intent(in) :: path
        r = int(c_rmdir(trim(path)//c_null_char))
    end function posix_rmdir

    ! --- fsync a directory so a rename into/within it survives a crash ------------------
    ! Open the path read-only, fsync the fd, close it. (O_RDONLY works on a regular file
    ! too, so this doubles as a file fsync if ever needed.) Returns 0 only on full success.
    integer function posix_fsync_dir(path) result(r)
        character(len=*), intent(in) :: path
        integer(c_int) :: fd, rc
        fd = c_open(trim(path)//c_null_char, O_RDONLY)
        if (fd < 0_c_int) then
            r = -1
            return
        end if
        rc = c_fsync(fd)
        if (c_close(fd) /= 0_c_int) rc = -1_c_int
        r = int(rc)
    end function posix_fsync_dir

    ! --- recursively delete a whole tree (a multi-file Zarr store dir) ------------------
    ! nftw(FTW_DEPTH|FTW_PHYS) drives a post-order walk; the callback remove()s each entry
    ! bottom-up (chunk files + .zarray/.zattrs, then the now-empty subdirs, then the root).
    ! Returns 0 only if nftw itself succeeded AND every remove() succeeded; a missing
    ! dirpath surfaces as the negative nftw return (ENOENT), which is the right signal for
    ! tmp-cleanup of a path that may not exist.
    integer function posix_rmtree(path) result(r)
        character(len=*), intent(in) :: path
        integer(c_int) :: rc
        rmtree_nerr = 0
        rc = c_nftw(trim(path)//c_null_char, c_funloc(rmtree_cb), 64_c_int, &
                    ior(FTW_DEPTH, FTW_PHYS))
        if (rc /= 0_c_int) then
            r = int(rc)            ! e.g. dirpath missing (ENOENT) -> -1
        else if (rmtree_nerr /= 0) then
            r = -1                 ! tree partially survived (a remove() failed)
        else
            r = 0
        end if
    end function posix_rmtree

    ! nftw callback. Signature MUST match the C type
    !   int (*)(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf)
    ! We use only fpath (the C string pointer nftw passes us). Always returns 0 so the walk
    ! continues and deletes as much as possible; a failed remove() bumps rmtree_nerr, which
    ! posix_rmtree turns into a -1 status. sb/typeflag/ftwbuf are deliberately unused.
    function rmtree_cb(fpath, sb, typeflag, ftwbuf) bind(C) result(r)
        type(c_ptr),    value :: fpath
        type(c_ptr),    value :: sb
        integer(c_int), value :: typeflag
        type(c_ptr),    value :: ftwbuf
        integer(c_int)        :: r
        if (c_remove_ptr(fpath) /= 0_c_int) rmtree_nerr = rmtree_nerr + 1
        r = 0_c_int
    end function rmtree_cb

end module mod_io_posix
