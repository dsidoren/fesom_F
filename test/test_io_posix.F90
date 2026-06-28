program test_io_posix
    ! Task 3.1 GATE for mod_io_posix: the POSIX filesystem shims the restart writer needs
    ! (rename / unlink / rmdir / fsync-dir / recursive rmtree) work WITHOUT forking — so no
    ! segfault after MPI_Init (the M9 execute_command_line lesson) and a clean exit 0.
    !
    ! Exercise: build a nested multi-file tree (scratch/top/a/b/c.bin + a sibling file),
    ! fsync the dir, atomically rename the top dir, then posix_rmtree the renamed tree.
    ! Assert every shim returns success, the old path is gone after rename, the file is
    ! gone after rmtree, and a second rmtree of the now-missing dir reports an error (proof
    ! the tree is truly gone — inquire(file=<dir>) is unreliable for directories in Fortran,
    ! so we verify dir-removal by the file disappearing + the ENOENT re-rmtree). Also a
    ! single-file unlink + empty-dir rmdir round-trip. Reaching MPI_Finalize with nfail==0
    ! IS the no-fork / exit-0 proof.
    use mpi
    use mod_io_zarr,  only: zarr_mkdir   ! shared mkdir -p (now public) — builds the tree
    use mod_io_posix
    implicit none

    character(len=*), parameter :: base = 'test_io_posix_scratch'
    character(len=*), parameter :: top  = base//'/top'
    character(len=*), parameter :: ren  = base//'/top_renamed'
    integer :: ierr, nfail, r
    logical :: ex

    nfail = 0
    call MPI_Init(ierr)

    ! Clean any leftovers from a previous run (ignore status — may not exist).
    r = posix_rmtree(base)

    ! 1) Build a nested multi-file tree: base/top/a/b/c.bin + base/top/a/d.bin (bytes).
    call zarr_mkdir(top//'/a/b')
    call write_bytes(top//'/a/b/c.bin')
    call write_bytes(top//'/a/d.bin')
    inquire(file=top//'/a/b/c.bin', exist=ex); call check(ex,      'tree built: c.bin exists')
    inquire(file=top//'/a/d.bin',   exist=ex); call check(ex,      'tree built: d.bin exists')

    ! 2) fsync the top directory (durability before rename) — must succeed.
    r = posix_fsync_dir(top); call check(r == 0, 'posix_fsync_dir(top) == 0')

    ! 3) atomic rename of the top dir (same filesystem) — must succeed; contents follow.
    r = posix_rename(top, ren); call check(r == 0, 'posix_rename(top -> top_renamed) == 0')
    inquire(file=top//'/a/b/c.bin', exist=ex); call check(.not. ex, 'old path gone after rename')
    inquire(file=ren//'/a/b/c.bin', exist=ex); call check(ex,       'new path present after rename')

    ! 4) single-file unlink + empty-dir rmdir round-trip.
    call zarr_mkdir(base//'/solo')
    call write_bytes(base//'/solo/f.bin')
    r = posix_unlink(base//'/solo/f.bin'); call check(r == 0, 'posix_unlink(f.bin) == 0')
    inquire(file=base//'/solo/f.bin', exist=ex); call check(.not. ex, 'file gone after unlink')
    r = posix_rmdir(base//'/solo');        call check(r == 0, 'posix_rmdir(solo) == 0')

    ! 5) recursive delete of the renamed multi-file tree — must succeed and leave nothing.
    r = posix_rmtree(ren); call check(r == 0, 'posix_rmtree(renamed) == 0')
    inquire(file=ren//'/a/b/c.bin', exist=ex); call check(.not. ex, 'c.bin gone after rmtree')
    inquire(file=ren//'/a/d.bin',   exist=ex); call check(.not. ex, 'd.bin gone after rmtree')
    ! Robust dir-gone check: re-rmtree the now-missing dir must report an error (ENOENT).
    r = posix_rmtree(ren); call check(r /= 0, 're-rmtree of missing dir errors (tree truly gone)')

    ! 6) clean the scratch base.
    r = posix_rmtree(base)

    if (nfail == 0) then
        write(*,'(a)') 'test_io_posix: ALL PASS (rename/unlink/rmdir/fsync/rmtree, no fork, exit 0)'
    else
        write(*,'(a,i0,a)') 'test_io_posix: ', nfail, ' FAILURE(S)'
    end if
    call MPI_Finalize(ierr)
    if (nfail /= 0) error stop 1

contains

    subroutine write_bytes(path)
        character(len=*), intent(in) :: path
        integer :: u, ios
        open(newunit=u, file=path, form='unformatted', access='stream', &
             status='replace', action='write', iostat=ios)
        call check(ios == 0, 'open(write) '//path)
        if (ios == 0) then
            write(u) 'FESOM3-restart-shim-bytes'//char(10)
            close(u)
        end if
    end subroutine write_bytes

    subroutine check(cond, name)
        logical,          intent(in) :: cond
        character(len=*), intent(in) :: name
        if (.not. cond) then
            nfail = nfail + 1
            write(*,'(a)') '  FAIL: '//name
        end if
    end subroutine check

end program test_io_posix
