program fesom_restartcrash
    ! Restart Stage 3 (Task 3.3) crash-safety GATE: prove the atomic finalize + restart.latest pointer
    ! + keep-N prune behave correctly WITHOUT a real crash, by injecting the artifacts a crash leaves
    ! behind. The restart analog of fesom_restartsmoke; reuses the same pi-mesh, no-ocean scaffold.
    !
    !   1. Write checkpoint C1 (tag fesom.2000.001.03600) atomically. Assert restart.latest NAMES C1
    !      and restart_resolve_latest returns C1.
    !   2. Inject, by hand, the two things a crashed/interrupted writer can leave in the restart dir:
    !        - a stray  fesom.2000.001.05000.tmp/  (an unfinished staging dir from a crashed write), and
    !        - a finalized-but-UNPOINTED  fesom.2000.001.09000/  (a later checkpoint nothing points at).
    !      Assert restart_resolve_latest STILL returns C1 — it only ever follows restart.latest and never
    !      scans the dir, so both injected artifacts are ignored.
    !   3. Write checkpoint C2 (tag fesom.2000.001.07200) atomically with restart_keep=1. Assert the
    !      pointer FLIPPED to C2 (resolve returns C2), that the keep-N=1 prune removed C1 (its
    !      checkpoint.json is gone) while C2 remains, that the unpointed LATER folder survives (the prune
    !      protects the active pointer target even though 09000 sorts after C2 — so C2 itself is never
    !      pruned), and that the stray .tmp is left untouched (documented: never resolved, never pruned).
    !
    ! Prints PASS/<label> per assertion plus the pointer contents and resolved paths (before/after), then
    ! RESTARTCRASH OK iff every assertion passed (else error stop 1). No FESOM2 oracle needed.
    !
    !   FESOM3_MESH_DIR           mesh dir (default: pi)
    !   FESOM3_RESTARTCRASH_DIR   restart dir (default: ./restartcrash)
    use mpi
    use, intrinsic :: iso_fortran_env, only: real64
    use mod_precision,    only: WP
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex, set_partition
    use mod_mesh_read,    only: read_mesh
    use mod_mesh_areas,   only: compute_geometry
    use mod_part_bounds,  only: is_multirank, local_dims
    use mod_io_zarr,      only: zarr_mkdir
    use mod_io_decomp,    only: DECOMP_NODE
    use mod_io_restart
    implicit none
    character(len=512) :: mesh_dir, out_dir
    type(t_partit) :: partit
    type(t_mesh)   :: mesh
    type(t_restart) :: R
    real(WP), allocatable, target :: eta(:)
    real(WP), pointer :: peta(:)
    character(len=:), allocatable :: got
    character(len=512) :: ptr
    logical :: mr, ok, ex
    integer :: nsw, i, gid, ierr, mype, nfail
    integer :: nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF
    character(len=*), parameter :: C1 = 'fesom.2000.001.03600'
    character(len=*), parameter :: C2 = 'fesom.2000.001.07200'
    character(len=*), parameter :: STRAY_TMP = 'fesom.2000.001.05000.tmp'   ! crashed-write staging dir
    character(len=*), parameter :: LATER_UNP = 'fesom.2000.001.09000'       ! finalized but unpointed (later)

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call get_environment_variable('FESOM3_RESTARTCRASH_DIR', out_dir)
    if (len_trim(out_dir) == 0) out_dir = 'restartcrash'

    call par_init(partit)
    call set_partition(partit, trim(mesh_dir))
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)

    mype = partit%mype
    mr   = is_multirank(partit)
    call local_dims(mesh, partit, nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF)

    ! keep-N = 1: a 2nd checkpoint must prune the 1st down to the newest one.
    call restart_init(R, mesh, partit, chunk_horiz=1000, restart_keep=1)

    ! one synthetic node field (value(g)==g), associated as a live pointer (same as fesom_restartsmoke).
    allocate(eta(max(1,nNodO)))
    do i = 1, nNodO
        gid = i
        if (mr) gid = partit%myList_nod2D(i)
        eta(i) = real(gid, WP)
    end do
    peta => eta
    call restart_register_field(R, 'eta_n', 'm', DECOMP_NODE, p2d=peta)

    nfail = 0

    ! ---- 1) write C1 atomically; assert the pointer names it and resolve follows it ----------
    call restart_write(R, trim(out_dir), 2000, 1, 3600.0_real64, globalstep=1)
    if (mype == 0) then
        call read_pointer(trim(out_dir), ptr)
        write(*,'(a)') '  pointer after C1     = "'//trim(ptr)//'"'
        call check(trim(ptr) == C1, 'restart.latest names C1 after write')
        call restart_resolve_latest(trim(out_dir), got, ok)
        write(*,'(a,l1)') '  resolve after C1     = "'//trim(got)//'" ok=', ok
        call check(ok .and. trim(got) == trim(out_dir)//'/'//C1, 'restart_resolve_latest returns C1')
    end if

    ! ---- 2) inject a stray .tmp + an unpointed later finalized dir; resolve must STILL give C1 ----
    if (mype == 0) then
        call zarr_mkdir(trim(out_dir)//'/'//STRAY_TMP)                 ! crashed-write staging dir...
        call touch_file(trim(out_dir)//'/'//STRAY_TMP//'/partial.bin') ! ...with a partial artifact
        call zarr_mkdir(trim(out_dir)//'/'//LATER_UNP)                 ! finalized but UNPOINTED (later tag)
        call touch_file(trim(out_dir)//'/'//LATER_UNP//'/checkpoint.json')
        call restart_resolve_latest(trim(out_dir), got, ok)
        write(*,'(a,l1)') '  resolve after inject = "'//trim(got)//'" ok=', ok
        call check(ok .and. trim(got) == trim(out_dir)//'/'//C1, &
                   'resolve STILL returns C1 (ignores stray .tmp + unpointed later folder)')
    end if
    if (mr) call MPI_Barrier(partit%MPI_COMM_FESOM, ierr)

    ! ---- 3) write C2 atomically (keep-N=1): pointer flips, C1 pruned, C2 + protected later dir kept --
    call restart_write(R, trim(out_dir), 2000, 1, 7200.0_real64, globalstep=2)
    if (mype == 0) then
        call read_pointer(trim(out_dir), ptr)
        write(*,'(a)') '  pointer after C2     = "'//trim(ptr)//'"'
        call check(trim(ptr) == C2, 'restart.latest names C2 after 2nd write (atomic flip)')
        call restart_resolve_latest(trim(out_dir), got, ok)
        write(*,'(a,l1)') '  resolve after C2     = "'//trim(got)//'" ok=', ok
        call check(ok .and. trim(got) == trim(out_dir)//'/'//C2, 'restart_resolve_latest returns C2')
        inquire(file=trim(out_dir)//'/'//C1//'/checkpoint.json', exist=ex)
        call check(.not. ex, 'keep-N=1 pruned C1 (its checkpoint.json is gone)')
        inquire(file=trim(out_dir)//'/'//C2//'/checkpoint.json', exist=ex)
        call check(ex, 'C2 remains after prune (the active pointer target is never pruned)')
        inquire(file=trim(out_dir)//'/'//LATER_UNP//'/checkpoint.json', exist=ex)
        call check(ex, 'unpointed LATER folder survives (protect-target keeps C2 despite 09000>07200)')
        inquire(file=trim(out_dir)//'/'//STRAY_TMP//'/partial.bin', exist=ex)
        call check(ex, 'stray .tmp left untouched (never resolved, never pruned — documented)')
    end if

    ! agree across ranks, then a clean exit 0 on success / nonzero on any failure.
    call MPI_Bcast(nfail, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierr)
    if (mype == 0) then
        write(*,'(a,i0,a)') 'fesom_restartcrash: nod2D=', mesh%nod2D, ' done'
        if (nfail == 0) then
            write(*,'(a)') 'RESTARTCRASH OK'
        else
            write(*,'(a,i0,a)') 'RESTARTCRASH FAIL (', nfail, ' assertion(s) failed)'
        end if
    end if
    if (nfail /= 0) error stop 1
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)

contains

    subroutine check(cond, label)
        logical,          intent(in) :: cond
        character(len=*), intent(in) :: label
        if (cond) then
            write(*,'(a)') '  PASS  '//label
        else
            write(*,'(a)') '  FAIL  '//label
            nfail = nfail + 1
        end if
    end subroutine check

    ! Read the first line of <dir>/restart.latest into `line` (blank if missing/unreadable).
    subroutine read_pointer(dir, line)
        character(len=*), intent(in)  :: dir
        character(len=*), intent(out) :: line
        integer :: u, ios
        line = ''
        open(newunit=u, file=trim(dir)//'/restart.latest', status='old', action='read', &
             form='formatted', iostat=ios)
        if (ios /= 0) return
        read(u, '(a)', iostat=ios) line
        close(u)
        line = adjustl(line)
    end subroutine read_pointer

    ! Create an (empty-ish) file, mkdir-ing parents — used to fake crash artifacts.
    subroutine touch_file(path)
        character(len=*), intent(in) :: path
        integer :: u, ios
        open(newunit=u, file=trim(path), status='replace', action='write', &
             form='unformatted', access='stream', iostat=ios)
        if (ios /= 0) return
        write(u) 'x'//char(10)
        close(u)
    end subroutine touch_file

end program fesom_restartcrash
