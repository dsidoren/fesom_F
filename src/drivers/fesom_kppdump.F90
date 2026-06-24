program fesom_kppdump
    ! M5a-1 gate driver: build the KPP rank-independent constants + lookup tables
    ! (oce_mixing_kpp_init) and the wscale sweep, dumping them in the EXACT format of the
    ! oracle's FESOM_KPP_DUMP_DIR instrumentation (oce_ale_mixing_kpp.F90:221-266) so they
    ! diff byte-for-byte (es24.16 = exact double round-trip).
    !
    ! Output dir from env FESOM3_KPP_OUT (default "."). Writes:
    !   <dir>/kpp_init_rank0.txt    — Vtc/cg/deltaz/deltau + the wmt/wst tables
    !   <dir>/kpp_wscale_rank0.txt  — wscale(zehat,us) over the same 201x101 sweep grid
    !
    ! Config = work_core namelist.oce DOUBLES: Ricr=0.3, concv=1.6.
    use mpi
    use mod_precision,  only: WP
    use oce_mixing_kpp, only: oce_mixing_kpp_init, wscale, &
                              Vtc, cg, deltaz, deltau, wmt, wst, nni, nnj
    implicit none

    integer :: ierr, rank
    character(len=512) :: outdir
    character(len=640) :: path
    integer :: u, ios, ii, jj
    integer, parameter :: NZ = 201, NU = 101
    real(kind=WP) :: zehat, us, wm, ws

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    if (rank == 0) then
        call get_environment_variable('FESOM3_KPP_OUT', outdir, status=ios)
        if (ios /= 0 .or. len_trim(outdir) == 0) outdir = '.'

        ! work_core: Ricr=0.3, concv=1.6
        call oce_mixing_kpp_init(0.3_WP, 1.6_WP)

        ! --- K1: init constants + wmt/wst lookup tables --------------------------
        write(path,'(a,"/kpp_init_rank0.txt")') trim(outdir)
        open(newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,*) '[kppdump] cannot open ', trim(path); call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        write(u,'("# Vtc ",es24.16)')    Vtc
        write(u,'("# cg ",es24.16)')     cg
        write(u,'("# deltaz ",es24.16)') deltaz
        write(u,'("# deltau ",es24.16)') deltau
        write(u,'("# i j wmt wst  (nni=",i0," nnj=",i0,")")') nni, nnj
        do ii=0, nni+1
            do jj=0, nnj+1
                write(u,'(i0," ",i0," ",es24.16," ",es24.16)') ii, jj, wmt(ii,jj), wst(ii,jj)
            end do
        end do
        close(u)

        ! --- K2: wscale sweep (identical grid to the oracle) ---------------------
        write(path,'(a,"/kpp_wscale_rank0.txt")') trim(outdir)
        open(newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,*) '[kppdump] cannot open ', trim(path); call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        write(u,'("# i j zehat ustar wm ws  (sweep ",i0,"x",i0,")")') NZ, NU
        do ii=0, NZ-1
            zehat = -1.0e-6_WP + real(ii,WP) * (2.0e-6_WP / real(NZ-1,WP))
            do jj=0, NU-1
                us = real(jj,WP) * (0.05_WP / real(NU-1,WP))
                call wscale(zehat, us, wm, ws)
                write(u,'(i0," ",i0," ",es24.16," ",es24.16," ",es24.16," ",es24.16)') &
                     ii, jj, zehat, us, wm, ws
            end do
        end do
        close(u)

        write(*,'(a)') '[kppdump] wrote kpp_init_rank0.txt + kpp_wscale_rank0.txt to '//trim(outdir)
    end if

    call MPI_Finalize(ierr)
end program fesom_kppdump
