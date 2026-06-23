program fesom_analytic
    ! M0.8 end-to-end driver (no physics): MPI -> partition -> analytic mesh ->
    ! allocate state -> run N empty steps -> finalize. Proves the full wiring of
    ! the foundation. Exits 0 on success, 1 on any internal check failure.
    !
    ! 1-rank is the real path; on multiple ranks each PE independently builds the
    ! full analytic mesh (no decomposition yet — multi-rank analytic partitioning
    ! is deferred), which still runs clean for a no-physics step.
    use mod_precision, only: WP
    use mod_model, only: t_model, model_init_analytic, model_step, model_finalize
    implicit none

    integer, parameter :: NX = 17, NY = 13, NL = 6, NSTEPS = 10
    type(t_model) :: model
    integer :: s, nfail

    call model_init_analytic(model, NX, NY, NL, &
                             Lx=170000.0_WP, Ly=130000.0_WP, max_depth=200.0_WP)

    nfail = 0
    if (model%mesh%nod2D  /= NX*NY)          nfail = nfail + 1
    if (model%mesh%elem2D /= 2*(NX-1)*(NY-1)) nfail = nfail + 1
    if (model%mesh%nl     /= NL)             nfail = nfail + 1
    if (.not. allocated(model%dyn%uv))       nfail = nfail + 1
    if (.not. allocated(model%tracers%data)) nfail = nfail + 1
    if (.not. allocated(model%ice%data))     nfail = nfail + 1
    if (model%mesh%edge2D <= 0)              nfail = nfail + 1
    if (.not. allocated(model%mesh%gradient_sca)) nfail = nfail + 1

    do s = 1, NSTEPS
        call model_step(model)
    end do
    if (model%nsteps_done /= NSTEPS) nfail = nfail + 1

    if (model%partit%mype == 0) then
        if (nfail == 0) then
            write(*,'(a,i0,a,i0,a,i0,a,i0,a)') 'fesom_analytic: OK  (nod2D=', &
                model%mesh%nod2D, ' elem2D=', model%mesh%elem2D, ' edge2D=', &
                model%mesh%edge2D, ', ', model%nsteps_done, ' steps)'
        else
            write(*,'(a,i0,a)') 'fesom_analytic: ', nfail, ' CHECK FAILURE(S)'
        end if
    end if

    call model_finalize(model)
    if (nfail /= 0) error stop 1
end program fesom_analytic
