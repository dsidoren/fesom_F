module mod_model
    ! Orchestration aggregate + lifecycle (decisions D6). t_model is a convenience
    ! container used ONLY at driver/orchestration level; kernels take the specific
    ! types they touch (dependency injection). Lifecycle: model_init / model_step /
    ! model_finalize. model_step is empty in M0 (physics enters at M1).
    use mod_precision, only: WP
    use mod_mesh,      only: t_mesh
    use mod_partit,    only: t_partit
    use mod_dyn,       only: t_dyn
    use mod_tracer,    only: t_tracer
    use mod_ice,       only: t_ice
    use mod_partitioning, only: par_init, par_ex
    use mod_mesh_analytic, only: generate_analytic_mesh
    implicit none
    private
    public :: t_model, model_init_analytic, model_step, model_finalize, allocate_state

    type t_model
        type(t_partit) :: partit
        type(t_mesh)   :: mesh
        type(t_dyn)    :: dyn
        type(t_tracer) :: tracers
        type(t_ice)    :: ice
        integer :: nsteps_done = 0
    end type t_model

contains

    subroutine model_init_analytic(model, nx, ny, nl, Lx, Ly, max_depth)
        ! MPI/partition -> analytic mesh -> allocate evolving state. No physics.
        type(t_model), intent(inout) :: model
        integer,       intent(in)    :: nx, ny, nl
        real(kind=WP), intent(in)    :: Lx, Ly, max_depth
        call par_init(model%partit)
        call generate_analytic_mesh(model%mesh, model%partit, nx, ny, nl, Lx, Ly, max_depth)
        call allocate_state(model)
    end subroutine model_init_analytic

    subroutine allocate_state(model)
        ! Size the prognostic evolving-type arrays to the mesh. (No physics; this
        ! just proves the data model wires to the mesh. Work/aux arrays grow per
        ! milestone as kernels need them.)
        type(t_model), intent(inout) :: model
        integer :: ne, nn, nl, i
        ne = model%mesh%elem2D; nn = model%mesh%nod2D; nl = model%mesh%nl

        ! dynamics prognostic
        allocate(model%dyn%uv(2, nl-1, ne), model%dyn%uv_rhs(2, nl-1, ne))
        allocate(model%dyn%uv_rhsAB(2, nl-1, ne, model%dyn%AB_order))
        allocate(model%dyn%uvnode(2, nl-1, nn))
        allocate(model%dyn%w(nl, nn), model%dyn%w_e(nl, nn), model%dyn%w_i(nl, nn))
        allocate(model%dyn%eta_n(nn), model%dyn%d_eta(nn))
        allocate(model%dyn%ssh_rhs(nn), model%dyn%ssh_rhs_old(nn))
        model%dyn%uv = 0.0_WP; model%dyn%uv_rhs = 0.0_WP; model%dyn%uv_rhsAB = 0.0_WP
        model%dyn%uvnode = 0.0_WP
        model%dyn%w = 0.0_WP; model%dyn%w_e = 0.0_WP; model%dyn%w_i = 0.0_WP
        model%dyn%eta_n = 0.0_WP; model%dyn%d_eta = 0.0_WP
        model%dyn%ssh_rhs = 0.0_WP; model%dyn%ssh_rhs_old = 0.0_WP

        ! tracers (T, S)
        model%tracers%num_tracers = 2
        allocate(model%tracers%data(model%tracers%num_tracers))
        do i = 1, model%tracers%num_tracers
            allocate(model%tracers%data(i)%values(nl-1, nn))
            model%tracers%data(i)%values = 0.0_WP
            model%tracers%data(i)%ID = i
        end do

        ! ice prognostic
        allocate(model%ice%a_ice(nn), model%ice%m_ice(nn), model%ice%m_snow(nn))
        allocate(model%ice%uice(nn), model%ice%vice(nn))
        allocate(model%ice%sigma11(ne), model%ice%sigma12(ne), model%ice%sigma22(ne))
        model%ice%a_ice = 0.0_WP; model%ice%m_ice = 0.0_WP; model%ice%m_snow = 0.0_WP
        model%ice%uice = 0.0_WP; model%ice%vice = 0.0_WP
        model%ice%sigma11 = 0.0_WP; model%ice%sigma12 = 0.0_WP; model%ice%sigma22 = 0.0_WP
    end subroutine allocate_state

    subroutine model_step(model)
        ! Empty in M0 (the faithful ocean/ice sequence enters at M1/M2).
        type(t_model), intent(inout) :: model
        model%nsteps_done = model%nsteps_done + 1
    end subroutine model_step

    subroutine model_finalize(model)
        ! Explicit deallocation (derived-type allocatables would also auto-release
        ! at scope exit) + MPI finalize.
        type(t_model), intent(inout) :: model
        integer :: comm, mype
        comm = model%partit%MPI_COMM_FESOM; mype = model%partit%mype
        if (allocated(model%dyn%uv))      deallocate(model%dyn%uv)
        if (allocated(model%dyn%uv_rhs))  deallocate(model%dyn%uv_rhs)
        if (allocated(model%dyn%uv_rhsAB))deallocate(model%dyn%uv_rhsAB)
        if (allocated(model%dyn%uvnode))  deallocate(model%dyn%uvnode)
        if (allocated(model%dyn%w))       deallocate(model%dyn%w)
        if (allocated(model%tracers%data))deallocate(model%tracers%data)
        if (allocated(model%ice%a_ice))   deallocate(model%ice%a_ice)
        call par_ex(comm, mype)
    end subroutine model_finalize

end module mod_model
