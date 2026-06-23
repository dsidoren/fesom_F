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
    use mod_step_oce,  only: step_oce
    implicit none
    private
    public :: t_model, model_init_analytic, model_step, model_finalize, allocate_state
    ! M2.9b: model_ocean_step is the model-level entry point for the assembled ocean step
    ! (mod_step_oce::step_oce). The forcing/diffusivity inputs are still passed in (D7) —
    ! they are sourced by the model itself once forcing (M2.10) lands; the full lifecycle
    ! driver (M2.11) calls this each timestep. The M2.9b byte-gate drives step_oce directly
    ! (src/drivers/fesom_stepdump.F90) so it can prescribe the inputs for the gate.
    public :: model_ocean_step

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

        ! ice prognostic (M3: data(1:3) = a_ice/m_ice/m_snow; sigma is EVP stress on elems).
        ! Convenience skeleton wiring only — the real ice allocation is ice_allocate (M3a).
        model%ice%num_itracers = 3
        allocate(model%ice%data(3))
        do i = 1, 3
            allocate(model%ice%data(i)%values(nn))
            model%ice%data(i)%values = 0.0_WP
            model%ice%data(i)%ID = i
        end do
        allocate(model%ice%uice(nn), model%ice%vice(nn))
        allocate(model%ice%work%sigma11(ne), model%ice%work%sigma12(ne), model%ice%work%sigma22(ne))
        model%ice%uice = 0.0_WP; model%ice%vice = 0.0_WP
        model%ice%work%sigma11 = 0.0_WP; model%ice%work%sigma12 = 0.0_WP; model%ice%work%sigma22 = 0.0_WP
    end subroutine allocate_state

    subroutine model_step(model)
        ! Lifecycle step counter. The assembled ocean dynamics are run via
        ! model_ocean_step (M2.9b), which the full driver (M2.11) wires in once forcing
        ! (M2.10) sources its inputs; the analytic M0 driver leaves this as the counter.
        type(t_model), intent(inout) :: model
        model%nsteps_done = model%nsteps_done + 1
    end subroutine model_step

    subroutine model_ocean_step(model, dt, lfirst, Ki, heat_flux, water_flux, &
                                virtual_salt, relax_salt, real_salt_flux, is_nonlinfs, &
                                stress_surf)
        ! Model-level wrapper for the assembled ocean timestep (mod_step_oce::step_oce):
        ! one faithful FESOM2 oce_timestep_ale on the model's prognostic state. The
        ! diffusivity/forcing inputs are explicit (D7) until the core sources them (M4 / M2.10).
        type(t_model), intent(inout) :: model
        real(kind=WP), intent(in) :: dt
        logical,       intent(in) :: lfirst
        real(kind=WP), intent(in) :: Ki(:,:)
        real(kind=WP), intent(in) :: heat_flux(:), water_flux(:), virtual_salt(:)
        real(kind=WP), intent(in) :: relax_salt(:), real_salt_flux(:), is_nonlinfs
        real(kind=WP), intent(in) :: stress_surf(:,:)
        call step_oce(model%nsteps_done+1, dt, lfirst, model%dyn, model%tracers, model%mesh, &
                      Ki, heat_flux, water_flux, virtual_salt, relax_salt, &
                      real_salt_flux, is_nonlinfs, stress_surf)
        model%nsteps_done = model%nsteps_done + 1
    end subroutine model_ocean_step

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
        if (allocated(model%ice%data))    deallocate(model%ice%data)
        call par_ex(comm, mype)
    end subroutine model_finalize

end module mod_model
