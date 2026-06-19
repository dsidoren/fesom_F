module oce_ale_tracer
    ! Tracer step integration (the advection part), transcribed from FESOM2 v2.7.3
    ! oce_ale_tracer.F90:146-271 (SUBROUTINE adv_tracers_ale, advection portion).
    !
    ! adv_tracers_ale loops over tracers; per tracer it runs the advection body
    ! (advect_tracer):
    !   init_tracers_AB  -> zero del_ttf*, AB-interpolate valuesAB, rebuild edge gradient
    !   do_oce_adv_tra   -> del_ttf_advhoriz / del_ttf_advvert
    !   del_ttf += del_ttf_advhoriz + del_ttf_advvert
    ! velocities come from the dynamics state (UV at elements; full/implicit/explicit
    ! vertical velocity w/w_i/w_e at nodes), exactly as FESOM2 adv_tracers_ale:177-232.
    !
    ! SCOPE (M1.4): advection only. The rest of FESOM2 adv_tracers_ale's per-tracer
    ! body — diff_tracers_ale (horizontal+implicit-vertical diffusion + the tracer
    ! update values+=del_ttf), radioactive decay, relax_to_clim, the GM bolus-velocity
    ! add/subtract, SPP — enters at M2 (needs density/mixing/the SSH solve). The
    ! per-tracer del_ttf this routine assembles is the gated M1.4 output. 1-rank only;
    ! the FESOM2 exchange_nod(values) at the loop tail is a no-op here (M1.5).
    use mod_precision,      only: WP
    use mod_mesh,           only: t_mesh
    use mod_dyn,            only: t_dyn
    use mod_tracer,         only: t_tracer
    use oce_tracer_mod,     only: init_tracers_AB
    use oce_adv_tra_driver, only: do_oce_adv_tra
    implicit none
    private
    public :: adv_tracers_ale, advect_tracer

contains

    subroutine adv_tracers_ale(dt, dynamics, tracers, mesh)
        real(kind=WP),  intent(in)            :: dt
        type(t_mesh),   intent(in)            :: mesh
        type(t_dyn),    intent(inout), target :: dynamics
        type(t_tracer), intent(inout)         :: tracers
        integer :: tr_num
        do tr_num = 1, tracers%num_tracers
            call advect_tracer(dt, tr_num, dynamics, tracers, mesh)
            ! M2: diff_tracers_ale (incl. the implicit vertical solve + values+=del_ttf),
            !     decay, relax_to_clim, exchange_nod(values).
        end do
    end subroutine adv_tracers_ale

    subroutine advect_tracer(dt, tr_num, dynamics, tracers, mesh)
        ! One tracer's advection contribution to del_ttf (FESOM2 adv_tracers_ale body).
        real(kind=WP),  intent(in)            :: dt
        integer,        intent(in)            :: tr_num
        type(t_mesh),   intent(in)            :: mesh
        type(t_dyn),    intent(inout), target :: dynamics
        type(t_tracer), intent(inout)         :: tracers
        integer :: n

        call init_tracers_AB(tr_num, tracers, mesh)
        call do_oce_adv_tra(dt, dynamics%uv, dynamics%w, dynamics%w_i, dynamics%w_e, &
                            tr_num, dynamics, tracers, mesh)
        ! total tracer tendency = horizontal + vertical advection (del_ttf was zeroed
        ! in init_tracers_AB).
        do n = 1, mesh%nod2D
            tracers%work%del_ttf(:, n) = tracers%work%del_ttf(:, n) &
                                       + tracers%work%del_ttf_advhoriz(:, n) &
                                       + tracers%work%del_ttf_advvert(:, n)
        end do
    end subroutine advect_tracer

end module oce_ale_tracer
