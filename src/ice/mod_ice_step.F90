module mod_ice_step
    ! M3f: the assembled sea-ice time step. Faithful transcription of FESOM2 v2.7.3
    ! ice_setup_step.F90::ice_timestep (:96) — the per-step ice driver the main runloop
    ! calls (fesom_module.F90:719) between ocean2ice (:689) and oce_fluxes_mom/oce_fluxes
    ! (:727-728). It chains the M3b/M3c/M3d leaf kernels into one routine with LIVE data
    ! flow (each kernel reads the previous kernel's output):
    !
    !   EVPdynamics_solve   EVP/mEVP rheology -> uice/vice (+ prognostic sigma)   (M3b)
    !   ice_TG_rhs          Taylor-Galerkin advection rhs (ice_rhs_a/m/ms)        (M3c)
    !   ice_fct_solve       FCT solve (Zalesak limiter) -> advected a/m/m_snow     (M3c)
    !   cut_off             hmin/Armin clamp                                       (M3d)
    !   thermodynamics      0-layer Semtner growth -> a/m/m_snow + flx_h/flx_fw   (M3d)
    !
    ! Faithful-omission notes (NEITHER is consumed in the reduced M3 config, so omitting
    ! them is byte-neutral for the gated ocean substeps):
    !  * The CMIP6 dynamical growth-rate diagnostics (ice%thermo%dyngr*, ice_setup_step:203
    !    + :307) are write-only output diagnostics — nothing downstream reads them.
    !  * The post-thermo h_ice/h_snow effective-thickness diagnostics (ice_setup_step:320)
    !    are likewise output-only.
    !  * The cavity cleans (cavity_ice_clean_vel/_ma) are use_cavity-only (off).
    !
    ! ocean2ice (the ocean-surface read) and oce_fluxes_mom/oce_fluxes (the coupling-out)
    ! are SEPARATE runloop calls in FESOM2 — they stay separate in the lifecycle driver,
    ! exactly as the runloop has them; ice_timestep is only the 5-kernel chain above.
    !
    ! M2.12 optional-partit pattern: partit absent (or npes==1) -> each kernel runs its
    ! proven 1-rank path VERBATIM (the possibly-absent optional is passed straight through;
    ! Fortran treats a non-present optional actual as absent). Present+npes>1 -> the kernels
    ! use owned/halo bounds + their FESOM2 halo exchanges. No arithmetic here, only sequencing.
    use mod_mesh,        only: t_mesh
    use mod_partit,      only: t_partit
    use mod_ice,         only: t_ice
    use mod_ice_dyn,     only: EVPdynamics_solve
    use mod_ice_fct,     only: ice_TG_rhs, ice_fct_solve
    use mod_ice_thermo,  only: t_atmflux, cut_off, thermodynamics
    implicit none
    private
    public :: ice_timestep

contains

    subroutine ice_timestep(ice, mesh, atm, partit)
        type(t_ice),     intent(inout), target :: ice
        type(t_mesh),    intent(in),    target :: mesh
        type(t_atmflux), intent(inout), target :: atm
        type(t_partit),  intent(in),    optional :: partit

        !_______________________________________________________________________
        ! Dynamics: EVP (or mEVP) momentum solve -> uice/vice, prognostic sigma.
        call EVPdynamics_solve(ice, mesh, partit)

        !_______________________________________________________________________
        ! Advection: Taylor-Galerkin rhs + FCT solve -> advected a_ice/m_ice/m_snow.
        call ice_TG_rhs(ice, mesh, partit)
        call ice_fct_solve(ice, mesh, partit)

        !_______________________________________________________________________
        ! Clamp tiny ice/area to zero (hmin/Armin), then thermodynamic growth.
        call cut_off(ice, mesh, partit)
        call thermodynamics(ice, mesh, atm, partit)
    end subroutine ice_timestep

end module mod_ice_step
