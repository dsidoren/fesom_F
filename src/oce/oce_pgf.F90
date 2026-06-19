module oce_pgf
    ! Hydrostatic pressure gradient force (PGF) for which_ALE='linfs' with full
    ! cells. Transcribed from FESOM2 v2.7.3 oce_ale_pressure_bv.F90:
    !   pressure_force_4_linfs          (507-569)  dispatcher: branches on
    !                                              use_partial_cell/use_cavity/which_pgf
    !   pressure_force_4_linfs_fullcell (575-614)  the full-cell PGF kernel (ported)
    !
    ! M2.2 SCOPE — consumes the M2.1 hpressure (top-down hydrostatic integral) and
    ! contracts it with gradient_sca (the per-element linear shape-function
    ! derivatives, geometry-gated) into the element PGF pgf_x/pgf_y — the horizontal
    ! pressure term of the momentum RHS. NO hydrostatic integration here (that lives
    ! in oce_pressure_bv / M2.1). The dispatcher's non-fullcell branches (nemo,
    ! shchepetkin, cubicspline, easypgf, cavity, partial cells) are v1-dropped and
    ! UNGATED on pi (use_partial_cell=.false., no cavity ⇒ the fullcell branch); the
    ! dispatcher itself is pure namelist branching (no numerics) so it is not ported.
    !
    !   pgf_x(nz,elem) = Σ_k gradient_sca(k,  elem) * hpressure(nz, elnodes(k)) / density_0
    !   pgf_y(nz,elem) = Σ_k gradient_sca(3+k,elem) * hpressure(nz, elnodes(k)) / density_0
    !
    ! BIT-IDENTITY NOTES:
    !  - density_0 is a RUNTIME divisor, but byte-identical on both sides (= 1030),
    !    so the -no-prec-div reciprocal matches (L7/L10). Transcribe the expression
    !    VERBATIM: the /density_0 is INSIDE the sum (each of the 3 products divided,
    !    then summed) — NOT sum(...)/density_0.
    !  - gradient_sca and the elem2D_nodes order are geometry-gate-proven (L9); the
    !    same gradient_sca contraction already byte-matched in oce_tracer_grad (M1.1).
    !  - v1 = triangles only (the gradient_sca 1:3 / 4:6 slices are triangle-specific).
    !
    ! 1-rank only: the element loop runs to elem2D (== myDim_elem2D at 1-rank);
    ! pgf_x/pgf_y are element-local (no halo dependency, FESOM2 does not exchange
    ! them). The outputs are intent(inout) and NOT zeroed inside (FESOM2 leaves
    ! below-bottom entries at their allocation-zero) — the CALLER pre-zeros them.
    use mod_precision, only: WP
    use mod_mesh,      only: t_mesh
    use mod_constants, only: density_0
    implicit none
    private
    public :: pressure_force_4_linfs_fullcell

contains

    !===========================================================================
    subroutine pressure_force_4_linfs_fullcell(hpressure, mesh, pgf_x, pgf_y)
        ! hpressure: (nl, nod2D) in (M2.1 output; only levels 1..nlevels-1 read).
        ! pgf_x/pgf_y: (nl-1, elem2D) inout; the caller pre-zeros them (see header).
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: hpressure(mesh%nl, mesh%nod2D)
        real(kind=WP), intent(inout) :: pgf_x(mesh%nl-1, mesh%elem2D)
        real(kind=WP), intent(inout) :: pgf_y(mesh%nl-1, mesh%elem2D)
        integer :: elem, elnodes(3), nle, ule, nlz

        !_______________________________________________________________________
        ! loop over triangular elements
        do elem=1, mesh%elem2D
            !___________________________________________________________________
            ! number of levels at elem
            nle = mesh%nlevels(elem)-1
            ule = mesh%ulevels(elem)

            !___________________________________________________________________
            ! node indices of elem
            elnodes = mesh%elem2D_nodes(1:3, elem)

            !___________________________________________________________________
            ! loop over mid-depth levels to calculate the pressure gradient force
            ! (pgf) --> from top to bottom
            do nlz=ule,nle
                pgf_x(nlz,elem) = sum(mesh%gradient_sca(1:3,elem)*hpressure(nlz,elnodes)/density_0)
                pgf_y(nlz,elem) = sum(mesh%gradient_sca(4:6,elem)*hpressure(nlz,elnodes)/density_0)
            end do
        end do
    end subroutine pressure_force_4_linfs_fullcell

end module oce_pgf
