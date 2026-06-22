module oce_tracer_grad
    ! Elemental tracer gradient tr_xy, transcribed from FESOM2 v2.7.3
    ! oce_tracer_mod.F90:149-190 (SUBROUTINE tracer_gradient_elements).
    !
    ! tr_xy(1:2, nz, elem) is the horizontal gradient of the tracer field ttf,
    ! reconstructed per element from the linear shape-function coefficients
    ! gradient_sca (built in mod_mesh_areas). It feeds fill_up_dn_grad (MUSCL).
    !
    ! Clean-architecture change vs FESOM2 (D7, plan "USE-globals -> explicit type
    ! arguments"): FESOM2 writes the module-global tr_xy (o_ARRAYS); here tr_xy is
    ! an explicit intent(out) argument. 1-rank only (myDim_elem2D == elem2D); the
    ! multi-rank loop bounds + halo exchange of tr_xy are deferred to M1.5/M2.12.
    use mod_precision,    only: WP
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_part_bounds,  only: owned_bounds
    implicit none
    private
    public :: tracer_gradient_elements

contains

    subroutine tracer_gradient_elements(ttf, tr_xy, mesh, partit)
        ! computes elemental gradient of tracer  (oce_tracer_mod.F90:149-190).
        ! M2.12b: optional partit -> loop over OWNED elements (the caller exchanges
        ! tr_xy to the full halo afterwards). 1-rank (partit absent) loops all elements.
        type(t_mesh),  intent(in)  :: mesh
        real(kind=WP), intent(in)  :: ttf(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(out) :: tr_xy(2, mesh%nl-1, mesh%elem2D)
        type(t_partit), intent(in), optional :: partit
        integer :: elem, elnodes(3), nz, nzmin, nzmax
        integer :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        do elem = 1, nElemO
            elnodes = mesh%elem2D_nodes(1:3, elem)
            nzmin = mesh%ulevels(elem)
            nzmax = mesh%nlevels(elem)
            do nz = nzmin, nzmax-1
                tr_xy(1, nz, elem) = sum(mesh%gradient_sca(1:3, elem) * ttf(nz, elnodes))
                tr_xy(2, nz, elem) = sum(mesh%gradient_sca(4:6, elem) * ttf(nz, elnodes))
            end do
        end do
    end subroutine tracer_gradient_elements

end module oce_tracer_grad
