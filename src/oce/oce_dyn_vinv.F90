module oce_dyn_vinv
    ! Vector-invariant momentum advection (momadv_opt==1):
    !
    !     (curl u + f) x u  +  grad(u^2/2)  +  w du/dz
    !
    ! PROVENANCE. Transcribed from /home/a/a270029/qq/oce_vinv_mom_adv.F90 (routines
    ! relative_vorticity + v_inv_mom_adv). There is NO FESOM2 oracle for this path: the
    ! v2.7.3 branch exists but ABORTS --
    !   oce_ale_vel_rhs.F90:268  'in moment not adapted mom_adv advection typ for ALE'
    ! so FESOM2 never adapted the vector-invariant form to ALE and it was never tested
    ! there. Correctness therefore rests on analytic unit tests (test/test_vinv.F90) plus
    ! the no-work invariant u.[(f+zeta) x u] == 0 checked in fesom_conserve, NOT on a
    ! reference dump. See docs/plans/2026-10-01-momadv-vector-invariant.md.
    !
    ! INTEGRATION. momentum_adv_vinv is a pure ADDITION into UV_rhsAB, exactly like
    ! momentum_adv_scalar. compute_vel_rhs step (3) already wrote f x u into that slot, so
    ! this routine adds  zeta x u - grad(KE) + w du/dz, which completes (f+zeta) x u.
    ! compute_vel_rhs steps (1)/(2)/(4) are untouched and momadv_opt==2 is unreachable
    ! from here, so that path stays bit-identical.
    !
    ! FORCED DEVIATIONS from qq (it is pre-ALE research code):
    !   qq                                   -> here                  why
    !   dz = zbar(1:nl-1)-zbar(2:nl) global  -> mesh%helem(nz,elem)   ALE: per-elem/step
    !   area(nz,n)                           -> mesh%area(n)          collapsed to 1-D
    !   w_cv(:,elem)                         -> 1/3                   absent; triangles
    !   elnodes(4), gradient_sca(1:4)/(5:8)  -> (1:3), (1:3)/(4:6)    MAX_NV=4 (L15 trap)
    !   exchange_nod3D(x)                    -> exchange_nod(x,partit) halo API
    !   loops from 1                         -> ulevels..nlevels-1    cavity bounds
    !   Wvel                                 -> dynamics%w_e          match momadv_opt==2
    !   whole-RHS assembly                   -> advection blocks only rest is in velrhs
    ! qq's `else` branch (explicit vertical viscosity: wind stress + a hard-coded
    ! friction=0.005 "Soufflet" replacing the real bottom drag) is NOT transcribed; we
    ! take its i_vert_visc=.true. branch, whose own comment reads "Do only advection".
    use mod_precision,   only: WP
    use mod_mesh,        only: t_mesh
    use mod_dyn,         only: t_dyn
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: owned_bounds, is_multirank
    use mod_halo,        only: exchange_nod
    implicit none
    private
    public :: relative_vorticity, momentum_adv_vinv

contains

    subroutine relative_vorticity(dynamics, mesh, partit)
        ! Block A. zeta = curl(u) at nodes, as the circulation integral around the median-
        ! dual contour (qq oce_vinv_mom_adv.F90:5-61).
        !
        ! Per edge the contour crosses the two adjacent element centres, so the segment
        ! contribution is  d . u  with d = edge_cross_dxdy (METRES here, as qq assumes),
        ! taken with opposite signs for the two edge nodes. qq's THREE level ranges are
        ! kept verbatim: both elements wet, then el1 only, then el2 only. Those one-sided
        ! ranges are what makes the operator correct at a bathymetry step -- below the
        ! shallower element only the deeper one contributes, i.e. the dry element's velocity
        ! is treated as exactly zero while area(n) keeps its FULL weight. That is precisely
        ! the bottom-at-vertices rule (mod_mesh.F90), so area(nz,n) -> area(n) and the
        ! one-sided ranges are mutually consistent, not two independent choices.
        !
        ! Halo: every edge incident to an owned node is owned, and both triangles of an
        ! owned edge are owned (LESSONS L31), so the loop over 1..nEdgeO is complete for
        ! owned rows and UV needs NO element exchange. One exchange_nod publishes zeta to
        ! the halo, which the per-element (f+zeta) average then reads.
        !
        ! Cavity: the TOP asymmetry (one adjacent element shallower at the surface) is NOT
        ! handled -- the ranges below start at the deeper ulevels. Untested without a
        ! cavity mesh, like the rest of the tree's cavity paths.
        type(t_mesh),   intent(in)              :: mesh
        type(t_dyn),    intent(inout), target   :: dynamics
        type(t_partit), intent(in),    optional :: partit
        integer       :: n, nz, edge, el(2), enodes(2), ul1, ul2, nl1, nl2, nboth, ultop
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: dx1, dy1, dx2, dy2, c1
        real(kind=WP), dimension(:,:),   pointer :: vort
        real(kind=WP), dimension(:,:,:), pointer :: UV

        vort => dynamics%work%vorticity
        UV   => dynamics%uv
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        ! Zero the FULL local range -- owned AND halo, all levels. The scatter below writes
        ! into halo rows, and -init=zero does not cover allocatables, so accumulating onto
        ! fresh heap could trap under -fpe0 before the exchange made the value irrelevant.
        vort(:, 1:nNodL) = 0.0_WP

        do edge = 1, nEdgeO
            enodes = mesh%edges(:, edge)
            el     = mesh%edge_tri(:, edge)
            ul1    = mesh%ulevels(el(1))
            nl1    = mesh%nlevels(el(1)) - 1
            dx1    = real(mesh%edge_cross_dxdy(1, edge), WP)
            dy1    = real(mesh%edge_cross_dxdy(2, edge), WP)
            ul2    = ul1
            nl2    = 0
            dx2    = 0.0_WP
            dy2    = 0.0_WP
            if (el(2) > 0) then
                dx2 = real(mesh%edge_cross_dxdy(3, edge), WP)
                dy2 = real(mesh%edge_cross_dxdy(4, edge), WP)
                ul2 = mesh%ulevels(el(2))
                nl2 = mesh%nlevels(el(2)) - 1
            end if
            nboth = min(nl1, nl2)
            ultop = max(ul1, ul2)

            ! both elements wet
            do nz = ultop, nboth
                c1 = dx1*UV(1,nz,el(1)) + dy1*UV(2,nz,el(1)) &
                   - dx2*UV(1,nz,el(2)) - dy2*UV(2,nz,el(2))
                vort(nz, enodes(1)) = vort(nz, enodes(1)) + c1
                vort(nz, enodes(2)) = vort(nz, enodes(2)) - c1
            end do
            ! el1 only (el2 dry or absent): its velocity contributes zero
            do nz = max(nboth+1, ul1), nl1
                c1 = dx1*UV(1,nz,el(1)) + dy1*UV(2,nz,el(1))
                vort(nz, enodes(1)) = vort(nz, enodes(1)) + c1
                vort(nz, enodes(2)) = vort(nz, enodes(2)) - c1
            end do
            ! el2 only
            do nz = max(nboth+1, ul2), nl2
                c1 = -dx2*UV(1,nz,el(2)) - dy2*UV(2,nz,el(2))
                vort(nz, enodes(1)) = vort(nz, enodes(1)) + c1
                vort(nz, enodes(2)) = vort(nz, enodes(2)) - c1
            end do
        end do

        ! circulation -> vorticity.
        !   area(n): the full median-dual area, sum over the node's elements of
        !            elem_area/3 (DEVIATION: qq divides by area(nz,n)).
        !   SIGN   : qq's combination c1 = d1.u(el1) - d2.u(el2), added to enodes(1) and
        !            subtracted from enodes(2), traverses the dual contour in the sense that
        !            matches OUR edge_tri left/right convention -- no flip needed.
        !            VERIFIED: test_vinv V2 gives zeta == -alpha for u = alpha*y to
        !            relative 5e-13 over 81729 node-levels, i.e. exactly
        !            curl(u) = dv/dx - du/dy. That same exactness also PROVES the
        !            area(nz,n) -> area(n) deviation is right: zeta = circulation/area
        !            only reproduces -alpha if area(n) equals the sum of elem_area/3 over
        !            the node's WET adjacent elements.
        do n = 1, nNodO
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n)-1
                vort(nz,n) = vort(nz,n)/real(mesh%area(n), WP)
            end do
        end do
        if (is_multirank(partit)) call exchange_nod(vort, partit)
    end subroutine relative_vorticity

    subroutine momentum_adv_vinv(dynamics, mesh, partit)
        ! STUB (Tasks 3-5). Adds zeta x u - grad(KE) + w du/dz into UV_rhsAB(1,1:2,:).
        ! Deliberately a no-op until the blocks land, so Task 1 can prove that adding the
        ! momadv_opt==1 dispatch leaves momadv_opt==2 bit-identical.
        type(t_mesh),   intent(in)              :: mesh
        type(t_dyn),    intent(inout), target   :: dynamics
        type(t_partit), intent(in),    optional :: partit
    end subroutine momentum_adv_vinv

end module oce_dyn_vinv
