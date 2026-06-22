module oce_adv_tra_ver
    ! Vertical tracer advection — interface fluxes, transcribed from FESOM2 v2.7.3
    ! oce_adv_tra_ver.F90:
    !   adv_tra_ver_upw1 (244-328)  1st-order upwind (explicit)
    !   adv_tra_ver_qr4c (332-434)  QR 4th-order centered (num_ord = 4th-order fraction)
    !
    ! Each returns a flux given at the vertical interfaces of the scalar volumes,
    ! flux(nz,node) with nz = 1..nl interfaces (surface = ulevels_nod2D, zero at the
    ! bottom interface nlevels_nod2D). o_init_zero=.true. zeroes the flux first;
    ! .false. SUBTRACTS the new contribution from the input flux (so an HO call after
    ! an LO call yields the antidiffusive flux). flux is NOT multiplied by dt — the
    ! driver scatters it as (flux(nz)-flux(nz+1))*dt/areasvol (oce_adv_tra_flux).
    !
    ! area(nz,n) is the geom-proven node control-volume area. Z_3d_n/zbar_3d_n are the
    ! ALE per-node mid-depth / interface-depth arrays; at the initial state (no cavity,
    ! full cells) the caller builds them from FESOM2 init_ale (oce_ale.F90:531-566) and
    ! the gate verifies them. QR4C's qc/qu/qd divide by (Z_3d_n(k)-Z_3d_n(k+1)), a
    ! RUNTIME divisor — FESOM2 does the identical division on byte-identical operands,
    ! so the -no-prec-div reciprocal matches (cf. LESSONS L7: the trap is a runtime vs
    ! LITERAL divisor mismatch, not a runtime divisor per se).
    !
    ! 1-rank only (myDim_nod2D == nod2D). On pi the shallowest column is 4 layers, so
    ! neither the 1-layer ttf(0) read nor the QR4C 2-layer double-write (2nd-layer and
    ! bottom-1 both hitting interface nzmin+1) occurs; the statement order below still
    ! reproduces both faithfully for deeper-min meshes.
    use mod_precision,   only: WP
    use mod_mesh,        only: t_mesh
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: owned_bounds
    implicit none
    private
    public :: adv_tra_ver_upw1, adv_tra_ver_qr4c

contains

    !===========================================================================
    subroutine adv_tra_ver_upw1(w, ttf, mesh, flux, o_init_zero, partit)
        ! 1st-order upwind explicit vertical flux (oce_adv_tra_ver.F90:244-328).
        ! M2.12b: optional partit -> loop over OWNED nodes (myDim_nod2D).
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: ttf(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: w  (mesh%nl,   mesh%nod2D)
        real(kind=WP), intent(inout) :: flux(mesh%nl,  mesh%nod2D)
        logical, optional, intent(in) :: o_init_zero
        type(t_partit), intent(in), optional :: partit
        logical :: l_init_zero
        integer :: n, nz, nzmax, nzmin
        integer :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        l_init_zero = .true.
        if (present(o_init_zero)) l_init_zero = o_init_zero
        if (l_init_zero) then
            do n = 1, nNodO
                do nz = 1, mesh%nl
                    flux(nz, n) = 0.0_WP
                end do
            end do
        end if

        do n = 1, nNodO
            nzmax = mesh%nlevels_nod2D(n)
            nzmin = mesh%ulevels_nod2D(n)
            ! vert. flux at surface layer
            nz = nzmin
            flux(nz,n) = -w(nz,n)*ttf(nz,n)*mesh%area(nz,n) - flux(nz,n)
            ! vert. flux at bottom layer --> zero bottom flux
            nz = nzmax
            flux(nz,n) = 0.0_WP - flux(nz,n)
            ! vert. flux at remaining levels (upwind)
            do nz = nzmin+1, nzmax-1
                flux(nz,n) = -0.5*(                                            &
                              ttf(nz  ,n)*(w(nz,n)+abs(w(nz,n))) +             &
                              ttf(nz-1,n)*(w(nz,n)-abs(w(nz,n))))*mesh%area(nz,n) - flux(nz,n)
            end do
        end do
    end subroutine adv_tra_ver_upw1

    !===========================================================================
    subroutine adv_tra_ver_qr4c(w, ttf, mesh, num_ord, flux, o_init_zero, partit)
        ! QR 4th-order centered vertical flux (oce_adv_tra_ver.F90:332-434). num_ord =
        ! fraction of the 4th-order (centered) contribution; (1-num_ord) is upwind-QR.
        ! M2.12b: optional partit -> loop over OWNED nodes (myDim_nod2D).
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: num_ord
        real(kind=WP), intent(in)    :: ttf(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: w  (mesh%nl,   mesh%nod2D)
        real(kind=WP), intent(inout) :: flux(mesh%nl,  mesh%nod2D)
        logical, optional, intent(in) :: o_init_zero
        type(t_partit), intent(in), optional :: partit
        logical :: l_init_zero
        integer :: n, nz, nzmax, nzmin
        real(kind=WP) :: Tmean, Tmean1, Tmean2, qc, qu, qd
        integer :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        l_init_zero = .true.
        if (present(o_init_zero)) l_init_zero = o_init_zero
        if (l_init_zero) then
            do n = 1, nNodO
                do nz = 1, mesh%nl
                    flux(nz, n) = 0.0_WP
                end do
            end do
        end if

        do n = 1, nNodO
            nzmax = mesh%nlevels_nod2D(n)
            nzmin = mesh%ulevels_nod2D(n)
            ! vert. flux at surface layer
            nz = nzmin
            flux(nz,n) = -ttf(nz,n)*w(nz,n)*mesh%area(nz,n) - flux(nz,n)
            ! vert. flux 2nd layer --> centered differences
            nz = nzmin+1
            flux(nz,n) = -0.5_WP*(ttf(nz-1,n)+ttf(nz,n))*w(nz,n)*mesh%area(nz,n) - flux(nz,n)
            ! vert. flux at bottom-1 layer --> centered differences
            nz = nzmax-1
            flux(nz,n) = -0.5_WP*(ttf(nz-1,n)+ttf(nz,n))*w(nz,n)*mesh%area(nz,n) - flux(nz,n)
            ! vert. flux at bottom layer --> zero bottom flux
            nz = nzmax
            flux(nz,n) = 0.0_WP - flux(nz,n)
            ! vert. flux at remaining levels (4th order)
            do nz = nzmin+2, nzmax-2
                qc = (ttf(nz-1,n)-ttf(nz  ,n))/(mesh%Z_3d_n(nz-1,n)-mesh%Z_3d_n(nz  ,n))
                qu = (ttf(nz  ,n)-ttf(nz+1,n))/(mesh%Z_3d_n(nz  ,n)-mesh%Z_3d_n(nz+1,n))
                qd = (ttf(nz-2,n)-ttf(nz-1,n))/(mesh%Z_3d_n(nz-2,n)-mesh%Z_3d_n(nz-1,n))

                Tmean1 = ttf(nz  ,n)+(2*qc+qu)*(mesh%zbar_3d_n(nz,n)-mesh%Z_3d_n(nz  ,n))/3.0_WP
                Tmean2 = ttf(nz-1,n)+(2*qc+qd)*(mesh%zbar_3d_n(nz,n)-mesh%Z_3d_n(nz-1,n))/3.0_WP
                Tmean  = (w(nz,n)+abs(w(nz,n)))*Tmean1+(w(nz,n)-abs(w(nz,n)))*Tmean2
                flux(nz,n) = (-0.5_WP*(1.0_WP-num_ord)*Tmean - num_ord*(0.5_WP*(Tmean1+Tmean2))*w(nz,n))*mesh%area(nz,n) - flux(nz,n)
            end do
        end do
    end subroutine adv_tra_ver_qr4c

end module oce_adv_tra_ver
