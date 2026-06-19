module oce_adv_tra_hor
    ! Horizontal tracer advection — edge finite-volume fluxes, transcribed from
    ! FESOM2 v2.7.3 oce_adv_tra_hor.F90:
    !   adv_tra_hor_upw1  (64-257)  low-order upwind
    !   adv_tra_hor_muscl (261-542) MUSCL (3rd/4th-order, bottom-stable via nboundary_lay)
    !   adv_tra_hor_mfct  (546-834) MUSCL for the FCT path (no bottom-boundary clamp)
    !
    ! Each returns an EDGE flux that contributes with +sign to the 1st edge node and
    ! -sign to the 2nd (the driver scatters it to nodes). o_init_zero=.true. zeroes
    ! the flux first; .false. SUBTRACTS the new contribution from the input flux
    ! (so a HO call after an LO call yields the antidiffusive flux). flux is NOT
    ! multiplied by dt.
    !
    ! 1-rank only (myDim_edge2D == edge2D; el(2)==0 marks a boundary edge — its
    ! el(2)-using depth segments have empty ranges so el(2) is never indexed).
    ! Multi-rank loop bounds are deferred to M1.5. The volume flux scales the edge
    ! velocity by helem (layer thickness at element); a = mean r_earth*elem_cos over
    ! the edge's two elements (used by the MUSCL reconstruction, dead in upw1 — kept
    ! to mirror FESOM2).
    use mod_precision, only: WP
    use mod_constants, only: r_earth
    use mod_mesh,      only: t_mesh
    implicit none
    private
    public :: adv_tra_hor_upw1, adv_tra_hor_muscl, adv_tra_hor_mfct

contains

    !===========================================================================
    subroutine adv_tra_hor_upw1(vel, ttf, mesh, flux, o_init_zero)
        ! Low-order upwind horizontal flux (oce_adv_tra_hor.F90:64-257).
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: ttf(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: vel(2, mesh%nl-1, mesh%elem2D)
        real(kind=WP), intent(inout) :: flux(mesh%nl-1, mesh%edge2D)
        logical, optional, intent(in) :: o_init_zero
        logical       :: l_init_zero
        real(kind=WP) :: deltaX1, deltaY1, deltaX2, deltaY2, a, vflux
        integer       :: el(2), enodes(2), nz, edge, nu12, nl12, nl1, nl2, nu1, nu2

        l_init_zero = .true.
        if (present(o_init_zero)) l_init_zero = o_init_zero
        if (l_init_zero) then
            do edge = 1, mesh%edge2D
                do nz = 1, mesh%nl-1
                    flux(nz, edge) = 0.0_WP
                end do
            end do
        end if

        do edge = 1, mesh%edge2D
            enodes = mesh%edges(:, edge)
            el     = mesh%edge_tri(:, edge)
            nl1    = mesh%nlevels(el(1)) - 1
            nu1    = mesh%ulevels(el(1))
            deltaX1 = mesh%edge_cross_dxdy(1, edge)
            deltaY1 = mesh%edge_cross_dxdy(2, edge)
            a = r_earth * mesh%elem_cos(el(1))
            nl2 = 0; nu2 = 0
            if (el(2) > 0) then
                deltaX2 = mesh%edge_cross_dxdy(3, edge)
                deltaY2 = mesh%edge_cross_dxdy(4, edge)
                nl2 = mesh%nlevels(el(2)) - 1
                nu2 = mesh%ulevels(el(2))
                a = 0.5_WP * (a + r_earth * mesh%elem_cos(el(2)))
            end if
            nl12 = min(nl1, nl2)
            nu12 = max(nu1, nu2)
            ! (A) el(1)-only top segments (cavity surface)
            do nz = nu1, nu12-1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                call flux_lo(nz, vflux)
            end do
            ! (B) el(2)-only top segments (cavity surface)
            if (nu2 > 0) then
                do nz = nu2, nu12-1
                    vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                    call flux_lo(nz, vflux)
                end do
            end if
            ! (C) both segments
            do nz = nu12, nl12
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1)) &
                      + ( vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                call flux_lo(nz, vflux)
            end do
            ! (D) remaining el(1) segments
            do nz = nl12+1, nl1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                call flux_lo(nz, vflux)
            end do
            ! (E) remaining el(2) segments
            do nz = nl12+1, nl2
                vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                call flux_lo(nz, vflux)
            end do
        end do

    contains
        subroutine flux_lo(nz, vflux)
            integer,       intent(in) :: nz
            real(kind=WP), intent(in) :: vflux
            flux(nz, edge) = -0.5_WP * ( ttf(nz, enodes(1))*(vflux+abs(vflux)) &
                                       + ttf(nz, enodes(2))*(vflux-abs(vflux)) ) - flux(nz, edge)
        end subroutine flux_lo
    end subroutine adv_tra_hor_upw1

    !===========================================================================
    subroutine adv_tra_hor_muscl(vel, ttf, mesh, num_ord, flux, edge_up_dn_grad, nboundary_lay, o_init_zero)
        ! MUSCL horizontal flux (oce_adv_tra_hor.F90:261-542). num_ord = fraction of
        ! 4th-order (centered) contribution; (1-num_ord) is 3rd-order upwind. The
        ! per-node clamp c_lo = max(sign(1,nboundary_lay-nz),0) switches off the
        ! linear-reconstruction increment below a node's boundary layer.
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: num_ord
        real(kind=WP), intent(in)    :: ttf(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: vel(2, mesh%nl-1, mesh%elem2D)
        real(kind=WP), intent(inout) :: flux(mesh%nl-1, mesh%edge2D)
        integer,       intent(in)    :: nboundary_lay(mesh%nod2D)
        real(kind=WP), intent(in)    :: edge_up_dn_grad(4, mesh%nl-1, mesh%edge2D)
        logical, optional, intent(in) :: o_init_zero
        logical       :: l_init_zero
        real(kind=WP) :: deltaX1, deltaY1, deltaX2, deltaY2, a, vflux
        integer       :: el(2), enodes(2), nz, edge, nu12, nl12, nl1, nl2, nu1, nu2

        l_init_zero = .true.
        if (present(o_init_zero)) l_init_zero = o_init_zero
        if (l_init_zero) then
            do edge = 1, mesh%edge2D
                flux(:, edge) = 0.0_WP
            end do
        end if

        do edge = 1, mesh%edge2D
            enodes = mesh%edges(:, edge)
            el     = mesh%edge_tri(:, edge)
            nl1    = mesh%nlevels(el(1)) - 1
            nu1    = mesh%ulevels(el(1))
            deltaX1 = mesh%edge_cross_dxdy(1, edge)
            deltaY1 = mesh%edge_cross_dxdy(2, edge)
            a = r_earth * mesh%elem_cos(el(1))
            nl2 = 0; nu2 = 0
            if (el(2) > 0) then
                deltaX2 = mesh%edge_cross_dxdy(3, edge)
                deltaY2 = mesh%edge_cross_dxdy(4, edge)
                nl2 = mesh%nlevels(el(2)) - 1
                nu2 = mesh%ulevels(el(2))
                a = 0.5_WP * (a + r_earth * mesh%elem_cos(el(2)))
            end if
            nl12 = min(nl1, nl2)
            nu12 = max(nu1, nu2)
            ! (A)
            do nz = nu1, nu12-1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                call flux_ho(nz, clof(enodes(1),nz), clof(enodes(2),nz), vflux)
            end do
            ! (B)
            if (nu2 > 0) then
                do nz = nu2, nu12-1
                    vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                    call flux_ho(nz, clof(enodes(1),nz), clof(enodes(2),nz), vflux)
                end do
            end if
            ! (C)
            do nz = nu12, nl12
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1)) &
                      + ( vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                call flux_ho(nz, clof(enodes(1),nz), clof(enodes(2),nz), vflux)
            end do
            ! (D)
            do nz = nl12+1, nl1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                call flux_ho(nz, clof(enodes(1),nz), clof(enodes(2),nz), vflux)
            end do
            ! (E)
            do nz = nl12+1, nl2
                vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                call flux_ho(nz, clof(enodes(1),nz), clof(enodes(2),nz), vflux)
            end do
        end do

    contains
        real(kind=WP) function clof(node, nz)
            integer, intent(in) :: node, nz
            clof = real(max(sign(1, nboundary_lay(node)-nz), 0), WP)
        end function clof
        subroutine flux_ho(nz, clo1, clo2, vflux)
            integer,       intent(in) :: nz
            real(kind=WP), intent(in) :: clo1, clo2, vflux
            real(kind=WP) :: Tmean1, Tmean2, cHO
            Tmean2 = ttf(nz, enodes(2)) - &
                     (2.0_WP*(ttf(nz, enodes(2))-ttf(nz, enodes(1))) + &
                      mesh%edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge) + &
                      mesh%edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP*clo2
            Tmean1 = ttf(nz, enodes(1)) + &
                     (2.0_WP*(ttf(nz, enodes(2))-ttf(nz, enodes(1))) + &
                      mesh%edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge) + &
                      mesh%edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP*clo1
            cHO = (vflux+abs(vflux))*Tmean1 + (vflux-abs(vflux))*Tmean2
            flux(nz,edge) = -0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(Tmean1+Tmean2) - flux(nz,edge)
        end subroutine flux_ho
    end subroutine adv_tra_hor_muscl

    !===========================================================================
    subroutine adv_tra_hor_mfct(vel, ttf, mesh, num_ord, flux, edge_up_dn_grad, o_init_zero)
        ! MUSCL for the FCT path (oce_adv_tra_hor.F90:546-834). Same as
        ! adv_tra_hor_muscl but WITHOUT the c_lo bottom-boundary clamp (the
        ! reconstruction near bottom topography is not upwind; runs with FCT only).
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: num_ord
        real(kind=WP), intent(in)    :: ttf(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: vel(2, mesh%nl-1, mesh%elem2D)
        real(kind=WP), intent(inout) :: flux(mesh%nl-1, mesh%edge2D)
        real(kind=WP), intent(in)    :: edge_up_dn_grad(4, mesh%nl-1, mesh%edge2D)
        logical, optional, intent(in) :: o_init_zero
        logical       :: l_init_zero
        real(kind=WP) :: deltaX1, deltaY1, deltaX2, deltaY2, a, vflux
        integer       :: el(2), enodes(2), nz, edge, nu12, nl12, nl1, nl2, nu1, nu2

        l_init_zero = .true.
        if (present(o_init_zero)) l_init_zero = o_init_zero
        if (l_init_zero) then
            do edge = 1, mesh%edge2D
                flux(:, edge) = 0.0_WP
            end do
        end if

        do edge = 1, mesh%edge2D
            enodes = mesh%edges(:, edge)
            el     = mesh%edge_tri(:, edge)
            nl1    = mesh%nlevels(el(1)) - 1
            nu1    = mesh%ulevels(el(1))
            deltaX1 = mesh%edge_cross_dxdy(1, edge)
            deltaY1 = mesh%edge_cross_dxdy(2, edge)
            a = r_earth * mesh%elem_cos(el(1))
            nl2 = 0; nu2 = 0
            if (el(2) > 0) then
                deltaX2 = mesh%edge_cross_dxdy(3, edge)
                deltaY2 = mesh%edge_cross_dxdy(4, edge)
                nl2 = mesh%nlevels(el(2)) - 1
                nu2 = mesh%ulevels(el(2))
                a = 0.5_WP * (a + r_earth * mesh%elem_cos(el(2)))
            end if
            nl12 = min(nl1, nl2)
            nu12 = max(nu1, nu2)
            ! (A)
            do nz = nu1, nu12-1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                call flux_ho(nz, vflux)
            end do
            ! (B)
            if (nu2 > 0) then
                do nz = nu2, nu12-1
                    vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                    call flux_ho(nz, vflux)
                end do
            end if
            ! (C)
            do nz = nu12, nl12
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1)) &
                      + ( vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                call flux_ho(nz, vflux)
            end do
            ! (D)
            do nz = nl12+1, nl1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                call flux_ho(nz, vflux)
            end do
            ! (E)
            do nz = nl12+1, nl2
                vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                call flux_ho(nz, vflux)
            end do
        end do

    contains
        subroutine flux_ho(nz, vflux)
            integer,       intent(in) :: nz
            real(kind=WP), intent(in) :: vflux
            real(kind=WP) :: Tmean1, Tmean2, cHO
            Tmean2 = ttf(nz, enodes(2)) - &
                     (2.0_WP*(ttf(nz, enodes(2))-ttf(nz, enodes(1))) + &
                      mesh%edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge) + &
                      mesh%edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP
            Tmean1 = ttf(nz, enodes(1)) + &
                     (2.0_WP*(ttf(nz, enodes(2))-ttf(nz, enodes(1))) + &
                      mesh%edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge) + &
                      mesh%edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP
            cHO = (vflux+abs(vflux))*Tmean1 + (vflux-abs(vflux))*Tmean2
            flux(nz,edge) = -0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(Tmean1+Tmean2) - flux(nz,edge)
        end subroutine flux_ho
    end subroutine adv_tra_hor_mfct

end module oce_adv_tra_hor
