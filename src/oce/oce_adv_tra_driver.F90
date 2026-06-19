module oce_adv_tra_driver
    ! Tracer-advection driver, transcribed from FESOM2 v2.7.3
    ! oce_adv_tra_driver.F90:46-419 (subroutine do_oce_adv_tra).
    !
    ! Assembles the byte-proven kernels (M1.1 horizontal, M1.2 vertical, M1.3 FCT
    ! limiter, the edge->node scatter) into the per-tracer advection operator. For a
    ! tracer with tra_adv_lim=='FCT' it runs the Zalesak path:
    !   LO = UPW1(values) horiz + vert            -> fct_LO  (low-order solution)
    !   HO = {hor scheme, ver scheme}(valuesAB)   -> antidiffusive flux (HO-LO)
    !   oce_tra_adv_fct clips the antidiffusive flux
    !   oce_tra_adv_flux2dtracer(use_lo) scatters -> del_ttf_advhoriz/advvert
    ! For a non-FCT tracer (do_zero_flux=.true.) the HO scheme is applied directly to
    ! valuesAB (o_init_zero=.true.) and scattered without the low-order reconstruction.
    !
    ! Dispatch is on the per-tracer scheme selectors (tracers%data(tr_num)%tra_adv_*):
    !   tra_adv_hor : UPW1 / MUSCL / MFCT
    !   tra_adv_ver : UPW1 / QR4C            (CDIFF / PPM declared but not yet ported)
    !   tra_adv_lim : FCT or anything else (non-FCT)
    ! tra_adv_ph / tra_adv_pv are the horizontal / vertical 4th-order fractions
    ! (opth / optv) passed to MUSCL/MFCT/QR4C.
    !
    ! SCOPE / clean-architecture deviations (D7, "USE-globals -> explicit args"):
    !  * 1-rank only (myDim_* == global). The two FESOM2 halo exchanges
    !    (exchange_nod(fct_LO) after the LO build, the implicit ones in the kernels)
    !    are no-ops at 1 rank and are dropped here; multi-rank is M1.5.
    !  * use_wsplit (split implicit/explicit vertical velocity) needs adv_tra_vert_impl
    !    (not ported until M2); the branch is guarded with a clear error. pi runs
    !    use_wsplit=.false., so w == w_e and the LO/HO vertical both use the same field.
    !  * the DVD (downgradient-variance-decomposition) diagnostic blocks are dropped (not v1).
    !  * work arrays live in tracers%work (MP); kernels take WP. At the DP/SP anchor
    !    MP==WP so the pointers/args bind by kind value; revisit for FP16 (D3).
    use mod_precision,    only: WP
    use mod_mesh,         only: t_mesh
    use mod_dyn,          only: t_dyn
    use mod_tracer,       only: t_tracer
    use oce_adv_tra_hor,  only: adv_tra_hor_upw1, adv_tra_hor_muscl, adv_tra_hor_mfct
    use oce_adv_tra_ver,  only: adv_tra_ver_upw1, adv_tra_ver_qr4c
    use oce_adv_tra_fct,  only: oce_tra_adv_fct
    use oce_adv_tra_flux, only: oce_tra_adv_flux2dtracer
    implicit none
    private
    public :: do_oce_adv_tra

contains

    subroutine do_oce_adv_tra(dt, vel, w, wi, we, tr_num, dynamics, tracers, mesh)
        ! do_oce_adv_tra (oce_adv_tra_driver.F90:46-419). vel = horizontal velocity at
        ! elements; w / we / wi = full / explicit / implicit vertical velocity at nodes.
        real(kind=WP),  intent(in)            :: dt
        type(t_mesh),   intent(in)            :: mesh
        type(t_dyn),    intent(in)            :: dynamics
        type(t_tracer), intent(inout), target :: tracers
        integer,        intent(in)            :: tr_num
        real(kind=WP),  intent(in)            :: vel(2, mesh%nl-1, mesh%elem2D)
        real(kind=WP),  intent(in), target    :: w (mesh%nl, mesh%nod2D)
        real(kind=WP),  intent(in)            :: wi(mesh%nl, mesh%nod2D)
        real(kind=WP),  intent(in), target    :: we(mesh%nl, mesh%nod2D)

        ! pointers into the per-tracer data / shared work (mirrors FESOM2 lines 92-106)
        real(kind=WP), pointer :: ttf(:,:), ttfAB(:,:), fct_LO(:,:)
        real(kind=WP), pointer :: adv_flux_hor(:,:), adv_flux_ver(:,:), dttf_h(:,:), dttf_v(:,:)
        real(kind=WP), pointer :: fct_ttf_min(:,:), fct_ttf_max(:,:), fct_plus(:,:), fct_minus(:,:)
        real(kind=WP), pointer :: edge_up_dn_grad(:,:,:)
        integer,       pointer :: nboundary_lay(:)
        real(kind=WP), pointer :: pwvel(:,:)

        real(kind=WP) :: opth, optv
        logical       :: do_zero_flux
        integer       :: e, n, nz, enodes(2), el(2), nl1, nl2, nu1, nu2, nl12, nu12

        ttf             => tracers%data(tr_num)%values
        ttfAB           => tracers%data(tr_num)%valuesAB
        opth            =  tracers%data(tr_num)%tra_adv_ph
        optv            =  tracers%data(tr_num)%tra_adv_pv
        fct_LO          => tracers%work%fct_LO
        adv_flux_ver    => tracers%work%adv_flux_ver
        adv_flux_hor    => tracers%work%adv_flux_hor
        edge_up_dn_grad => tracers%work%edge_up_dn_grad
        nboundary_lay   => tracers%work%nboundary_lay
        fct_ttf_min     => tracers%work%fct_ttf_min
        fct_ttf_max     => tracers%work%fct_ttf_max
        fct_plus        => tracers%work%fct_plus
        fct_minus       => tracers%work%fct_minus
        dttf_h          => tracers%work%del_ttf_advhoriz
        dttf_v          => tracers%work%del_ttf_advvert

        !_______________________________________________________________________
        ! FCT: low-order horizontal+vertical solution + low-order antidiffusive part
        if (trim(tracers%data(tr_num)%tra_adv_lim) == 'FCT') then
            ! low-order upwind horizontal flux (zero the flux first)
            call adv_tra_hor_upw1(vel, ttf, mesh, adv_flux_hor, o_init_zero=.true.)

            ! fct_LO = scatter of the LO horizontal flux over edges
            do n = 1, mesh%nod2D
                do nz = 1, mesh%nl-1
                    fct_LO(nz,n) = 0.0_WP
                end do
            end do
            do e = 1, mesh%edge2D
                enodes = mesh%edges(:, e)
                el     = mesh%edge_tri(:, e)
                nl1 = mesh%nlevels(el(1))-1
                nu1 = mesh%ulevels(el(1))
                nl2 = 0; nu2 = 0
                if (el(2) > 0) then
                    nl2 = mesh%nlevels(el(2))-1
                    nu2 = mesh%ulevels(el(2))
                end if
                nl12 = max(nl1, nl2)
                nu12 = nu1
                if (nu2 > 0) nu12 = min(nu1, nu2)
                do nz = nu12, nl12
                    fct_LO(nz, enodes(1)) = fct_LO(nz, enodes(1)) + adv_flux_hor(nz, e)
                    fct_LO(nz, enodes(2)) = fct_LO(nz, enodes(2)) - adv_flux_hor(nz, e)
                end do
            end do

            ! low-order upwind vertical flux (explicit part: we), then finish the LO solution
            call adv_tra_ver_upw1(we, ttf, mesh, adv_flux_ver, o_init_zero=.true.)
            do n = 1, mesh%nod2D
                nu1 = mesh%ulevels_nod2D(n)
                nl1 = mesh%nlevels_nod2D(n)
                do nz = nu1, nl1-1
                    fct_LO(nz,n) = (ttf(nz,n)*mesh%hnode(nz,n) &
                                  + (fct_LO(nz,n) + (adv_flux_ver(nz,n)-adv_flux_ver(nz+1,n)))*dt/mesh%areasvol(nz,n)) &
                                  / mesh%hnode_new(nz,n)
                end do
            end do

            if (dynamics%use_wsplit) then
                ! implicit (w-split) vertical correction needs adv_tra_vert_impl — M2.
                error stop 'do_oce_adv_tra: use_wsplit=.true. not yet supported (needs adv_tra_vert_impl, M2)'
            end if
            ! exchange_nod(fct_LO) — no-op at 1 rank (M1.5 for multi-rank)
        end if

        do_zero_flux = .true.
        if (trim(tracers%data(tr_num)%tra_adv_lim) == 'FCT') do_zero_flux = .false.

        !_______________________________________________________________________
        ! horizontal advection (FCT: high-order antidiffusive; else: full high-order)
        select case (trim(tracers%data(tr_num)%tra_adv_hor))
        case ('MUSCL')
            call adv_tra_hor_muscl(vel, ttfAB, mesh, opth, adv_flux_hor, edge_up_dn_grad, nboundary_lay, o_init_zero=do_zero_flux)
        case ('MFCT')
            call adv_tra_hor_mfct (vel, ttfAB, mesh, opth, adv_flux_hor, edge_up_dn_grad,                o_init_zero=do_zero_flux)
        case ('UPW1')
            call adv_tra_hor_upw1 (vel, ttfAB, mesh,       adv_flux_hor,                                 o_init_zero=do_zero_flux)
        case default
            error stop 'do_oce_adv_tra: unknown tra_adv_hor (expected UPW1/MUSCL/MFCT)'
        end select

        if (trim(tracers%data(tr_num)%tra_adv_lim) == 'FCT') then
            pwvel => w
        else
            pwvel => we
        end if

        !_______________________________________________________________________
        ! vertical advection (FCT: high-order antidiffusive; else: full high-order)
        select case (trim(tracers%data(tr_num)%tra_adv_ver))
        case ('QR4C')
            call adv_tra_ver_qr4c(pwvel, ttfAB, mesh, optv, adv_flux_ver, o_init_zero=do_zero_flux)
        case ('UPW1')
            call adv_tra_ver_upw1(pwvel, ttfAB, mesh,       adv_flux_ver, o_init_zero=do_zero_flux)
        case ('CDIFF', 'PPM')
            error stop 'do_oce_adv_tra: tra_adv_ver CDIFF/PPM not yet ported (M2)'
        case default
            error stop 'do_oce_adv_tra: unknown tra_adv_ver (expected UPW1/QR4C)'
        end select

        !_______________________________________________________________________
        ! limit (FCT) and scatter to the per-node tendencies del_ttf_advhoriz/advvert
        if (trim(tracers%data(tr_num)%tra_adv_lim) == 'FCT') then
            call oce_tra_adv_fct(dt, ttf, fct_LO, adv_flux_hor, adv_flux_ver, &
                                 fct_ttf_min, fct_ttf_max, fct_plus, fct_minus, mesh)
            call oce_tra_adv_flux2dtracer(dt, dttf_h, dttf_v, adv_flux_hor, adv_flux_ver, mesh, &
                                          use_lo=.true., ttf=ttf, lo=fct_LO)
        else
            call oce_tra_adv_flux2dtracer(dt, dttf_h, dttf_v, adv_flux_hor, adv_flux_ver, mesh)
        end if
    end subroutine do_oce_adv_tra

end module oce_adv_tra_driver
