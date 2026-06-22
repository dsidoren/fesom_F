module oce_ale_tracer
    ! Tracer step integration, transcribed from FESOM2 v2.7.3 oce_ale_tracer.F90.
    !
    ! M1.4 (advection): adv_tracers_ale loops over tracers; per tracer it runs the
    ! advection body (advect_tracer):
    !   init_tracers_AB  -> zero del_ttf*, AB-interpolate valuesAB, rebuild edge gradient
    !   do_oce_adv_tra   -> del_ttf_advhoriz / del_ttf_advvert
    !   del_ttf += del_ttf_advhoriz + del_ttf_advvert
    !
    ! M2.9a (diffusion solve): diff_tracers_ale is the per-tracer diffusion + ALE
    ! reconstruct, transcribed from FESOM2 oce_ale_tracer.F90:335-491:
    !   diff_part_hor_redi       -> del_ttf += horizontal diffusion  (Redi=.false. branch)
    !   ALE reconstruct          -> del_ttf += T*(hnode-hnode_new) ; T += del_ttf/hnode_new
    !   diff_ver_part_impl_ale   -> implicit vertical diffusion (Thomas/TDMA), the FIRST
    !                               consumer of the M2.8 PP vertical diffusivity dyn%work%Kv
    ! Gated for the pi / reduced-M2 config:
    !   Redi=.false.        -> the isoneutral slope (slope_tapered/Ki) terms vanish (isredi=0)
    !   tra_adv_lim='FCT'   -> do_wimpl=.false. (implicit vertical advection off; the FCT
    !                          scheme carries vertical advection explicitly via the full w)
    !   i_vert_diff=.true.  -> the implicit (not explicit) vertical-diffusion path
    !   PP mixing (not KPP) -> use_kpp_nonlclflx=.false.; no sw-penetration / icebergs
    ! The branches deferred for that config (Redi explicit slopes, KPP nonlocal fluxes,
    ! shortwave penetration, icebergs, the explicit diff_ver_part_expl_ale, the biharmonic
    ! diff_part_bh) enter with their enabling features (M2.10+); the do_wimpl advection
    ! terms are transcribed but unexercised on FCT (gated .false.). 1-rank only; the FESOM2
    ! exchange_nod(values) loop tail is a no-op here (M2.12).
    use mod_precision,      only: WP, MP
    use mod_constants,      only: vcpw
    use mod_mesh,           only: t_mesh
    use mod_dyn,            only: t_dyn
    use mod_tracer,         only: t_tracer
    use mod_partit,         only: t_partit
    use mod_part_bounds,    only: owned_bounds, is_multirank, local_dims
    use mod_halo,           only: exchange_nod
    use oce_tracer_mod,     only: init_tracers_AB
    use oce_tracer_grad,    only: tracer_gradient_elements
    use oce_adv_tra_driver, only: do_oce_adv_tra
    implicit none
    private
    public :: adv_tracers_ale, advect_tracer, diff_tracers_ale, diff_ver_part_impl_ale
    ! M2.9b (step assembly): solve_tracers_ale is the full per-tracer wrapper
    ! (advection + diffusion solve + salinity clamp) the ocean step calls.
    public :: solve_tracers_ale

contains

    subroutine adv_tracers_ale(dt, dynamics, tracers, mesh)
        real(kind=WP),  intent(in)            :: dt
        type(t_mesh),   intent(in)            :: mesh
        type(t_dyn),    intent(inout), target :: dynamics
        type(t_tracer), intent(inout)         :: tracers
        integer :: tr_num
        do tr_num = 1, tracers%num_tracers
            call advect_tracer(dt, tr_num, dynamics, tracers, mesh)
            ! M2.9b (step assembly): diff_tracers_ale (this module) + decay + relax_to_clim
            !     + the salinity clamp + exchange_nod(values) wrap into solve_tracers_ale.
        end do
    end subroutine adv_tracers_ale

    subroutine advect_tracer(dt, tr_num, dynamics, tracers, mesh, partit)
        ! One tracer's advection contribution to del_ttf (FESOM2 adv_tracers_ale body).
        ! M2.12b: optional partit threaded to init_tracers_AB / do_oce_adv_tra; the
        ! del_ttf accumulation runs over OWNED nodes (the gated tendency).
        real(kind=WP),  intent(in)            :: dt
        integer,        intent(in)            :: tr_num
        type(t_mesh),   intent(in)            :: mesh
        type(t_dyn),    intent(inout), target :: dynamics
        type(t_tracer), intent(inout)         :: tracers
        type(t_partit), intent(in), optional  :: partit
        integer :: n
        integer :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        call init_tracers_AB(tr_num, tracers, mesh, partit)
        call do_oce_adv_tra(dt, dynamics%uv, dynamics%w, dynamics%w_i, dynamics%w_e, &
                            tr_num, dynamics, tracers, mesh, partit)
        ! total tracer tendency = horizontal + vertical advection (del_ttf was zeroed
        ! in init_tracers_AB).
        do n = 1, nNodO
            tracers%work%del_ttf(:, n) = tracers%work%del_ttf(:, n) &
                                       + tracers%work%del_ttf_advhoriz(:, n) &
                                       + tracers%work%del_ttf_advvert(:, n)
        end do
    end subroutine advect_tracer

    !===========================================================================
    subroutine solve_tracers_ale(dt, dynamics, tracers, mesh, Ki, &
                                 heat_flux, water_flux, virtual_salt, relax_salt, &
                                 real_salt_flux, is_nonlinfs, partit)
        ! FESOM2 oce_ale_tracer.F90:135-331, the pi / reduced-M2 path: the full per-tracer
        ! solve the ocean step (M2.9b) calls. Per tracer, in order:
        !   advect_tracer    -> init_tracers_AB (zeros del_ttf, AB-interpolate valuesAB,
        !                       rebuild edge gradient) + do_oce_adv_tra + del_ttf +=
        !                       del_ttf_advhoriz + del_ttf_advvert
        !   tr_xy            -> elemental gradient of the still-un-updated tracer (T^n),
        !                       computed here and passed in (FESOM2 builds it inside
        !                       diff_tracers_ale; same value, same timing)
        !   diff_tracers_ale -> horizontal diffusion + ALE reconstruct + implicit vertical
        !                       diffusion TDMA (consumes the LIVE dyn%work%Kv from M2.8)
        !   relax_to_clim    -> clim_relax=0 on pi -> short-circuits (3D restoring; M2.10)
        !   exchange_nod     -> 1-rank no-op (lifted at M2.12)
        ! then the salinity clamp S in [3, 45] over all nodes.
        !
        ! Gated OFF by the reduced-M2 namelist (transcribed-deferred): SPP rejected-salt,
        ! Fer_GM bolus add/subtract, the toy relaxations (soufflet/neverworld2/dbgyre),
        ! radioactive decay (14C/39Ar), the ptracers 3D restore, the age-tracer clamp.
        real(kind=WP),  intent(in)            :: dt
        type(t_dyn),    intent(inout), target :: dynamics
        type(t_tracer), intent(inout), target :: tracers
        type(t_mesh),   intent(in),    target :: mesh
        real(kind=WP),  intent(in)            :: Ki(mesh%nl-1, mesh%nod2D)
        real(kind=WP),  intent(in)            :: heat_flux(mesh%nod2D), water_flux(mesh%nod2D)
        real(kind=WP),  intent(in)            :: virtual_salt(mesh%nod2D), relax_salt(mesh%nod2D)
        real(kind=WP),  intent(in)            :: real_salt_flux(mesh%nod2D), is_nonlinfs
        type(t_partit), intent(in), optional  :: partit
        integer :: tr_num, node, nzmin, nzmax
        integer :: nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF
        real(kind=WP), allocatable :: tr_xy(:,:,:)
        real(kind=WP), dimension(:,:), pointer :: Svalues

        ! M2.12c-3: tr_xy is sized to the LOCAL element count (mesh%elem2D holds the GLOBAL
        ! count in the partitioned mesh). tracer_gradient_elements writes OWNED elements;
        ! diff_part_hor_redi reads tr_xy only at the (owned) triangles of owned edges.
        call local_dims(mesh, partit, nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF)
        allocate(tr_xy(2, mesh%nl-1, nElemF))

        do tr_num = 1, tracers%num_tracers
            ! advection: del_ttf = advhoriz + advvert (del_ttf zeroed in init_tracers_AB)
            call advect_tracer(dt, tr_num, dynamics, tracers, mesh, partit)
            ! elemental gradient of the pre-diffusion tracer (advection left values = T^n)
            call tracer_gradient_elements(tracers%data(tr_num)%values, tr_xy, mesh, partit)
            ! horizontal diffusion + ALE reconstruct + implicit vertical-diffusion TDMA
            call diff_tracers_ale(tr_num, dt, dynamics, tracers, mesh, tr_xy, Ki, &
                                  heat_flux, water_flux, virtual_salt, relax_salt, &
                                  real_salt_flux, is_nonlinfs, partit)
            ! relax_to_clim (clim_relax=0): no-op. exchange_nod(values) (FESOM2 :268): the
            ! owned values are complete (diff over owned nodes/edges); this fills the halo
            ! for the NEXT step's init_tracers_AB / advection (and for the salinity clamp's
            ! owned+halo loop below).
            if (is_multirank(partit)) call exchange_nod(tracers%data(tr_num)%values, partit)
        end do

        ! salinity clamp (tracer 2 = salinity, FESOM2 :304-316): S in [3, 45], owned+halo
        ! (the clamp is per-node idempotent, so the owned values match the 1-rank result).
        Svalues => tracers%data(2)%values
        do node = 1, nNodL
            nzmax = mesh%nlevels_nod2D(node) - 1
            nzmin = mesh%ulevels_nod2D(node)
            where (Svalues(nzmin:nzmax,node) > 45.0_WP) Svalues(nzmin:nzmax,node) = 45.0_WP
            where (Svalues(nzmin:nzmax,node) <  3.0_WP) Svalues(nzmin:nzmax,node) =  3.0_WP
        end do

        deallocate(tr_xy)
    end subroutine solve_tracers_ale

    !===========================================================================
    subroutine diff_tracers_ale(tr_num, dt, dynamics, tracers, mesh, tr_xy, Ki, &
                                heat_flux, water_flux, virtual_salt, relax_salt, &
                                real_salt_flux, is_nonlinfs, partit)
        ! Per-tracer diffusion + ALE tracer reconstruct (FESOM2 oce_ale_tracer.F90:335-491,
        ! Redi=.false. / i_vert_diff=.true. / PP path). del_ttf enters with the advection
        ! tendency; this routine adds horizontal diffusion, reconstructs the new tracer
        ! T* = (h^n*T^n + del_ttf)/h^{n+1}, then applies the implicit vertical diffusion.
        !
        ! tr_xy / Ki / the surface flux arrays are explicit arguments here (D7; they are
        ! external inputs not yet sourced in FESOM3: Ki needs mesh_resolution [M4], the
        ! surface fluxes come from forcing [M2.10] — the M2.5 prescribe-the-input precedent).
        integer,        intent(in)            :: tr_num
        real(kind=WP),  intent(in)            :: dt
        type(t_dyn),    intent(inout), target :: dynamics
        type(t_tracer), intent(inout), target :: tracers
        type(t_mesh),   intent(in),    target :: mesh
        real(kind=WP),  intent(in)            :: tr_xy(2, mesh%nl-1, mesh%elem2D)
        real(kind=WP),  intent(in)            :: Ki(mesh%nl-1, mesh%nod2D)
        real(kind=WP),  intent(in)            :: heat_flux(mesh%nod2D), water_flux(mesh%nod2D)
        real(kind=WP),  intent(in)            :: virtual_salt(mesh%nod2D), relax_salt(mesh%nod2D)
        real(kind=WP),  intent(in)            :: real_salt_flux(mesh%nod2D), is_nonlinfs
        type(t_partit), intent(in), optional  :: partit
        integer :: n, nzmin, nzmax
        integer :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP), dimension(:,:), pointer :: trarr
        real(kind=MP), dimension(:,:), pointer :: del_ttf

        trarr   => tracers%data(tr_num)%values
        del_ttf => tracers%work%del_ttf
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        ! horizontal diffusion: del_ttf += R_T^n (Redi=.false. -> plain Laplacian)
        call diff_part_hor_redi(tr_xy, Ki, dt, tracers, mesh, partit)

        ! explicit vertical diffusion (diff_ver_part_expl_ale) is skipped: i_vert_diff=.true.
        ! Redi vertical projection (diff_ver_part_redi_expl) is skipped: Redi=.false.

        !_______________________________________________________________________
        ! ALE tracer reconstruct: T* = (dt*R_T^n + h^{n-0.5}*T^{n-0.5})/h^{n+0.5}
        ! (FESOM2 :464-477). For linfs hnode_new==hnode so the (hnode-hnode_new) term
        ! vanishes and this collapses to T* = T + del_ttf/hnode.
        ! M2.12c-3: OWNED node loop (FESOM2 :465 do n=1, myDim_nod2D); del_ttf at owned nodes
        ! is complete (advection over owned nodes + diff_part_hor_redi over owned edges, both
        ! invariant-i complete), so T* is correct at owned nodes — the halo is filled by
        ! solve_tracers_ale's exchange_nod(values).
        do n = 1, nNodO
            nzmax = mesh%nlevels_nod2D(n) - 1
            nzmin = mesh%ulevels_nod2D(n)
            del_ttf(nzmin:nzmax,n) = del_ttf(nzmin:nzmax,n) + trarr(nzmin:nzmax,n)* &
                                     (mesh%hnode(nzmin:nzmax,n)-mesh%hnode_new(nzmin:nzmax,n))
            trarr(nzmin:nzmax,n)   = trarr(nzmin:nzmax,n) + &
                                     del_ttf(nzmin:nzmax,n)/mesh%hnode_new(nzmin:nzmax,n)
        end do

        !_______________________________________________________________________
        ! implicit vertical diffusion (i_vert_diff=.true.): the TDMA, consumes dyn%work%Kv
        call diff_ver_part_impl_ale(tr_num, dt, dynamics, tracers, mesh, &
                                    heat_flux, water_flux, virtual_salt, relax_salt, &
                                    real_salt_flux, is_nonlinfs, partit)

        ! biharmonic tracer diffusion (diff_part_bh) is skipped: smooth_bh_tra=.false.
    end subroutine diff_tracers_ale

    !===========================================================================
    subroutine diff_part_hor_redi(tr_xy, Ki, dt, tracers, mesh, partit)
        ! Horizontal tracer diffusion, Redi=.false. branch of FESOM2 oce_ale_tracer.F90:1173.
        ! Edge-based flux form: across each edge the diffusive flux Kh*(Tx,Ty) (Kh = mean of
        ! the two edge-node diffusivities Ki, (Tx,Ty) = elemental tracer gradient tr_xy)
        ! crossing the edge mid-faces (edge_cross_dxdy) is scattered, with opposite sign, into
        ! del_ttf at the two edge nodes. The Redi isoneutral-slope terms (slope_tapered/tr_z,
        ! x isredi) vanish for Redi=off; they are deferred to a future Redi gate.
        ! M2.12c-3: optional partit -> OWNED edge loop (FESOM2 :1208 do edge=1, myDim_edge2D).
        ! For owned edges both triangles are owned (invariant ii) so tr_xy/helem read in the
        ! owned-only element range; Ki/areasvol are read at the edge's two nodes (one may be a
        ! halo node — Ki is prescribed owned+halo, areasvol is M2.12b-exchanged). The owned-
        ! node del_ttf scatter is complete (invariant i).
        type(t_tracer), intent(inout), target :: tracers
        type(t_mesh),   intent(in),    target :: mesh
        real(kind=WP),  intent(in)            :: tr_xy(2, mesh%nl-1, mesh%elem2D)
        real(kind=WP),  intent(in)            :: Ki(mesh%nl-1, mesh%nod2D)
        real(kind=WP),  intent(in)            :: dt
        type(t_partit), intent(in), optional  :: partit
        integer :: edge, nz, el(2), enodes(2)
        integer :: nl1, ul1, nl2, ul2, nl12, ul12
        integer :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: deltaX1, deltaY1, deltaX2, deltaY2, c, Fx, Fy, Tx, Ty, Kh, dz
        real(kind=WP) :: rhs1(mesh%nl-1), rhs2(mesh%nl-1)
        real(kind=MP), dimension(:,:), pointer :: del_ttf

        del_ttf => tracers%work%del_ttf
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        do edge = 1, nEdgeO
            rhs1 = 0.0_WP
            rhs2 = 0.0_WP
            deltaX1 = mesh%edge_cross_dxdy(1,edge)
            deltaY1 = mesh%edge_cross_dxdy(2,edge)
            el      = mesh%edge_tri(:,edge)
            enodes  = mesh%edges(:,edge)
            nl1     = mesh%nlevels(el(1)) - 1
            ul1     = mesh%ulevels(el(1))
            nl2 = 0
            ul2 = 0
            if (el(2) > 0) then
                nl2 = mesh%nlevels(el(2)) - 1
                ul2 = mesh%ulevels(el(2))
                deltaX2 = mesh%edge_cross_dxdy(3,edge)
                deltaY2 = mesh%edge_cross_dxdy(4,edge)
            end if
            nl12 = min(nl1, nl2)
            ul12 = max(ul1, ul2)
            ! (A) levels present only in el(1) above the shared range
            do nz = ul1, ul12-1
                Kh = sum(Ki(nz, enodes))/2.0_WP
                dz = mesh%helem(nz, el(1))
                Tx = tr_xy(1,nz,el(1)); Ty = tr_xy(2,nz,el(1))
                Fx = Kh*Tx; Fy = Kh*Ty
                c  = (-deltaX1*Fy + deltaY1*Fx)*dz
                rhs1(nz) = rhs1(nz) + c
                rhs2(nz) = rhs2(nz) - c
            end do
            ! (B) levels present only in el(2) above the shared range
            if (ul2 > 0) then
                do nz = ul2, ul12-1
                    Kh = sum(Ki(nz, enodes))/2.0_WP
                    dz = mesh%helem(nz, el(2))
                    Tx = tr_xy(1,nz,el(2)); Ty = tr_xy(2,nz,el(2))
                    Fx = Kh*Tx; Fy = Kh*Ty
                    c  = (deltaX2*Fy - deltaY2*Fx)*dz
                    rhs1(nz) = rhs1(nz) + c
                    rhs2(nz) = rhs2(nz) - c
                end do
            end if
            ! (C) shared range: both elements contribute (averaged gradient + full cross)
            do nz = ul12, nl12
                Kh = sum(Ki(nz, enodes))/2.0_WP
                dz = sum(mesh%helem(nz, el))/2.0_WP
                Tx = 0.5_WP*(tr_xy(1,nz,el(1))+tr_xy(1,nz,el(2)))
                Ty = 0.5_WP*(tr_xy(2,nz,el(1))+tr_xy(2,nz,el(2)))
                Fx = Kh*Tx; Fy = Kh*Ty
                c  = ((deltaX2-deltaX1)*Fy - (deltaY2-deltaY1)*Fx)*dz
                rhs1(nz) = rhs1(nz) + c
                rhs2(nz) = rhs2(nz) - c
            end do
            ! (D) deeper levels present only in el(1)
            do nz = nl12+1, nl1
                Kh = sum(Ki(nz, enodes))/2.0_WP
                dz = mesh%helem(nz, el(1))
                Tx = tr_xy(1,nz,el(1)); Ty = tr_xy(2,nz,el(1))
                Fx = Kh*Tx; Fy = Kh*Ty
                c  = (-deltaX1*Fy + deltaY1*Fx)*dz
                rhs1(nz) = rhs1(nz) + c
                rhs2(nz) = rhs2(nz) - c
            end do
            ! (E) deeper levels present only in el(2)
            do nz = nl12+1, nl2
                Kh = sum(Ki(nz, enodes))/2.0_WP
                dz = mesh%helem(nz, el(2))
                Tx = tr_xy(1,nz,el(2)); Ty = tr_xy(2,nz,el(2))
                Fx = Kh*Tx; Fy = Kh*Ty
                c  = (deltaX2*Fy - deltaY2*Fx)*dz
                rhs1(nz) = rhs1(nz) + c
                rhs2(nz) = rhs2(nz) - c
            end do
            ! scatter into del_ttf at both edge nodes over their full joint level range
            nl12 = max(nl1, nl2)
            ul12 = ul1
            if (ul2 > 0) ul12 = min(ul1, ul2)
            del_ttf(ul12:nl12,enodes(1)) = del_ttf(ul12:nl12,enodes(1)) &
                + rhs1(ul12:nl12)*dt/mesh%areasvol(ul12:nl12,enodes(1))
            del_ttf(ul12:nl12,enodes(2)) = del_ttf(ul12:nl12,enodes(2)) &
                + rhs2(ul12:nl12)*dt/mesh%areasvol(ul12:nl12,enodes(2))
        end do
    end subroutine diff_part_hor_redi

    !===========================================================================
    subroutine diff_ver_part_impl_ale(tr_num, dt, dynamics, tracers, mesh, &
                                      heat_flux, water_flux, virtual_salt, relax_salt, &
                                      real_salt_flux, is_nonlinfs, partit)
        ! Implicit vertical tracer diffusion (FESOM2 oce_ale_tracer.F90:562-1082), the
        ! Redi=.false. / PP path. Per node it builds the tridiagonal system for the implicit
        ! vertical-diffusion increment dTnew = T^{n+0.5} - T*, with the vertical diffusivity
        ! K_33 = Kv (the M2.8 PP/convective coefficient, here dyn%work%Kv — the FIRST consumer
        ! of a PP output), the surface boundary flux from bc_surface (heat / virtual-salt /
        ! relaxation), and solves it with the Thomas algorithm; then T += dTnew.
        !
        ! Redi (isredi=0) zeroes the isoneutral Ty/Ty1 terms; do_wimpl (use_wsplit .and. not
        ! FCT) adds the implicit vertical advection, transcribed but .false. on the pi FCT
        ! tracers. KPP nonlocal fluxes / shortwave penetration / icebergs are deferred with
        ! their features (M2.10+).
        integer,        intent(in)            :: tr_num
        real(kind=WP),  intent(in)            :: dt
        type(t_dyn),    intent(inout), target :: dynamics
        type(t_tracer), intent(inout), target :: tracers
        type(t_mesh),   intent(in),    target :: mesh
        real(kind=WP),  intent(in)            :: heat_flux(mesh%nod2D), water_flux(mesh%nod2D)
        real(kind=WP),  intent(in)            :: virtual_salt(mesh%nod2D), relax_salt(mesh%nod2D)
        real(kind=WP),  intent(in)            :: real_salt_flux(mesh%nod2D), is_nonlinfs
        type(t_partit), intent(in), optional  :: partit
        !
        real(kind=WP) :: a(mesh%nl), b(mesh%nl), c(mesh%nl), tr(mesh%nl)
        real(kind=WP) :: cp(mesh%nl), tp(mesh%nl)
        real(kind=WP) :: zbar_n(mesh%nl), Z_n(mesh%nl-1)
        integer       :: nz, n, nzmax, nzmin, id
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: m, zinv, dz, zinv1, zinv2, v_adv
        logical       :: do_wimpl
        real(kind=WP), dimension(:,:), pointer :: trarr, Wvel_i

        trarr  => tracers%data(tr_num)%values
        Wvel_i => dynamics%w_i
        id     =  tracers%data(tr_num)%ID
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        ! FCT tracers (or no wsplit) -> implicit vertical advection off
        do_wimpl = .true.
        if ((trim(tracers%data(tr_num)%tra_adv_lim) == 'FCT') .or. (.not. dynamics%use_wsplit)) &
            do_wimpl = .false.

        ! M2.12c-3: OWNED node loop (FESOM2 :726 do n=1, myDim_nod2D). The TDMA is per-column
        ! (no halo coupling); each owned column reads its own Kv/hnode_new/area/areasvol + the
        ! prescribed surface fluxes at n, and updates trarr at owned nodes. The halo values
        ! are filled by solve_tracers_ale's exchange_nod(values).
        do n = 1, nNodO
            a  = 0.0_WP; b = 0.0_WP; c = 0.0_WP; tr = 0.0_WP; tp = 0.0_WP; cp = 0.0_WP
            nzmax = mesh%nlevels_nod2D(n)
            nzmin = mesh%ulevels_nod2D(n)
            !___________________________________________________________________
            ! rebuild zbar_n (layer interfaces) / Z_n (mid-depths) from the NEW
            ! thickness hnode_new, bottom-up (FESOM2 :743-751). zbar_n_bot(n) is the
            ! node bottom depth = zbar_3d_n(nzmax,n) (full cells).
            zbar_n = 0.0_WP; Z_n = 0.0_WP
            zbar_n(nzmax)   = mesh%zbar_3d_n(nzmax, n)
            Z_n(nzmax-1)    = zbar_n(nzmax) + mesh%hnode_new(nzmax-1,n)/2.0_WP
            do nz = nzmax-1, nzmin+1, -1
                zbar_n(nz)  = zbar_n(nz+1) + mesh%hnode_new(nz,n)
                Z_n(nz-1)   = zbar_n(nz)   + mesh%hnode_new(nz-1,n)/2.0_WP
            end do
            zbar_n(nzmin)   = zbar_n(nzmin+1) + mesh%hnode_new(nzmin,n)
            !___________________________________________________________________
            ! surface layer coefficients (isredi=0 -> no Ty1 term)
            nz = nzmin
            zinv2 = 1.0_WP/(Z_n(nz)-Z_n(nz+1))
            zinv  = 1.0_WP*dt
            a(nz) = 0.0_WP
            c(nz) = -dynamics%work%Kv(nz+1,n)*zinv2*zinv * mesh%area(nz+1,n)/mesh%areasvol(nz,n)
            b(nz) = -c(nz) + mesh%hnode_new(nz,n)
            if (do_wimpl) then
                v_adv = zinv * ( mesh%area(nz  ,n)/mesh%areasvol(nz,n) )
                b(nz) = b(nz) + Wvel_i(nz, n)*v_adv
                v_adv = zinv * mesh%area(nz+1,n)/mesh%areasvol(nz,n)
                b(nz) = b(nz) - min(0._WP, Wvel_i(nz+1, n))*v_adv
                c(nz) = c(nz) - max(0._WP, Wvel_i(nz+1, n))*v_adv
            end if
            zinv1 = zinv2   ! backup 1/dz for the next (deeper) level, as FESOM2
            !___________________________________________________________________
            ! interior layers
            do nz = nzmin+1, nzmax-2
                zinv2 = 1.0_WP/(Z_n(nz)-Z_n(nz+1))
                a(nz) = -dynamics%work%Kv(nz,n)  *zinv1*zinv * ( mesh%area(nz  ,n)/mesh%areasvol(nz,n) )
                c(nz) = -dynamics%work%Kv(nz+1,n)*zinv2*zinv *   mesh%area(nz+1,n)/mesh%areasvol(nz,n)
                b(nz) = -a(nz)-c(nz) + mesh%hnode_new(nz,n)
                zinv1 = zinv2
                if (do_wimpl) then
                    v_adv = zinv * ( mesh%area(nz  ,n)/mesh%areasvol(nz,n) )
                    a(nz) = a(nz) + min(0._WP, Wvel_i(nz, n))*v_adv
                    b(nz) = b(nz) + max(0._WP, Wvel_i(nz, n))*v_adv
                    v_adv = zinv * mesh%area(nz+1,n)/mesh%areasvol(nz,n)
                    b(nz) = b(nz) - min(0._WP, Wvel_i(nz+1, n))*v_adv
                    c(nz) = c(nz) - max(0._WP, Wvel_i(nz+1, n))*v_adv
                end if
            end do
            !___________________________________________________________________
            ! bottom layer (nz = nzmax-1)
            nz = nzmax-1
            zinv = 1.0_WP*dt
            a(nz) = -dynamics%work%Kv(nz,n)*zinv1*zinv * ( mesh%area(nz  ,n)/mesh%areasvol(nz,n) )
            c(nz) = 0.0_WP
            b(nz) = -a(nz) + mesh%hnode_new(nz,n)
            if (do_wimpl) then
                v_adv = zinv * ( mesh%area(nz  ,n)/mesh%areasvol(nz,n) )
                a(nz) = a(nz) + min(0._WP, Wvel_i(nz, n))*v_adv
                b(nz) = b(nz) + max(0._WP, Wvel_i(nz, n))*v_adv
            end if
            !___________________________________________________________________
            ! rhs = K_33*dt*d/dz*T* in flux form (= -a*T*_{i-1} -(b-dz)*T*_i -c*T*_{i+1})
            nz = nzmin
            dz = mesh%hnode_new(nz,n)
            tr(nz) = -(b(nz)-dz)*trarr(nz,n) - c(nz)*trarr(nz+1,n)
            do nz = nzmin+1, nzmax-2
                dz = mesh%hnode_new(nz,n)
                tr(nz) = -a(nz)*trarr(nz-1,n) - (b(nz)-dz)*trarr(nz,n) - c(nz)*trarr(nz+1,n)
            end do
            nz = nzmax-1
            dz = mesh%hnode_new(nz,n)
            tr(nz) = -a(nz)*trarr(nz-1,n) - (b(nz)-dz)*trarr(nz,n)
            ! (KPP nonlocal / shortwave penetration / iceberg rhs contributions: deferred)
            !___________________________________________________________________
            ! surface boundary flux (heat / virtual-salt / relaxation)
            tr(nzmin) = tr(nzmin) + bc_surface(id, trarr(nzmin,n), dt, heat_flux(n), &
                            water_flux(n), virtual_salt(n), relax_salt(n), &
                            real_salt_flux(n), is_nonlinfs)
            !___________________________________________________________________
            ! Thomas algorithm: forward elimination then back substitution
            cp(nzmin) = c(nzmin)/b(nzmin)
            tp(nzmin) = tr(nzmin)/b(nzmin)
            do nz = nzmin+1, nzmax-1
                m      = b(nz) - cp(nz-1)*a(nz)
                cp(nz) = c(nz)/m
                tp(nz) = (tr(nz)-tp(nz-1)*a(nz))/m
            end do
            tr(nzmax-1) = tp(nzmax-1)
            do nz = nzmax-2, nzmin, -1
                tr(nz) = tp(nz) - cp(nz)*tr(nz+1)
            end do
            !___________________________________________________________________
            ! update tracer: trarr (=T*) + dTnew
            do nz = nzmin, nzmax-1
                trarr(nz,n) = trarr(nz,n) + tr(nz)
            end do
        end do
    end subroutine diff_ver_part_impl_ale

    !===========================================================================
    real(kind=WP) function bc_surface(id, sval, dt, heat_flux, water_flux, &
                                      virtual_salt, relax_salt, real_salt_flux, is_nonlinfs)
        ! Surface boundary flux added to the top tridiagonal row (FESOM2 oce_ale_tracer.F90
        ! :1475-1680). Only temperature (id=1) and salinity (id=2) are reached on pi; the
        ! transient-tracer / recom / wiso / age cases are deferred with their tracers.
        integer,       intent(in) :: id
        real(kind=WP), intent(in) :: sval, dt, heat_flux, water_flux
        real(kind=WP), intent(in) :: virtual_salt, relax_salt, real_salt_flux, is_nonlinfs
        select case (id)
        case (1)  ! temperature
            bc_surface = -dt*(heat_flux/vcpw + sval*water_flux*is_nonlinfs)
        case (2)  ! salinity (virtual_salt is 0 for zlevel/zstar; nonzero for linfs)
            bc_surface =  dt*(virtual_salt + relax_salt + real_salt_flux*is_nonlinfs)
        case default
            bc_surface = 0.0_WP
        end select
    end function bc_surface

end module oce_ale_tracer
