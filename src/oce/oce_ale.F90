module oce_ale
    ! M2.7 — ALE (linfs) velocity / SSH / thickness-W update: the post-CG tail of the
    ! FESOM2 timestep (oce_ale.F90:3946-4081 + the inline eta_n blend at :3973-3977).
    !
    ! After the free-surface CG solve (M2.6 -> d_eta) the timestep does, in order:
    !   1. update_vel        — u^(n+1) = u* + UV_rhs + [-g*theta*dt*grad(d_eta)]
    !                          (FESOM2 oce_dyn.F90:88-173)
    !   2. compute_hbar_ale  — hbar_old=hbar; hbar += dt/areasvol * div(UV); dhe
    !                          (FESOM2 oce_ale.F90:2165-2307)
    !   3. update_eta_n      — eta_n = alpha*hbar + (1-alpha)*hbar_old (inline :3973-3977)
    !   4. vert_vel_ale      — W = -cumsum(div(UV*h))/area; linfs leaves hnode_new=hnode;
    !                          then compute_CFLz + compute_Wvel_split
    !                          (FESOM2 oce_ale.F90:2323-2884, the which_ale='linfs' path)
    !
    ! Scope reductions (all the gated-config branches the pi/reduced-M2 namelist takes):
    !  - which_ale='linfs' -> the zlevel/zstar thickness redistribution in vert_vel_ale is
    !    NOT taken, so hnode_new stays = hnode (its init value) and there is no surface
    !    Wvel/hnode correction; the compute_hbar_ale water_flux term ((1-alpha) / non-linfs)
    !    vanishes. The full-free-surface thickness evolution gets its own later gate.
    !  - Fer_GM=.false. -> the fer_UV/fer_Wvel transport in vert_vel_ale is dropped.
    !  - ldiag_ke=.false. -> the update_vel ke_* energy diagnostics are dropped (no ke_*
    !    in t_dyn, as in compute_vel_rhs / impl_vert_visc_ale).
    !  - 1-rank: myDim_nod2D+eDim_nod2D == myDim_nod2D == nod2D, myDim_edge2D == edge2D;
    !    every exchange_nod/exchange_elem is a no-op (lifted at M2.12). The OpenMP-off
    !    oracle (L16/L17/L19) compiles the omp locks/ORDERED out -> serial edge order, which
    !    matches these serial loops, so the divergence/W scatters byte-match by L9.
    !
    ! Byte-gate (M2.7, tools/run_pressure_gate.sh): every operand is already pinned —
    ! d_eta (M2.6), UV/UV_rhs (the post-TDMA M2.5 value), gradient_sca/edge_cross_dxdy/
    ! helem/areasvol/area/levels (geometry + M2.5), g/density_0/dt/alpha/theta — so the
    ! whole chain is max|delta|=0 by transitivity (L9). hbar is prescribed identically.
    use mod_precision,  only: WP, MP
    use mod_constants,  only: g
    use mod_param_phys, only: alpha, theta
    use mod_mesh,       only: t_mesh
    use mod_dyn,        only: t_dyn
    use mod_partit,     only: t_partit
    use mod_part_bounds, only: owned_bounds, is_multirank
    use mod_halo,       only: exchange_nod, exchange_elem, exchange_elem_full
    implicit none
    private
    public :: update_vel, compute_hbar_ale, update_eta_n, vert_vel_ale
    ! M2.9b (step assembly): compute_vel_nodes is the REAL dyn%uvnode source (prescribed
    ! at M2.8); update_thickness_ale is the linfs thickness commit (a no-op for which_ale=
    ! 'linfs'). Both join the faithful oce_timestep sequence assembled in mod_step_oce.
    public :: compute_vel_nodes, update_thickness_ale

contains

    !===========================================================================
    subroutine compute_vel_nodes(dynamics, mesh, partit)
        ! FESOM2 oce_dyn.F90:177-226. Horizontal velocity averaged from elements to
        ! nodes (area-weighted over the node's element star) -> dyn%uvnode, the input
        ! PP mixing (M2.8 oce_mixing_pp) reads for the vertical shear. Called at the top
        ! of the ocean step (FESOM2 fesom_module.F90:654, before oce_timestep_ale).
        ! M2.12c: optional partit -> owned node loop (FESOM2 :202 do n=1,myDim_nod2D) +
        ! the trailing exchange_nod(UVnode). The per-node accumulation over nod_in_elem2D
        ! is order-sensitive but byte-pinned by the geometry gate (L9), so it matches
        ! FESOM2's local element order on the same partition.
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(in),    target :: mesh
        type(t_partit), intent(in), optional :: partit
        !______________________________________________________________________
        integer       :: n, nz, k, elem, nln, uln, nle, ule
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: tx, ty, tvol
        real(kind=WP), dimension(:,:,:), pointer :: UV, UVnode

        UV     => dynamics%uv
        UVnode => dynamics%uvnode
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        do n = 1, nNodO
            uln = mesh%ulevels_nod2D(n)
            nln = mesh%nlevels_nod2D(n)
            do nz = uln, nln-1
                tvol = 0.0_WP; tx = 0.0_WP; ty = 0.0_WP
                do k = 1, mesh%nod_in_elem2D_num(n)
                    elem = mesh%nod_in_elem2D(k, n)
                    ule  = mesh%ulevels(elem)
                    nle  = mesh%nlevels(elem)
                    if (nle-1 < nz .or. nz < ule) cycle
                    tvol = tvol + mesh%elem_area(elem)
                    tx   = tx   + UV(1,nz,elem)*mesh%elem_area(elem)
                    ty   = ty   + UV(2,nz,elem)*mesh%elem_area(elem)
                end do
                UVnode(1,nz,n) = tx/tvol
                UVnode(2,nz,n) = ty/tvol
            end do
        end do
        if (is_multirank(partit)) call exchange_nod(UVnode, partit)   ! FESOM2 :225
    end subroutine compute_vel_nodes

    !===========================================================================
    subroutine update_vel(dynamics, mesh, dt, partit)
        ! FESOM2 oce_dyn.F90:88-173. New horizontal velocity after the SSH solve:
        !   u^(n+1) = u* + UV_rhs + [-g*theta*dt*grad(d_eta)]
        ! The SSH-gradient correction is a gradient_sca contraction of -g*theta*dt*d_eta
        ! over the element's three nodes (same shape as the M2.3 vel_rhs SSH-gradient and
        ! the M2.2 PGF). UV_rhs here is the post-TDMA value (M2.5; the CG solve does not
        ! touch it).
        ! M2.12c-3: optional partit -> OWNED element loop (FESOM2 :105 do elem=1,
        ! myDim_elem2D), then exchange_elem(UV) (:172). d_eta is read at the owned element's
        ! 3 nodes (within owned+halo, invariant iii) — valid from the CG's exchange_nod(x).
        ! The UV halo is refreshed owner->halo for the NEXT step (compute_vel_nodes/advection);
        ! the owned UV the post-update kernels read is computed directly here. exchange_elem
        ! is rank-1/2 only in mod_halo, so the rank-3 UV uses exchange_elem_full (eDim+eXDim,
        ! a superset of FESOM2's eDim exchange_elem; the owned values are identical and the
        ! owned-node dumps never depend on the eXDim halo).
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(in),    target :: mesh
        real(kind=WP), intent(in) :: dt
        type(t_partit), intent(in), optional :: partit
        !______________________________________________________________________
        integer       :: elem, nz, nzmin, nzmax, elnodes(3)
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: eta(3), Fx, Fy
        real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
        real(kind=WP), dimension(:),     pointer :: d_eta

        UV     => dynamics%uv
        UV_rhs => dynamics%uv_rhs
        d_eta  => dynamics%d_eta
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        do elem = 1, nElemO
            elnodes = mesh%elem2D_nodes(1:3, elem)        ! triangles (elem2D_nodes is MAX_NV=4) — L15
            eta = -g*theta*dt*d_eta(elnodes)
            Fx  = sum(mesh%gradient_sca(1:3, elem)*eta)
            Fy  = sum(mesh%gradient_sca(4:6, elem)*eta)
            nzmin = mesh%ulevels(elem)
            nzmax = mesh%nlevels(elem)
            do nz = nzmin, nzmax-1
                UV(1,nz,elem) = UV(1,nz,elem) + UV_rhs(1,nz,elem) + Fx
                UV(2,nz,elem) = UV(2,nz,elem) + UV_rhs(2,nz,elem) + Fy
            end do
        end do
        if (is_multirank(partit)) call exchange_elem_full(UV, partit)   ! FESOM2 :172
    end subroutine update_vel

    !===========================================================================
    subroutine compute_hbar_ale(dynamics, mesh, dt, partit)
        ! FESOM2 oce_ale.F90:2165-2307. Advance hbar (the elevation on semi-integer
        ! timesteps) by the depth-integrated horizontal divergence of the UPDATED UV:
        !   ssh_rhs_old = div_h( int(UV) dz )        (edge-scatter, NO alpha here)
        !   hbar_old    = hbar
        !   hbar        = hbar_old + ssh_rhs_old * dt / areasvol
        !   dhe         = sum(hbar - hbar_old)/3     (per element; feeds update_stiff,
        !                 which linfs never calls — gated anyway)
        ! Same edge-divergence structure as M2.6 compute_ssh_rhs_ale (same edge order,
        ! gated). linfs skips the water_flux term + its exchange_nod. Writes mesh%hbar/
        ! hbar_old/dhe and dynamics%ssh_rhs_old.
        ! M2.12c-3: optional partit -> FESOM2 :2197-2306 bounds: ssh_rhs_old zero owned+halo,
        ! the edge-divergence over OWNED edges (owned-node ssh_rhs_old complete, invariant i),
        ! hbar_old=hbar owned+halo, the hbar update OWNED, then exchange_nod(hbar) (:2293,
        ! ALWAYS — outside the linfs water_flux guard) so the dhe loop over OWNED elements can
        ! read hbar at the element's 3 nodes (within owned+halo). The linfs water_flux term +
        ! its exchange_nod(ssh_rhs_old) stay skipped.
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(inout), target :: mesh
        real(kind=WP), intent(in) :: dt
        type(t_partit), intent(in), optional :: partit
        !______________________________________________________________________
        integer       :: ed, el(2), enodes(2), elem, elnodes(3), n, nz, nzmin, nzmax
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: c1, c2, deltaX1, deltaX2, deltaY1, deltaY2
        real(kind=WP), dimension(:,:,:), pointer :: UV
        real(kind=WP), dimension(:),     pointer :: ssh_rhs_old

        UV          => dynamics%uv
        ssh_rhs_old => dynamics%ssh_rhs_old
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        do n = 1, nNodL
            ssh_rhs_old(n) = 0.0_WP
        end do

        do ed = 1, nEdgeO
            enodes = mesh%edges(:, ed)
            el     = mesh%edge_tri(:, ed)

            !__________________________________________________________________
            ! depth integral div(int(U_n)dz) for el(1)
            c1      = 0.0_WP
            deltaX1 = mesh%edge_cross_dxdy(1, ed)
            deltaY1 = mesh%edge_cross_dxdy(2, ed)
            nzmin   = mesh%ulevels(el(1))
            nzmax   = mesh%nlevels(el(1)) - 1
            do nz = nzmin, nzmax
                c1 = c1 + (UV(2,nz,el(1))*deltaX1 - UV(1,nz,el(1))*deltaY1)*mesh%helem(nz,el(1))
            end do

            !__________________________________________________________________
            ! ... and for el(2), unless ed is a boundary edge
            c2 = 0.0_WP
            if (el(2) > 0) then
                deltaX2 = mesh%edge_cross_dxdy(3, ed)
                deltaY2 = mesh%edge_cross_dxdy(4, ed)
                nzmin   = mesh%ulevels(el(2))
                nzmax   = mesh%nlevels(el(2)) - 1
                do nz = nzmin, nzmax
                    c2 = c2 - (UV(2,nz,el(2))*deltaX2 - UV(1,nz,el(2))*deltaY2)*mesh%helem(nz,el(2))
                end do
            end if

            ssh_rhs_old(enodes(1)) = ssh_rhs_old(enodes(1)) + (c1+c2)
            ssh_rhs_old(enodes(2)) = ssh_rhs_old(enodes(2)) - (c1+c2)
        end do

        ! linfs: the water_flux term + exchange_nod(ssh_rhs_old) are skipped (.not. linfs)

        do n = 1, nNodL
            mesh%hbar_old(n) = mesh%hbar(n)
        end do

        do n = 1, nNodO
            if (mesh%ulevels_nod2D(n) > 1) cycle          ! cavity node: hbar == hbar_old
            mesh%hbar(n) = mesh%hbar_old(n) + ssh_rhs_old(n)*dt/mesh%areasvol(mesh%ulevels_nod2D(n), n)
        end do
        if (is_multirank(partit)) call exchange_nod(mesh%hbar, partit)   ! FESOM2 :2293

        do elem = 1, nElemO
            elnodes = mesh%elem2D_nodes(1:3, elem)         ! L15
            if (mesh%ulevels(elem) > 1) then
                mesh%dhe(elem) = 0.0_WP
            else
                mesh%dhe(elem) = sum(mesh%hbar(elnodes) - mesh%hbar_old(elnodes))/3.0_WP
            end if
        end do
    end subroutine compute_hbar_ale

    !===========================================================================
    subroutine update_eta_n(dynamics, mesh, partit)
        ! FESOM2 oce_ale.F90:3973-3977 (inline in the step). Current dynamic elevation:
        !   eta_n = alpha*hbar + (1-alpha)*hbar_old   (only where ulevels_nod2D==1, i.e.
        ! no rigid-lid cavity). On pi alpha=1 -> eta_n = hbar (the new elevation).
        ! M2.12c-3: optional partit -> the FESOM2 owned+HALO node loop (:3974 do node=1,
        ! myDim_nod2D+eDim_nod2D); both hbar (compute_hbar_ale's exchange_nod) and hbar_old
        ! (its owned+halo copy) are halo-valid, so eta_n is computed at owned+halo with no
        ! exchange (matches FESOM2 — no exchange_nod(eta_n) after this loop).
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(in),    target :: mesh
        type(t_partit), intent(in), optional :: partit
        !______________________________________________________________________
        integer :: node, nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        do node = 1, nNodL
            if (mesh%ulevels_nod2D(node) == 1) &
                dynamics%eta_n(node) = alpha*mesh%hbar(node) + (1.0_WP-alpha)*mesh%hbar_old(node)
        end do
    end subroutine update_eta_n

    !===========================================================================
    subroutine vert_vel_ale(dynamics, mesh, dt, partit)
        ! FESOM2 oce_ale.F90:2323-2884, the which_ale='linfs' path. Vertical velocity from
        ! the layer-thickness divergence:
        !   Wvel(:,n)  = 0
        !   edge loop  -> Wvel = div(UV*h) per layer (edge-scatter, both signs)
        !   cumsum up  -> Wvel(nz) = Wvel(nz) + Wvel(nz+1)        (W_k = W_{k+1} - div(h_k u_k))
        !   /area      -> physical m/s
        ! linfs: dh/dt=0, so NO zlevel/zstar surface correction and hnode_new stays = hnode.
        ! Then the CFL diagnostic + the explicit/implicit Wvel split.
        ! M2.12c-3: optional partit -> FESOM2 :2381-2537 bounds: Wvel zero owned+halo, the
        ! div(UV*h) edge-scatter over OWNED edges (owned-node Wvel complete, invariant i),
        ! the cumulative sum + the area divide over OWNED nodes; then exchange_nod(Wvel) +
        ! exchange_nod(hnode_new) (:2870-2871) so compute_CFLz/compute_Wvel_split (which loop
        ! owned+halo) read a complete halo. The owned Wvel the ALE dump emits is computed
        ! directly (no halo dependence). hnode_new stays = hnode (linfs); its exchange is
        ! value-neutral but kept for faithfulness.
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(inout), target :: mesh
        real(kind=WP), intent(in) :: dt
        type(t_partit), intent(in), optional :: partit
        !______________________________________________________________________
        integer       :: ed, el(2), enodes(2), n, nz, nzmin, nzmax
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: deltaX1, deltaY1, deltaX2, deltaY2
        ! c1 is an array over levels (FESOM2 keeps it an array — a deadlock-avoidance note
        ! under OpenMP; here it just carries the per-level edge flux for the slice scatter).
        real(kind=WP) :: c1(mesh%nl-1)
        real(kind=WP), dimension(:,:,:), pointer :: UV
        real(kind=WP), dimension(:,:),   pointer :: Wvel

        UV   => dynamics%uv
        Wvel => dynamics%w
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        !______________________________________________________________________
        ! zero the vertical velocity over the whole column (below-bottom stays 0)
        do n = 1, nNodL
            Wvel(:, n) = 0.0_WP
        end do

        !______________________________________________________________________
        ! contributions from levels in the divergence (edge-scatter of div(UV*h))
        do ed = 1, nEdgeO
            enodes = mesh%edges(:, ed)
            el     = mesh%edge_tri(:, ed)

            deltaX1 = mesh%edge_cross_dxdy(1, ed)
            deltaY1 = mesh%edge_cross_dxdy(2, ed)
            nzmin   = mesh%ulevels(el(1))
            nzmax   = mesh%nlevels(el(1)) - 1
            do nz = nzmax, nzmin, -1
                c1(nz) = (UV(2,nz,el(1))*deltaX1 - UV(1,nz,el(1))*deltaY1)*mesh%helem(nz,el(1))
            end do
            Wvel(nzmin:nzmax, enodes(1)) = Wvel(nzmin:nzmax, enodes(1)) + c1(nzmin:nzmax)
            Wvel(nzmin:nzmax, enodes(2)) = Wvel(nzmin:nzmax, enodes(2)) - c1(nzmin:nzmax)

            c1 = 0.0_WP
            if (el(2) > 0) then
                deltaX2 = mesh%edge_cross_dxdy(3, ed)
                deltaY2 = mesh%edge_cross_dxdy(4, ed)
                nzmin   = mesh%ulevels(el(2))
                nzmax   = mesh%nlevels(el(2)) - 1
                do nz = nzmax, nzmin, -1
                    c1(nz) = -(UV(2,nz,el(2))*deltaX2 - UV(1,nz,el(2))*deltaY2)*mesh%helem(nz,el(2))
                end do
                Wvel(nzmin:nzmax, enodes(1)) = Wvel(nzmin:nzmax, enodes(1)) + c1(nzmin:nzmax)
                Wvel(nzmin:nzmax, enodes(2)) = Wvel(nzmin:nzmax, enodes(2)) - c1(nzmin:nzmax)
            end if
        end do

        !______________________________________________________________________
        ! cumulative summation of div(UV*h) vertically: W_k = W_k + W_{k+1}
        do n = 1, nNodO
            nzmin = mesh%ulevels_nod2D(n)
            nzmax = mesh%nlevels_nod2D(n) - 1
            do nz = nzmax, nzmin, -1
                Wvel(nz, n) = Wvel(nz, n) + Wvel(nz+1, n)
            end do
        end do

        !______________________________________________________________________
        ! divide by depth-dependent cell area -> physical vertical velocity (m/s)
        do n = 1, nNodO
            nzmin = mesh%ulevels_nod2D(n)
            nzmax = mesh%nlevels_nod2D(n) - 1
            do nz = nzmin, nzmax
                Wvel(nz, n) = Wvel(nz, n)/mesh%area(nz, n)
            end do
        end do

        ! linfs: no zlevel/zstar free-surface correction; hnode_new unchanged (= hnode).
        if (is_multirank(partit)) then
            call exchange_nod(Wvel, partit)             ! FESOM2 :2870
            call exchange_nod(mesh%hnode_new, partit)   ! FESOM2 :2871
        end if

        call compute_CFLz(dynamics, mesh, dt, partit)
        call compute_Wvel_split(dynamics, mesh, partit)
    end subroutine vert_vel_ale

    !===========================================================================
    subroutine compute_CFLz(dynamics, mesh, dt, partit)
        ! FESOM2 oce_ale.F90:3126-3213. Vertical CFL per layer (used by the Wvel split):
        !   CFL_z(nz)   += |Wvel(nz)  *dt/hnode_new(nz)|
        !   CFL_z(nz+1)  = |Wvel(nz+1)*dt/hnode_new(nz)|
        ! The two assignments are kept SEPARATE (not folded into one accumulate) — FESOM2
        ! notes this exact form is "for the sake of reproducibility ... (rounding error)"
        ! (L19: the reduction/accumulation form is part of the bits). cflmax + the warning
        ! prints are diagnostics, dropped.
        ! M2.12c-3: optional partit -> the FESOM2 owned+HALO node loops (:3154/:3161 do node=
        ! 1, myDim_nod2D+eDim_nod2D); Wvel/hnode_new are halo-valid from vert_vel_ale's two
        ! exchange_nod's, so CFL_z is built over owned+halo (feeds compute_Wvel_split's
        ! owned+halo Wvel_e/Wvel_i, which the NEXT step's momentum advection reads at halo).
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(in),    target :: mesh
        real(kind=WP), intent(in) :: dt
        type(t_partit), intent(in), optional :: partit
        !______________________________________________________________________
        integer       :: node, nz, nzmin, nzmax
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: c1, c2
        real(kind=WP), dimension(:,:), pointer :: Wvel, CFL_z

        Wvel  => dynamics%w
        CFL_z => dynamics%cfl_z
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        do node = 1, nNodL
            CFL_z(:, node) = 0._WP
        end do

        do node = 1, nNodL
            nzmin = mesh%ulevels_nod2D(node)
            nzmax = mesh%nlevels_nod2D(node) - 1
            do nz = nzmin, nzmax
                c1 = abs(Wvel(nz,   node)*dt/mesh%hnode_new(nz, node))
                c2 = abs(Wvel(nz+1, node)*dt/mesh%hnode_new(nz, node))
                CFL_z(nz,   node) = CFL_z(nz, node) + c1
                CFL_z(nz+1, node) = c2
            end do
        end do
    end subroutine compute_CFLz

    !===========================================================================
    subroutine compute_Wvel_split(dynamics, mesh, partit)
        ! FESOM2 oce_ale.F90:3217-3265. Split Wvel into explicit (Wvel_e) and implicit
        ! (Wvel_i) parts according to the vertical CFL, so the explicit part is capped at
        ! wsplit_maxcfl. use_wsplit=.true. + wsplit_maxcfl=1.0 on pi. OVERWRITES dynamics%
        ! w_e/w_i (the same arrays prescribed as M2.4/M2.5 inputs — fully consumed before
        ! this runs; in the real timestep they feed the NEXT step's momadv/ivertvisc).
        ! M2.12c-3: optional partit -> the FESOM2 owned+HALO node loop (:3251 do node=1,
        ! myDim_nod2D+eDim_nod2D); Wvel/CFL_z are halo-valid so Wvel_e/Wvel_i are produced
        ! at owned+halo with no trailing exchange (matches FESOM2).
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(in),    target :: mesh
        type(t_partit), intent(in), optional :: partit
        !______________________________________________________________________
        integer       :: node, nz, nzmin, nzmax
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: dd
        real(kind=WP), dimension(:,:), pointer :: Wvel, Wvel_e, Wvel_i, CFL_z

        Wvel   => dynamics%w
        Wvel_e => dynamics%w_e
        Wvel_i => dynamics%w_i
        CFL_z  => dynamics%cfl_z
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        do node = 1, nNodL
            nzmin = mesh%ulevels_nod2D(node)
            nzmax = mesh%nlevels_nod2D(node)
            do nz = nzmin, nzmax
                Wvel_e(nz, node) = Wvel(nz, node)
                Wvel_i(nz, node) = 0.0_WP
                if (dynamics%use_wsplit .and. (CFL_z(nz, node) > dynamics%wsplit_maxcfl)) then
                    dd = max((CFL_z(nz, node)-dynamics%wsplit_maxcfl), 0.0_WP)/max(dynamics%wsplit_maxcfl, 1.e-12_WP)
                    Wvel_e(nz, node) = (1.0_WP/(1.0_WP+dd))*Wvel(nz, node)   ! explicit (=1 if dd=0)
                    Wvel_i(nz, node) = (dd    /(1.0_WP+dd))*Wvel(nz, node)   ! implicit (=1 if dd=inf)
                end if
            end do
        end do
    end subroutine compute_Wvel_split

    !===========================================================================
    subroutine update_thickness_ale(mesh, partit)
        ! FESOM2 oce_ale.F90:1226-1444. Commit the new layer thicknesses at the end of
        ! the step: hnode=hnode_new and recompute helem / zbar_3d_n / Z_3d_n. For
        ! which_ale='linfs' (the M2 path) the layer thickness is fixed (dh/dt=0, hnode_new
        ! stays = hnode), so NEITHER the zlevel NOR the zstar redistribution branch is
        ! taken — hnode/helem/zbar_3d_n/Z_3d_n are unchanged and the only remaining action
        ! is exchange_elem(helem) (:1440). The full zlevel/zstar thickness commit (the
        ! rescue_hnode_old DVD bookkeeping + the per-element helem average) enters with
        ! non-linfs ALE in a later gate.
        ! M2.12c-3: optional partit accepted for the step_oce threading; for linfs helem is
        ! unchanged AND already full-halo-valid (built from nlevels over the full local range
        ! at setup), so FESOM2's exchange_elem(helem) is value-neutral and is SKIPPED here
        ! (the L33 "don't add machinery you don't need" rule — no kernel reads a stale helem
        ! halo that this would fix). It enters with the real zlevel/zstar thickness commit.
        type(t_mesh), intent(in) :: mesh
        type(t_partit), intent(in), optional :: partit
        if (mesh%nod2D < 0 .and. present(partit)) return   ! silence unused-args; linfs = no-op
        ! linfs: no thickness redistribution; exchange_elem(helem) value-neutral (see above).
    end subroutine update_thickness_ale

end module oce_ale
