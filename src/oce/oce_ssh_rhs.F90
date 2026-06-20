module oce_ssh_rhs
    ! SSH (free-surface) implicit-solve assembly: the elevation stiffness matrix and
    ! the elevation right-hand side. Run in the timestep AFTER impl_vert_visc_ale and
    ! BEFORE the CG solve (FESOM2 oce_ale.F90:3920-3930, the use_ssh_se_subcycl=.false.
    ! branch). Transcribed from FESOM2 v2.7.3 oce_ale.F90:
    !   init_stiff_mat_ale   (:1584-1875) — builds the CSR stiffness operator once
    !   compute_ssh_rhs_ale  (:2012-2147) — assembles ssh_rhs each step
    !
    ! The split-implicit free-surface equation (eq. 18 of "FESOM2: from finite elements
    ! to finite volumes") is
    !       [ M/dt - alpha*theta*g*dt*div(H grad) ] d_eta = ssh_rhs
    ! a symmetric sparse linear system in the SSH increment d_eta, solved by the
    ! preconditioned CG in oce_ssh_solve. This module builds the left-hand-side matrix
    ! (init_stiff_mat_ale) and the right-hand side (compute_ssh_rhs_ale).
    !
    ! STIFFNESS MATRIX (init_stiff_mat_ale):
    !  - CSR sparsity from the node neighbourhood (self + edge-connected nodes); the
    !    diagonal (self) is the FIRST entry of every row (n_pos(1,n)=n).
    !  - Stiffness part (loop over edges, both adjacent elements): the H*div operator
    !    factor*(zbar_e_bot-zbar_e_srf)*(gradient_sca_x*edge_cross_dy - gradient_sca_y*
    !    edge_cross_dx), scattered +/- into the two edge-node rows, factor=g*dt*alpha*theta.
    !  - Mass part: + areasvol(surf)/dt added to each row diagonal.
    !  - For which_ale=='linfs' the matrix is built ONCE here with the unperturbed depth
    !    and never updated (FESOM2 oce_ale.F90:3921 calls update_stiff_mat_ale only for
    !    NON-linfs). So this single assembly IS the production matrix on pi/reduced-M2.
    !
    ! SSH RHS (compute_ssh_rhs_ale): the depth-integrated horizontal-divergence of the
    ! provisional transport alpha*(UV+UV_rhs) over the scalar control volume, as an
    ! edge flux scattered +/- into the two edge nodes; plus the linfs (1-alpha)*ssh_rhs_old
    ! term (which vanishes here since alpha=1, theta=1 on pi/reduced-M2). UV_rhs is the
    ! post-impl_vert_visc_ale velocity rhs (the M2.5-gated uv_rhs_ivv).
    !
    ! BIT-IDENTITY NOTES (the L9 transitive-gate pattern again — every operand already
    ! byte-pinned, so faithful transcription byte-matches like M1.1-M2.5):
    !  - The stiffness value scatter (+= fy*factor over edges) and the ssh_rhs scatter
    !    (+= (c1+c2) over edges) are FP-order-sensitive, but the 1-rank global edge order
    !    + edge_tri + edges + gradient_sca + edge_cross_dxdy + zbar_e_bot + areasvol +
    !    helem are all geometry/M2.5-proven (max|delta|=0), so the accumulation order
    !    matches FESOM2 exactly. The FESOM2 oracle is ENABLE_OPENMP=OFF (L16), so its
    !    !$OMP ORDERED / omp_set_lock edge-scatter compiles to a plain serial loop in
    !    edge order — identical to this serial transcription.
    !  - factor = g*dt*alpha*theta and areasvol/dt are runtime products/quotients with
    !    byte-identical operands on both sides; the Intel -no-prec-div reciprocal matches
    !    (L7/L10/L14). dt MUST be the pi namelist dt=86400/36 (the value FESOM2's
    !    init_stiff_mat_ale saw at ocean_setup), computed the SAME way (not a 2400.0
    !    literal — see the driver) so the bits agree.
    !  - zbar_e_srf = zbar(ulevels(elem)) = zbar(1) = 0 on pi (no cavity; FESOM2 default
    !    oce_ale.F90:525). Computed inline (no new mesh field): faithful to the non-cavity
    !    formula, and 0 on pi.
    !  - elem2D_nodes(1:3,el) sliced per the L15 MAX_NV=4-vs-3 trap.
    !
    ! 1-RANK SCOPE (multi-rank lifted at M2.12):
    !  - The CSR is built in LOCAL numbering (rowptr_loc/colind_loc, used by the CG); at
    !    1-rank local==global natural, so the FESOM2 global-contiguous remap (rpart.out
    !    mapping, the per-PE nza offset) is the identity and is dropped. rowptr/colind are
    !    set equal to rowptr_loc/colind_loc (their 1-rank values).
    !  - myDim_nod2D -> nod2D, myDim_edge2D -> edge2D, eDim_nod2D -> 0. The
    !    row<=myDim_nod2D / n1<=myDim_nod2D guards are always true at 1-rank (kept for
    !    multi-rank faithfulness). exchange_nod(ssh_rhs) at the tail is a no-op, dropped.
    use mod_precision,  only: WP, MP
    use mod_constants,  only: g
    use mod_param_phys, only: alpha, theta
    use mod_mesh,       only: t_mesh
    use mod_dyn,        only: t_dyn
    implicit none
    private
    public :: init_stiff_mat_ale, compute_ssh_rhs_ale

contains

    !===========================================================================
    subroutine init_stiff_mat_ale(mesh, dt)
        ! Build mesh%ssh_stiff (CSR sparsity + the linfs stiffness/mass values). Call
        ! ONCE; dt is the pi namelist timestep (see module/driver notes). Idempotent
        ! it is NOT (allocates ssh_stiff) — call exactly once per mesh.
        type(t_mesh),  intent(inout), target :: mesh
        real(kind=WP), intent(in)            :: dt
        !______________________________________________________________________
        integer              :: n, n1, n2, i, row, ed
        integer              :: elnodes(3), el(2), npos(3)
        integer              :: offset, nini, nend
        real(kind=WP)        :: factor, fy(3), zsrf
        integer, allocatable :: n_num(:), n_pos(:,:)

        associate(ssh_stiff => mesh%ssh_stiff)
        !__________________________________________________________________
        ! a) neighbourhood: n_num(n) = #neighbours of node n (incl. self),
        !    n_pos(:,n) = their local indices, n_pos(1,n)=n (self -> diagonal first)
        allocate(n_num(mesh%nod2D), n_pos(12, mesh%nod2D))
        n_pos = 0
        do n = 1, mesh%nod2D
            n_num(n)   = 1
            n_pos(1,n) = n
        end do
        do n = 1, mesh%edge2D
            n1 = mesh%edges(1,n)
            n2 = mesh%edges(2,n)
            if (n1 <= mesh%nod2D) then
                n_pos(n_num(n1)+1, n1) = n2
                n_num(n1) = n_num(n1)+1
            end if
            if (n2 <= mesh%nod2D) then
                n_pos(n_num(n2)+1, n2) = n1
                n_num(n2) = n_num(n2)+1
            end if
        end do

        !__________________________________________________________________
        ! b) CSR row pointers + nonzero count
        ssh_stiff%dim = mesh%nod2D
        allocate(ssh_stiff%rowptr(mesh%nod2D+1), ssh_stiff%rowptr_loc(mesh%nod2D+1))
        ssh_stiff%rowptr_loc(1) = 1
        do n = 1, mesh%nod2D
            ssh_stiff%rowptr_loc(n+1) = ssh_stiff%rowptr_loc(n) + n_num(n)
        end do
        ssh_stiff%nza = ssh_stiff%rowptr_loc(mesh%nod2D+1) - 1

        !__________________________________________________________________
        ! c) CSR column indices (local) + zero the values
        allocate(ssh_stiff%colind(ssh_stiff%nza), ssh_stiff%colind_loc(ssh_stiff%nza))
        allocate(ssh_stiff%values(ssh_stiff%nza))
        ssh_stiff%values = 0.0_WP
        do n = 1, mesh%nod2D
            nini = ssh_stiff%rowptr_loc(n)
            nend = ssh_stiff%rowptr_loc(n+1) - 1
            ssh_stiff%colind_loc(nini:nend) = n_pos(1:n_num(n), n)
        end do
        ! 1-rank: global == local natural numbering
        ssh_stiff%rowptr = ssh_stiff%rowptr_loc
        ssh_stiff%colind = ssh_stiff%colind_loc

        !__________________________________________________________________
        ! d) stiffness part: factor * H * div, scattered over edges. n_num is reused
        !    as the reverse-map (local node index -> its CSR position within the row).
        n_num  = 0
        factor = g*dt*alpha*theta
        do ed = 1, mesh%edge2D
            el = mesh%edge_tri(:, ed)
            do i = 1, 2   ! the two triangles sharing edge ed
                if (el(i) < 1) cycle   ! boundary edge has only one triangle
                elnodes = mesh%elem2D_nodes(1:3, el(i))
                zsrf    = mesh%zbar(mesh%ulevels(el(i)))   ! zbar_e_srf (=0, no cavity)
                fy(1:3) = (mesh%zbar_e_bot(el(i)) - zsrf) * &
                          ( mesh%gradient_sca(1:3,el(i)) * mesh%edge_cross_dxdy(2*i  ,ed)  &
                           -mesh%gradient_sca(4:6,el(i)) * mesh%edge_cross_dxdy(2*i-1,ed) )
                if (i==2) fy = -fy

                row = mesh%edges(1, ed)
                if (row <= mesh%nod2D) then
                    do n = ssh_stiff%rowptr_loc(row), ssh_stiff%rowptr_loc(row+1)-1
                        n_num(ssh_stiff%colind_loc(n)) = n
                    end do
                    npos = n_num(elnodes)
                    ssh_stiff%values(npos) = ssh_stiff%values(npos) + fy*factor
                end if

                row = mesh%edges(2, ed)
                if (row <= mesh%nod2D) then
                    do n = ssh_stiff%rowptr_loc(row), ssh_stiff%rowptr_loc(row+1)-1
                        n_num(ssh_stiff%colind_loc(n)) = n
                    end do
                    npos = n_num(elnodes)
                    ssh_stiff%values(npos) = ssh_stiff%values(npos) - fy*factor
                end if
            end do
        end do

        !__________________________________________________________________
        ! e) mass part: + areasvol(surf)/dt on the row diagonal (the first CSR entry).
        !    Skip cavity nodes (ulevels_nod2D>1: rigid-lid, no eta time-derivative).
        do row = 1, mesh%nod2D
            if (mesh%ulevels_nod2D(row) > 1) cycle
            offset = ssh_stiff%rowptr_loc(row)
            ssh_stiff%values(offset) = ssh_stiff%values(offset) &
                                     + mesh%areasvol(mesh%ulevels_nod2D(row),row)/dt
        end do

        deallocate(n_pos, n_num)
        end associate
    end subroutine init_stiff_mat_ale

    !===========================================================================
    subroutine compute_ssh_rhs_ale(dynamics, mesh)
        ! Assemble dynamics%ssh_rhs = depth-integrated horizontal divergence of
        ! alpha*(UV+UV_rhs), scattered as an edge flux into the two edge nodes
        ! (+ the linfs (1-alpha)*ssh_rhs_old term, = 0 on pi since alpha=1).
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(in),    target :: mesh
        !______________________________________________________________________
        integer       :: ed, el(2), enodes(2), nz, n, nzmin, nzmax
        real(kind=WP) :: c1, c2, deltaX1, deltaX2, deltaY1, deltaY2
        real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
        real(kind=WP), dimension(:),     pointer :: ssh_rhs, ssh_rhs_old

        UV          => dynamics%uv
        UV_rhs      => dynamics%uv_rhs
        ssh_rhs     => dynamics%ssh_rhs
        ssh_rhs_old => dynamics%ssh_rhs_old

        do n = 1, mesh%nod2D
            ssh_rhs(n) = 0.0_WP
        end do

        do ed = 1, mesh%edge2D
            enodes = mesh%edges(:, ed)
            el     = mesh%edge_tri(:, ed)

            !__________________________________________________________________
            ! alpha*nabla*int(U_n + U_rhs) dz contribution from el(1)
            c1      = 0.0_WP
            deltaX1 = mesh%edge_cross_dxdy(1, ed)
            deltaY1 = mesh%edge_cross_dxdy(2, ed)
            nzmin   = mesh%ulevels(el(1))
            nzmax   = mesh%nlevels(el(1)) - 1
            do nz = nzmin, nzmax
                c1 = c1 + alpha*((UV(2,nz,el(1))+UV_rhs(2,nz,el(1)))*deltaX1 -  &
                                 (UV(1,nz,el(1))+UV_rhs(1,nz,el(1)))*deltaY1)*mesh%helem(nz,el(1))
            end do

            !__________________________________________________________________
            ! ... and from el(2), unless ed is a boundary edge
            c2 = 0.0_WP
            if (el(2) > 0) then
                deltaX2 = mesh%edge_cross_dxdy(3, ed)
                deltaY2 = mesh%edge_cross_dxdy(4, ed)
                nzmin   = mesh%ulevels(el(2))
                nzmax   = mesh%nlevels(el(2)) - 1
                do nz = nzmin, nzmax
                    c2 = c2 - alpha*((UV(2,nz,el(2))+UV_rhs(2,nz,el(2)))*deltaX2 -  &
                                     (UV(1,nz,el(2))+UV_rhs(1,nz,el(2)))*deltaY2)*mesh%helem(nz,el(2))
                end do
            end if

            !__________________________________________________________________
            ! net "flux" scattered into the two edge nodes
            ssh_rhs(enodes(1)) = ssh_rhs(enodes(1)) + (c1+c2)
            ssh_rhs(enodes(2)) = ssh_rhs(enodes(2)) - (c1+c2)
        end do

        ! linfs water-flux term: ssh_rhs += (1-alpha)*ssh_rhs_old (= 0 since alpha=1).
        ! The non-linfs water_flux branch (zstar) is deferred with its own ALE gate.
        do n = 1, mesh%nod2D
            ssh_rhs(n) = ssh_rhs(n) + (1.0_WP-alpha)*ssh_rhs_old(n)
        end do
        ! exchange_nod(ssh_rhs) — 1-rank no-op (lifted at M2.12)
    end subroutine compute_ssh_rhs_ale

end module oce_ssh_rhs
