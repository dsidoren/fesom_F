module oce_ssh_solve
    ! Preconditioned conjugate-gradient solve of the free-surface system
    !       ssh_stiff * d_eta = ssh_rhs
    ! built by oce_ssh_rhs. Run in the timestep after compute_ssh_rhs_ale (FESOM2
    ! oce_ale.F90:3930 -> solve_ssh_ale). Transcribed from FESOM2 v2.7.3:
    !   solve_ssh_ale            (oce_ale.F90:3272-3311) — driver (precond once, then CG)
    !   ssh_solve_preconditioner (solver.F90:31-95)      — MITgcm Jacobi-symmetrised M^-1
    !   ssh_solve_cg             (solver.F90:98-281)      — the CG iteration
    !
    ! PRECONDITIONER (MITgcm JGR 102,5753-5766,1997): the inverse-preconditioner matrix
    ! K (one symmetrised Jacobi sweep) shares the ssh_stiff CSR pattern; row r has
    !   diagonal  K_rr = 1/a_rr
    !   off-diag  K_ri = -0.5*(a_ri/a_rr)/(a_rr + a_ii)
    ! where a_** are ssh_stiff entries and a_ii the diagonal in the neighbour row i.
    ! Symmetrisation (the 0.5 + the a_ii) makes K SPD so CG applies.
    !
    ! CG (standard preconditioned CG; the FESOM2 comment notes it follows Wikipedia,
    ! not the MITgcm beta): solves A x = b for x = d_eta, x0 = the incoming d_eta.
    !   r = b - A x ; z = M^-1 r ; p = z ; rho = r.z
    !   loop: Ap = A p ; al = rho/(p.Ap) ; x += al p ; r -= al Ap ; z = M^-1 r
    !         rho_new = r.z ; if sqrt((r.r)/nod2D) < rtol exit
    !         be = rho_new/rho ; rho = rho_new ; p = z + be p
    !   rtol = soltol*sqrt((b.b)/nod2D), soltol=1e-5.
    !
    ! BIT-IDENTITY NOTES:
    !  - Given byte-identical operands (the M2.6-gated ssh_stiff + ssh_rhs) the CG is a
    !    deterministic recurrence: same matrix+rhs+x0 -> same per-iteration scalars ->
    !    same iteration count (the convergence test fires at the same iter) -> same x.
    !    So the gate reduces to "are A and b byte-identical?" — the L9 transitive pattern.
    !  - REDUCTION FORM matters for the bits. The FESOM2 oracle is ENABLE_OPENMP=OFF and
    !    __openmp_reproducible is NOT defined, so its dot-products compile to the explicit
    !    serial `DO row; s=s+...; END DO` (the !$OMP REDUCTION clause is an inert comment).
    !    This file transcribes that serial DO form for s_old/s_aux/sprod — NOT a sum()
    !    intrinsic (which the compiler may reduce in a different order). The matrix-vector
    !    products DO use the sum() intrinsic over a CSR slice, exactly as FESOM2 — both
    !    codes use sum() there, same flags, byte-identical operands -> byte-identical.
    !  - The many runtime divisors (1/a_rr in the precond, al=rho/s_aux, be) are
    !    byte-identical on both sides (same operands), so -no-prec-div matches (L7/L14).
    !
    ! 1-RANK SCOPE (multi-rank lifted at M2.12):
    !  - exchange_nod(diag_values/rr/pp/x/d_eta) are halo broadcasts -> no-ops at 1-rank,
    !    dropped. The MPI_Allreduce of s_old/s_aux/sprod over 1 rank is the identity (the
    !    local serial sum IS the global sum), dropped. So the CG here is purely local +
    !    deterministic, matching the FESOM2 1-rank reduction order. The new multi-rank
    !    bit-identity risk (the cross-rank dot-product reduction order) is an M2.12 concern.
    !  - rr/zz/pp/App sized nod2D (no eDim halo); CSR via colind_loc/rowptr_loc (local).
    use mod_precision, only: WP, MP
    use mod_mesh,      only: t_mesh
    use mod_dyn,       only: t_dyn
    implicit none
    private
    public :: solve_ssh_ale

contains

    !===========================================================================
    subroutine solve_ssh_ale(dynamics, mesh, n_iter)
        ! Build the preconditioner on the first call, then CG-solve for d_eta.
        ! n_iter (optional, NOT part of FESOM2's signature) returns the CG iteration
        ! count for the driver's non-vacuity diagnostic only.
        type(t_dyn),  intent(inout), target  :: dynamics
        type(t_mesh), intent(inout), target  :: mesh
        integer, intent(out), optional       :: n_iter
        logical, save :: lfirst = .true.

        if (lfirst) call ssh_solve_preconditioner(dynamics, mesh)
        call ssh_solve_cg(dynamics%d_eta, dynamics%ssh_rhs, dynamics, mesh, n_iter)
        ! exchange_nod(d_eta) — 1-rank no-op
        lfirst = .false.
    end subroutine solve_ssh_ale

    !===========================================================================
    subroutine ssh_solve_preconditioner(dynamics, mesh)
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(inout), target :: mesh
        integer                    :: nend, row, node, n, offset
        real(kind=WP), allocatable :: diag_values(:)

        associate(ssh_stiff => mesh%ssh_stiff)
        nend = ssh_stiff%rowptr_loc(mesh%nod2D+1) - ssh_stiff%rowptr_loc(1)
        allocate(ssh_stiff%pr_values(nend))
        allocate(diag_values(mesh%nod2D))   ! 1-rank: no eDim halo

        ! diagonal of A (the first CSR entry of each row)
        do row = 1, mesh%nod2D
            offset = ssh_stiff%rowptr_loc(row) - ssh_stiff%rowptr_loc(1) + 1
            diag_values(row) = ssh_stiff%values(offset)
        end do
        ! exchange_nod(diag_values) — 1-rank no-op

        ! fill the inverse-preconditioner values
        do row = 1, mesh%nod2D
            offset = ssh_stiff%rowptr_loc(row) - ssh_stiff%rowptr_loc(1)
            nend   = ssh_stiff%rowptr_loc(row+1) - ssh_stiff%rowptr_loc(row)
            ssh_stiff%pr_values(offset+1) = 1.0_WP/ssh_stiff%values(offset+1)
            do n = 2, nend
                node = ssh_stiff%colind_loc(offset+n)
                ssh_stiff%pr_values(n+offset) =                                          &
                    -0.5_WP*(ssh_stiff%values(n+offset)/ssh_stiff%values(1+offset))      &
                    / (ssh_stiff%values(1+offset) + diag_values(node))
            end do
        end do
        deallocate(diag_values)

        n = mesh%nod2D   ! 1-rank: no eDim halo
        allocate(dynamics%solverinfo%rr(n),  dynamics%solverinfo%zz(n),  &
                 dynamics%solverinfo%pp(n),  dynamics%solverinfo%App(n))
        dynamics%solverinfo%rr  = 0.0_WP
        dynamics%solverinfo%zz  = 0.0_WP
        dynamics%solverinfo%pp  = 0.0_WP
        dynamics%solverinfo%App = 0.0_WP
        end associate
    end subroutine ssh_solve_preconditioner

    !===========================================================================
    subroutine ssh_solve_cg(x, rhs, dynamics, mesh, n_iter)
        real(kind=WP), intent(inout)         :: x(:)     ! d_eta (in: x0, out: solution)
        real(kind=WP), intent(in)            :: rhs(:)   ! ssh_rhs
        type(t_dyn),  intent(inout), target  :: dynamics
        type(t_mesh), intent(inout), target  :: mesh
        integer, intent(out), optional       :: n_iter
        integer                  :: row, iter
        real(kind=WP)            :: sprod(2), s_old, s_aux, al, be, rtol
        real(kind=MP), pointer   :: values(:), pr_values(:)
        real(kind=WP), pointer   :: rr(:), zz(:), pp(:), App(:)
        integer,       pointer   :: rptr(:), cind(:)

        values    => mesh%ssh_stiff%values
        pr_values => mesh%ssh_stiff%pr_values
        cind      => mesh%ssh_stiff%colind_loc
        rptr      => mesh%ssh_stiff%rowptr_loc
        rr  => dynamics%solverinfo%rr
        zz  => dynamics%solverinfo%zz
        pp  => dynamics%solverinfo%pp
        App => dynamics%solverinfo%App

        !__________________________________________________________________
        ! working tolerance rtol = soltol*sqrt((b.b)/nod2D)
        s_old = 0.0_WP
        do row = 1, mesh%nod2D
            s_old = s_old + rhs(row)*rhs(row)
        end do
        ! MPI_Allreduce(s_old) — 1-rank identity
        rtol = dynamics%solverinfo%soltol*sqrt(s_old/real(mesh%nod2D,WP))

        !__________________________________________________________________
        ! r0 = b - A x0
        do row = 1, mesh%nod2D
            rr(row) = rhs(row) - sum(values(rptr(row):rptr(row+1)-1)*x(cind(rptr(row):rptr(row+1)-1)))
        end do
        ! exchange_nod(rr) — 1-rank no-op

        !__________________________________________________________________
        ! z0 = M^-1 r0 ; search direction pp = z0
        do row = 1, mesh%nod2D
            zz(row) = sum(pr_values(rptr(row):rptr(row+1)-1)*rr(cind(rptr(row):rptr(row+1)-1)))
            pp(row) = zz(row)
        end do

        ! rho = r0.z0
        s_old = 0.0_WP
        do row = 1, mesh%nod2D
            s_old = s_old + rr(row)*zz(row)
        end do
        ! MPI_Allreduce(s_old) — 1-rank identity

        !__________________________________________________________________
        ! iterations
        if (present(n_iter)) n_iter = dynamics%solverinfo%maxiter   ! fallback: no convergence
        do iter = 1, dynamics%solverinfo%maxiter
            ! exchange_nod(pp) — 1-rank no-op
            do row = 1, mesh%nod2D
                App(row) = sum(values(rptr(row):rptr(row+1)-1)*pp(cind(rptr(row):rptr(row+1)-1)))
            end do

            s_aux = 0.0_WP
            do row = 1, mesh%nod2D
                s_aux = s_aux + pp(row)*App(row)
            end do
            ! MPI_Allreduce(s_aux) — 1-rank identity
            al = s_old/s_aux

            do row = 1, mesh%nod2D
                x(row)  = x(row)  + al*pp(row)
                rr(row) = rr(row) - al*App(row)
            end do
            ! exchange_nod(rr) — 1-rank no-op

            do row = 1, mesh%nod2D
                zz(row) = sum(pr_values(rptr(row):rptr(row+1)-1)*rr(cind(rptr(row):rptr(row+1)-1)))
            end do

            sprod(1:2) = 0.0_WP
            do row = 1, mesh%nod2D
                sprod(1) = sprod(1) + rr(row)*zz(row)
                sprod(2) = sprod(2) + rr(row)*rr(row)
            end do
            ! MPI_Allreduce(sprod,2) — 1-rank identity

            if (sqrt(sprod(2)/mesh%nod2D) < rtol) then
                if (present(n_iter)) n_iter = iter
                exit
            end if

            be    = sprod(1)/s_old
            s_old = sprod(1)
            do row = 1, mesh%nod2D
                pp(row) = zz(row) + be*pp(row)
            end do
        end do
        ! exchange_nod(x) — 1-rank no-op
    end subroutine ssh_solve_cg

end module oce_ssh_solve
