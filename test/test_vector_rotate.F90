program test_vector_rotate
    ! M9 Task 2.5 GATE for vector_r2g (mod_mesh_rotate): prove the newly-ported rotated->geographic
    ! vector transform is the exact inverse of the already-byte-gated geographic->rotated vector_g2r
    ! (used since M2.10 forcing). vector_r2g is transcribed VERBATIM from FESOM2
    ! gen_modules_rotate_grid.F90:164-202, so it is bit-faithful to FESOM2 by construction; this test
    ! is the invertibility + magnitude-preservation + flag_coord sanity check (round-off tolerance,
    ! since the round-trip passes through sin/cos/asin/atan2 — NOT a state byte-gate).
    !
    ! Uses a NON-identity displaced-pole rotation (50,15,-90) = the classic FESOM rotated grid (and the
    ! one fesom_outputsmoke drives), so the transform is genuinely exercised (not a no-op identity).
    use mpi
    use mod_precision,   only: WP
    use mod_partit,      only: t_partit
    use mod_partitioning, only: par_init, par_ex
    use mod_mesh_rotate, only: init_mesh_rotation, g2r, r2g, vector_g2r, vector_r2g
    implicit none

    type(t_partit) :: partit
    real(kind=WP), parameter :: D2R = 3.141592653589793_WP/180.0_WP
    real(kind=WP), parameter :: TOL = 1.0e-11_WP   ! round-trip / magnitude round-off bar
    real(kind=WP) :: max_rt, max_mag, max_flag, max_id
    integer :: nfail

    call par_init(partit)
    nfail = 0

    ! ---- non-identity rotation: round-trip + magnitude + flag_coord equivalence -------------------
    call init_mesh_rotation(50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP)
    call sweep(.false., max_rt, max_mag, max_flag)
    if (partit%mype == 0) then
        write(*,'(a,es10.3,a,es10.3,a,es10.3)') 'test_vector_rotate: round-trip max|Δ|=', max_rt, &
            '  magnitude max|Δ|=', max_mag, '  flag0==flag1 max|Δ|=', max_flag
        if (max_rt  > TOL) nfail = nfail + 1
        if (max_mag > TOL) nfail = nfail + 1
        if (max_flag> TOL) nfail = nfail + 1
    end if

    ! ---- identity rotation (0,0,0): vector_r2g must be a no-op (geographic == native) -------------
    call init_mesh_rotation(0.0_WP, 0.0_WP, 0.0_WP, 360.0_WP)
    call sweep(.true., max_id, max_mag, max_flag)
    if (partit%mype == 0) then
        write(*,'(a,es10.3)') 'test_vector_rotate: identity-rotation no-op max|Δ|=', max_id
        if (max_id > TOL) nfail = nfail + 1
    end if

    if (partit%mype == 0) then
        if (nfail == 0) then
            write(*,'(a)') 'test_vector_rotate: PASS (vector_r2g inverts vector_g2r, magnitude-preserving)'
        else
            write(*,'(a,i0,a)') 'test_vector_rotate: FAIL (', nfail, ' check(s))'
        end if
    end if
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    if (nfail /= 0) error stop 1

contains

    ! Sweep a grid of rotated positions x vectors. check_identity=.true. additionally asserts that
    ! vector_r2g leaves the vector unchanged (only valid for the identity rotation matrix).
    subroutine sweep(check_identity, mx_rt, mx_mag, mx_flag)
        logical,       intent(in)  :: check_identity
        real(kind=WP), intent(out) :: mx_rt, mx_mag, mx_flag
        real(kind=WP) :: rlon, rlat, glon, glat
        real(kind=WP) :: u0, v0, ug, vg, u1, v1, u2, v2, mag0, mag1
        integer :: ilon, ilat, iv
        real(kind=WP), dimension(5) :: uu = [1.0_WP, 0.0_WP,  3.0_WP, -5.0_WP, 10.0_WP]
        real(kind=WP), dimension(5) :: vv = [0.0_WP, 1.0_WP, -2.0_WP,  4.0_WP, 10.0_WP]
        mx_rt = 0.0_WP; mx_mag = 0.0_WP; mx_flag = 0.0_WP
        do ilat = -8, 8                     ! rotated lat in [-80,80] deg (avoid the poles)
            rlat = real(ilat, WP)*10.0_WP*D2R
            do ilon = -17, 17               ! rotated lon in [-170,170] deg
                rlon = real(ilon, WP)*10.0_WP*D2R
                call r2g(glon, glat, rlon, rlat)   ! geographic coords of this rotated position
                do iv = 1, 5
                    u0 = uu(iv); v0 = vv(iv)
                    ! rotated -> geographic, flag_coord=0 (position given in rotated coords)
                    ug = u0; vg = v0
                    call vector_r2g(ug, vg, rlon, rlat, 0)
                    ! geographic -> rotated (vector_g2r flag_coord=0 also takes rotated coords): inverse
                    u1 = ug; v1 = vg
                    call vector_g2r(u1, v1, rlon, rlat, 0)
                    mx_rt  = max(mx_rt,  abs(u1-u0), abs(v1-v0))
                    ! magnitude preserved by the rotation
                    mag0 = sqrt(u0*u0 + v0*v0); mag1 = sqrt(ug*ug + vg*vg)
                    mx_mag = max(mx_mag, abs(mag1-mag0))
                    ! flag_coord=1 at the geographic coords must equal flag_coord=0 at the rotated coords
                    u2 = u0; v2 = v0
                    call vector_r2g(u2, v2, glon, glat, 1)
                    mx_flag = max(mx_flag, abs(u2-ug), abs(v2-vg))
                    if (check_identity) mx_rt = max(mx_rt, abs(ug-u0), abs(vg-v0))
                end do
            end do
        end do
    end subroutine sweep

end program test_vector_rotate
