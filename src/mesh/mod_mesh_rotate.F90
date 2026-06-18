module mod_mesh_rotate
    ! Rotated-grid transforms (displaced pole). Transcribed from FESOM2 v2.7.3 /
    ! tracer_dwarf gen_modules_rotate_grid.F90. Euler convention (alpha,beta,gamma):
    ! rotate about z, then new x, then new z. Angles are passed in DEGREES and
    ! converted here (FESOM2 converts alphaEuler*rad at setup); cyclic_length too.
    use mod_precision, only: WP, MP
    use mod_constants, only: rad
    implicit none
    private
    public :: init_mesh_rotation, g2r, r2g, trim_cyclic

    real(kind=WP), save :: r2g_matrix(3,3) = 0.0_WP
    real(kind=WP), save :: cyclic_length_rad = 6.283185307179586_WP  ! 2*pi default

contains

    subroutine init_mesh_rotation(alpha_deg, beta_deg, gamma_deg, cyclic_deg)
        real(kind=WP), intent(in) :: alpha_deg, beta_deg, gamma_deg, cyclic_deg
        real(kind=WP) :: al, be, ga
        al = alpha_deg * rad; be = beta_deg * rad; ga = gamma_deg * rad
        r2g_matrix(1,1) = cos(ga)*cos(al) - sin(ga)*cos(be)*sin(al)
        r2g_matrix(1,2) = cos(ga)*sin(al) + sin(ga)*cos(be)*cos(al)
        r2g_matrix(1,3) = sin(ga)*sin(be)
        r2g_matrix(2,1) = -sin(ga)*cos(al) - cos(ga)*cos(be)*sin(al)
        r2g_matrix(2,2) = -sin(ga)*sin(al) + cos(ga)*cos(be)*cos(al)
        r2g_matrix(2,3) = cos(ga)*sin(be)
        r2g_matrix(3,1) = sin(be)*sin(al)
        r2g_matrix(3,2) = -sin(be)*cos(al)
        r2g_matrix(3,3) = cos(be)
        cyclic_length_rad = cyclic_deg * rad
    end subroutine init_mesh_rotation

    subroutine g2r(glon, glat, rlon, rlat)
        ! geographical -> rotated (all radians)
        real(kind=WP), intent(in)  :: glon, glat
        real(kind=WP), intent(out) :: rlon, rlat
        real(kind=WP) :: xr, yr, zr, xg, yg, zg
        xg = cos(glat)*cos(glon); yg = cos(glat)*sin(glon); zg = sin(glat)
        xr = r2g_matrix(1,1)*xg + r2g_matrix(1,2)*yg + r2g_matrix(1,3)*zg
        yr = r2g_matrix(2,1)*xg + r2g_matrix(2,2)*yg + r2g_matrix(2,3)*zg
        zr = r2g_matrix(3,1)*xg + r2g_matrix(3,2)*yg + r2g_matrix(3,3)*zg
        rlat = asin(zr)
        if (yr == 0.0_WP .and. xr == 0.0_WP) then
            rlon = 0.0_WP
        else
            rlon = atan2(yr, xr)
        end if
    end subroutine g2r

    subroutine r2g(glon, glat, rlon, rlat)
        ! rotated -> geographical (all radians)
        real(kind=WP), intent(out) :: glon, glat
        real(kind=WP), intent(in)  :: rlon, rlat
        real(kind=WP) :: xr, yr, zr, xg, yg, zg
        xr = cos(rlat)*cos(rlon); yr = cos(rlat)*sin(rlon); zr = sin(rlat)
        xg = r2g_matrix(1,1)*xr + r2g_matrix(2,1)*yr + r2g_matrix(3,1)*zr
        yg = r2g_matrix(1,2)*xr + r2g_matrix(2,2)*yr + r2g_matrix(3,2)*zr
        zg = r2g_matrix(1,3)*xr + r2g_matrix(2,3)*yr + r2g_matrix(3,3)*zr
        glat = asin(zg)
        if (yg == 0.0_WP .and. xg == 0.0_WP) then
            glon = 0.0_WP
        else
            glon = atan2(yg, xg)
        end if
    end subroutine r2g

    elemental subroutine trim_cyclic(b)
        ! wrap a longitude difference into (-L/2, L/2). cyclic_length in radians.
        real(kind=WP), intent(inout) :: b
        if (b >  0.5_WP*cyclic_length_rad) b = b - cyclic_length_rad
        if (b < -0.5_WP*cyclic_length_rad) b = b + cyclic_length_rad
    end subroutine trim_cyclic

end module mod_mesh_rotate
