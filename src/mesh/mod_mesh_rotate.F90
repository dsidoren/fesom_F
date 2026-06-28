module mod_mesh_rotate
    ! Rotated-grid transforms (displaced pole). Transcribed from FESOM2 v2.7.3 /
    ! tracer_dwarf gen_modules_rotate_grid.F90. Euler convention (alpha,beta,gamma):
    ! rotate about z, then new x, then new z. Angles are passed in DEGREES and
    ! converted here (FESOM2 converts alphaEuler*rad at setup); cyclic_length too.
    use mod_precision, only: WP, MP
    use mod_constants, only: rad
    implicit none
    private
    public :: init_mesh_rotation, g2r, r2g, vector_g2r, vector_r2g, trim_cyclic, get_cyclic_length

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

    subroutine vector_g2r(tlon, tlat, lon, lat, flag_coord)
        ! Rotate a 2-D vector (tlon,tlat) from geographic to rotated-mesh components.
        ! Transcribed VERBATIM from FESOM2 gen_modules_rotate_grid.F90:120-160. Used
        ! by the M2.10 forcing read to rotate the wind interpolation coefficients
        ! (flag_coord=0 => (lon,lat) are the ROTATED node coords, so r2g recovers the
        ! geographic (glon,glat); magnitude-preserving). All angles in radians.
        integer,       intent(in)    :: flag_coord
        real(kind=WP), intent(inout) :: tlon, tlat
        real(kind=WP), intent(in)    :: lon, lat
        real(kind=WP) :: rlon, rlat, glon, glat
        real(kind=WP) :: txg, tyg, tzg, txr, tyr, tzr
        if (flag_coord == 1) then  ! input is in geographical coordinates
           glon = lon; glat = lat
           call g2r(glon, glat, rlon, rlat)
        else                       ! input is in rotated coordinates
           rlon = lon; rlat = lat
           call r2g(glon, glat, rlon, rlat)
        end if
        ! vector in Cartesian geo. coordinates
        txg = -tlat*sin(glat)*cos(glon) - tlon*sin(glon)
        tyg = -tlat*sin(glat)*sin(glon) + tlon*cos(glon)
        tzg =  tlat*cos(glat)
        ! vector in rotated Cartesian coordinates
        txr = r2g_matrix(1,1)*txg + r2g_matrix(1,2)*tyg + r2g_matrix(1,3)*tzg
        tyr = r2g_matrix(2,1)*txg + r2g_matrix(2,2)*tyg + r2g_matrix(2,3)*tzg
        tzr = r2g_matrix(3,1)*txg + r2g_matrix(3,2)*tyg + r2g_matrix(3,3)*tzg
        ! vector in rotated coordinates
        tlat = -sin(rlat)*cos(rlon)*txr - sin(rlat)*sin(rlon)*tyr + cos(rlat)*tzr
        tlon = -sin(rlon)*txr + cos(rlon)*tyr
    end subroutine vector_g2r

    subroutine vector_r2g(tlon, tlat, lon, lat, flag_coord)
        ! Rotate a 2-D vector (tlon,tlat) from rotated-mesh to geographic components — the OUTPUT
        ! rotation that "unrotates" velocities/winds before writing. Transcribed VERBATIM from FESOM2
        ! gen_modules_rotate_grid.F90:164-202; the exact inverse of vector_g2r (build the Cartesian
        ! vector from the ROTATED angles, apply the TRANSPOSED r2g_matrix, project onto the GEOGRAPHIC
        ! angles). flag_coord=0 => (lon,lat) are the ROTATED coords of the vector position — what the
        ! output engine passes (io_meandata.F90:io_r2g feeds coord_nod2D with flag_coord=0). Magnitude-
        ! preserving. All angles in radians. NOTE the oracle's variable names: txg/tyg/tzg here hold the
        ! ROTATED Cartesian vector and txr/tyr/tzr the GEOGRAPHIC one (opposite of vector_g2r) — kept
        ! verbatim for a byte-faithful transcription.
        integer,       intent(in)    :: flag_coord
        real(kind=WP), intent(inout) :: tlon, tlat
        real(kind=WP), intent(in)    :: lon, lat
        real(kind=WP) :: rlon, rlat, glon, glat
        real(kind=WP) :: txg, tyg, tzg, txr, tyr, tzr
        if (flag_coord == 1) then  ! input is in geographical coordinates
           glon = lon; glat = lat
           call g2r(glon, glat, rlon, rlat)
        else                       ! input is in rotated coordinates
           rlon = lon; rlat = lat
           call r2g(glon, glat, rlon, rlat)
        end if
        ! vector in rotated Cartesian coordinates
        txg = -tlat*sin(rlat)*cos(rlon) - tlon*sin(rlon)
        tyg = -tlat*sin(rlat)*sin(rlon) + tlon*cos(rlon)
        tzg =  tlat*cos(rlat)
        ! vector in geographic Cartesian coordinates (TRANSPOSED r2g_matrix)
        txr = r2g_matrix(1,1)*txg + r2g_matrix(2,1)*tyg + r2g_matrix(3,1)*tzg
        tyr = r2g_matrix(1,2)*txg + r2g_matrix(2,2)*tyg + r2g_matrix(3,2)*tzg
        tzr = r2g_matrix(1,3)*txg + r2g_matrix(2,3)*tyg + r2g_matrix(3,3)*tzg
        ! vector in geographic coordinates
        tlat = -sin(glat)*cos(glon)*txr - sin(glat)*sin(glon)*tyr + cos(glat)*tzr
        tlon = -sin(glon)*txr + cos(glon)*tyr
    end subroutine vector_r2g

    elemental subroutine trim_cyclic(b)
        ! wrap a longitude difference into (-L/2, L/2). cyclic_length in radians.
        real(kind=WP), intent(inout) :: b
        if (b >  0.5_WP*cyclic_length_rad) b = b - cyclic_length_rad
        if (b < -0.5_WP*cyclic_length_rad) b = b + cyclic_length_rad
    end subroutine trim_cyclic

    pure real(kind=WP) function get_cyclic_length()
        ! cyclic_length (radians) for the elem_center/edge_center wraps, which use
        ! FESOM2's own >=/< comparisons against cyclic_length/2 (NOT trim_cyclic).
        get_cyclic_length = cyclic_length_rad
    end function get_cyclic_length

end module mod_mesh_rotate
