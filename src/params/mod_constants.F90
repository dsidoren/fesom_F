module mod_constants
    ! Physical + mathematical constants, transcribed verbatim from FESOM2 v2.7.3
    ! src/oce_modules.F90  MODULE o_PARAM  (line cites per value below).
    !
    ! Precision split per D3: geometry/coordinate/math constants are MP (mesh
    ! precision); physical quantities entering the WP arithmetic are WP. At the
    ! double-precision anchor MP == WP == 8, so every value below is bit-identical
    ! to FESOM2 (which declares them all `_WP` with WP=8). The `_MP` tags exist
    ! only so a future single/half WP build keeps geometry at >= single.
    use mod_precision, only: WP, MP
    implicit none
    public

    ! --- mathematical / geometric (MP) ---
    real(kind=MP), parameter :: pi      = 3.14159265358979_MP          ! oce_modules.F90:11 (truncated pi)
    real(kind=MP), parameter :: rad     = pi / 180.0_MP                ! oce_modules.F90:12
    real(kind=MP), parameter :: r_earth = 6367500.0_MP                 ! oce_modules.F90:16 [m]
    real(kind=MP), parameter :: omega   = 2 * pi / (3600.0_MP * 24.0_MP) ! oce_modules.F90:17 [rad/s]

    ! --- physical (WP) ---
    real(kind=WP), parameter :: density_0   = 1030.0_WP               ! oce_modules.F90:13 [kg/m^3]
    real(kind=WP), parameter :: density_0_r = 1.0_WP / density_0      ! oce_modules.F90:14 [m^3/kg]
    real(kind=WP), parameter :: g           = 9.81_WP                 ! oce_modules.F90:15 [m/s^2]
    real(kind=WP), parameter :: vcpw        = 4.2e6_WP                ! oce_modules.F90:18 [J/m^3/K] water heat capacity
    real(kind=WP), parameter :: inv_vcpw    = 1.0_WP / vcpw           ! oce_modules.F90:19
    real(kind=WP), parameter :: small       = 1.0e-8_WP              ! oce_modules.F90:20

end module mod_constants
