module mod_ice
    ! Sea-ice type (decision D5: EVOLVING; built fresh — the dwarf had deleted
    ! MOD_ICE). Minimal v1 skeleton: the EVP stress tensor (sigma11/12/22) is
    ! PROGNOSTIC (elastic memory persists across steps); concentration/mass/
    ! velocity are prognostic; rhs/strain are work. Fleshed out from FESOM2
    ! MOD_ICE at M3 (EVP).
    use mod_precision, only: WP
    use mod_binary_arrays, only: write_bin_array, read_bin_array
    implicit none
    save

    type t_ice_work
        real(kind=WP), allocatable, dimension(:) :: rhs_a, rhs_m, rhs_ms  ! advection rhs (nod2D)
        real(kind=WP), allocatable, dimension(:) :: eps11, eps12, eps22   ! strain rate (elem2D)
    contains
        procedure :: write_iw => write_t_ice_work
        procedure :: read_iw  => read_t_ice_work
    end type t_ice_work

    type t_ice
        ! ---- prognostic ----
        real(kind=WP), allocatable, dimension(:) :: a_ice, m_ice, m_snow  ! conc/ice mass/snow mass (nod2D)
        real(kind=WP), allocatable, dimension(:) :: uice, vice            ! ice velocity (nod2D)
        real(kind=WP), allocatable, dimension(:) :: sigma11, sigma12, sigma22 ! EVP stress (elem2D)
        real(kind=WP) :: ice_dt = 0.0_WP
        type(t_ice_work) :: work
    contains
        procedure :: write_unformatted => write_t_ice
        procedure :: read_unformatted  => read_t_ice
        generic   :: write(unformatted) => write_unformatted
        generic   :: read(unformatted)  => read_unformatted
    end type t_ice

contains

    subroutine write_t_ice_work(iw, unit)
        class(t_ice_work), intent(in) :: iw
        integer,           intent(in) :: unit
        integer :: iostat
        character(len=1024) :: iomsg
        call write_bin_array(iw%rhs_a,  unit, iostat, iomsg)
        call write_bin_array(iw%rhs_m,  unit, iostat, iomsg)
        call write_bin_array(iw%rhs_ms, unit, iostat, iomsg)
        call write_bin_array(iw%eps11,  unit, iostat, iomsg)
        call write_bin_array(iw%eps12,  unit, iostat, iomsg)
        call write_bin_array(iw%eps22,  unit, iostat, iomsg)
    end subroutine write_t_ice_work

    subroutine read_t_ice_work(iw, unit)
        class(t_ice_work), intent(inout) :: iw
        integer,           intent(in)    :: unit
        integer :: iostat
        character(len=1024) :: iomsg
        call read_bin_array(iw%rhs_a,  unit, iostat, iomsg)
        call read_bin_array(iw%rhs_m,  unit, iostat, iomsg)
        call read_bin_array(iw%rhs_ms, unit, iostat, iomsg)
        call read_bin_array(iw%eps11,  unit, iostat, iomsg)
        call read_bin_array(iw%eps12,  unit, iostat, iomsg)
        call read_bin_array(iw%eps22,  unit, iostat, iomsg)
    end subroutine read_t_ice_work

    subroutine write_t_ice(ice, unit, iostat, iomsg)
        class(t_ice), intent(in)    :: ice
        integer,      intent(in)    :: unit
        integer,      intent(out)   :: iostat
        character(*), intent(inout) :: iomsg
        write(unit, iostat=iostat, iomsg=iomsg) ice%ice_dt
        call write_bin_array(ice%a_ice,   unit, iostat, iomsg)
        call write_bin_array(ice%m_ice,   unit, iostat, iomsg)
        call write_bin_array(ice%m_snow,  unit, iostat, iomsg)
        call write_bin_array(ice%uice,    unit, iostat, iomsg)
        call write_bin_array(ice%vice,    unit, iostat, iomsg)
        call write_bin_array(ice%sigma11, unit, iostat, iomsg)
        call write_bin_array(ice%sigma12, unit, iostat, iomsg)
        call write_bin_array(ice%sigma22, unit, iostat, iomsg)
        call ice%work%write_iw(unit)
    end subroutine write_t_ice

    subroutine read_t_ice(ice, unit, iostat, iomsg)
        class(t_ice), intent(inout) :: ice
        integer,      intent(in)    :: unit
        integer,      intent(out)   :: iostat
        character(*), intent(inout) :: iomsg
        read(unit, iostat=iostat, iomsg=iomsg) ice%ice_dt
        call read_bin_array(ice%a_ice,   unit, iostat, iomsg)
        call read_bin_array(ice%m_ice,   unit, iostat, iomsg)
        call read_bin_array(ice%m_snow,  unit, iostat, iomsg)
        call read_bin_array(ice%uice,    unit, iostat, iomsg)
        call read_bin_array(ice%vice,    unit, iostat, iomsg)
        call read_bin_array(ice%sigma11, unit, iostat, iomsg)
        call read_bin_array(ice%sigma12, unit, iostat, iomsg)
        call read_bin_array(ice%sigma22, unit, iostat, iomsg)
        call ice%work%read_iw(unit)
    end subroutine read_t_ice

end module mod_ice
