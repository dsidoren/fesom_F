module mod_dyn
    ! Ocean dynamics type (decision D5: EVOLVING; prognostic state split from
    ! work/scratch). Structure from tracer_dwarf MOD_DYN (== FESOM2 v2.7.3 MOD_DYN),
    ! trimmed to v1 scope: split-explicit SSH (se_*), backscatter/UKE, energy
    ! diagnostics (ke_*) and iceberg arrays are dropped. fer_uv/fer_w are declared
    ! (GM, M4) but not yet serialized.
    !
    ! Aux 3D WORK fields (density/N^2/Kv/Av/hpressure/sw_alpha-beta/GM slopes, D5)
    ! are added to t_dyn_work at M2.1 with verified FESOM2 names.
    use mod_precision, only: WP
    use, intrinsic :: iso_fortran_env, only: int32
    use mod_binary_arrays, only: write_bin_array, read_bin_array
    implicit none
    save

    type t_solverinfo
        integer       :: ident   = 1
        integer       :: maxiter = 2000
        integer       :: restart = 15
        integer       :: fillin  = 3
        integer       :: lutype  = 2
        real(kind=WP) :: droptol = 1.e-8_WP
        real(kind=WP) :: soltol  = 1.e-5_WP
        real(kind=WP), allocatable :: rr(:), zz(:), pp(:), App(:)
    contains
        procedure :: write_si => write_t_solverinfo
        procedure :: read_si  => read_t_solverinfo
    end type t_solverinfo

    type t_dyn_work
        real(kind=WP), allocatable, dimension(:,:,:) :: uvnode_rhs
        real(kind=WP), allocatable, dimension(:,:)   :: u_c, v_c
        ! Aux 3D fields recomputed each step from T/S by the EOS/pressure pass
        ! (M2.1 oce_pressure_bv); NOT serialized — they are diagnostics of the
        ! prognostic tracer state, restored by recomputation, not from the restart.
        ! FESOM2 o_ARRAYS names. density_ref = density_0 unless use_density_ref.
        real(kind=WP), allocatable, dimension(:,:)   :: density_m_rho0  ! (nl-1,nod2D) in-situ density - density_ref (PGF)
        real(kind=WP), allocatable, dimension(:,:)   :: density_ref     ! (nl-1,nod2D) reference density
        real(kind=WP), allocatable, dimension(:,:)   :: hpressure       ! (nl,  nod2D) hydrostatic pressure
        real(kind=WP), allocatable, dimension(:,:)   :: bvfreq          ! (nl,  nod2D) N^2 (squared Brunt-Vaisala)
        real(kind=WP), allocatable, dimension(:,:)   :: pgf_x, pgf_y    ! (nl-1,elem2D) PGF (M2.2 oce_pgf), from hpressure
        ! M2.8 PP vertical mixing coefficients (oce_ale_mixing_pp). FESOM2 o_ARRAYS
        ! names Kv/Av; recomputed each step from N^2 + uvnode shear, NOT serialized.
        real(kind=WP), allocatable, dimension(:,:)   :: Kv              ! (nl,  nod2D) vertical diffusivity (tracers)
        real(kind=WP), allocatable, dimension(:,:)   :: Av              ! (nl,  elem2D) vertical viscosity (momentum)
    contains
        procedure :: write_dw => write_t_dyn_work
        procedure :: read_dw  => read_t_dyn_work
    end type t_dyn_work

    type t_dyn
        ! ---- prognostic / state ----
        real(kind=WP), allocatable, dimension(:,:,:)   :: uv, uv_rhs, fer_uv  ! (2,nl-1,elem2D)
        real(kind=WP), allocatable, dimension(:,:,:,:) :: uv_rhsAB            ! (AB_order-1,2,nl-1,elem2D) — FESOM2 order (ab_lvl,comp,nz,elem)
        real(kind=WP), allocatable, dimension(:,:,:)   :: uvnode              ! (2,nl-1,nod2D)
        real(kind=WP), allocatable, dimension(:,:)     :: w, w_e, w_i, w_old, cfl_z, fer_w ! (nl,nod2D)
        real(kind=WP), allocatable, dimension(:)       :: eta_n, d_eta, ssh_rhs, ssh_rhs_old ! (nod2D)
        integer :: AB_order = 2

        ! ---- run config (rides in the type, D6) ----
        logical       :: check_opt_visc = .true.
        integer       :: opt_visc       = 5
        real(kind=WP) :: visc_gamma0    = 0.03_WP
        real(kind=WP) :: visc_gamma1    = 0.1_WP
        real(kind=WP) :: visc_gamma2    = 0.285_WP
        ! harmonic (Laplacian) viscosity coefficients for opt_visc=7 (FESOM2 MOD_DYN:
        ! visc_gamma{0,1}_h). Default 0 -> pure biharmonic (the pi/reduced-M2 config).
        real(kind=WP) :: visc_gamma0_h  = 0.0_WP
        real(kind=WP) :: visc_gamma1_h  = 0.0_WP
        logical       :: use_ivertvisc  = .true.
        integer       :: momadv_opt     = 2
        logical       :: use_freeslip   = .false.
        logical       :: use_wsplit     = .false.
        real(kind=WP) :: wsplit_maxcfl  = 1.0_WP

        ! ---- sub-objects ----
        type(t_solverinfo) :: solverinfo
        type(t_dyn_work)   :: work
    contains
        procedure :: write_unformatted => write_t_dyn
        procedure :: read_unformatted  => read_t_dyn
        generic   :: write(unformatted) => write_unformatted
        generic   :: read(unformatted)  => read_unformatted
    end type t_dyn

contains

    subroutine write_t_solverinfo(si, unit)
        class(t_solverinfo), intent(in) :: si
        integer,             intent(in) :: unit
        integer :: iostat
        character(len=1024) :: iomsg
        write(unit, iostat=iostat, iomsg=iomsg) si%ident, si%maxiter, si%restart, &
            si%fillin, si%lutype
        write(unit, iostat=iostat, iomsg=iomsg) si%droptol, si%soltol
        call write_bin_array(si%rr,  unit, iostat, iomsg)
        call write_bin_array(si%zz,  unit, iostat, iomsg)
        call write_bin_array(si%pp,  unit, iostat, iomsg)
        call write_bin_array(si%App, unit, iostat, iomsg)
    end subroutine write_t_solverinfo

    subroutine read_t_solverinfo(si, unit)
        class(t_solverinfo), intent(inout) :: si
        integer,             intent(in)    :: unit
        integer :: iostat
        character(len=1024) :: iomsg
        read(unit, iostat=iostat, iomsg=iomsg) si%ident, si%maxiter, si%restart, &
            si%fillin, si%lutype
        read(unit, iostat=iostat, iomsg=iomsg) si%droptol, si%soltol
        call read_bin_array(si%rr,  unit, iostat, iomsg)
        call read_bin_array(si%zz,  unit, iostat, iomsg)
        call read_bin_array(si%pp,  unit, iostat, iomsg)
        call read_bin_array(si%App, unit, iostat, iomsg)
    end subroutine read_t_solverinfo

    subroutine write_t_dyn_work(dw, unit)
        class(t_dyn_work), intent(in) :: dw
        integer,           intent(in) :: unit
        integer :: iostat
        character(len=1024) :: iomsg
        call write_bin_array(dw%uvnode_rhs, unit, iostat, iomsg)
        call write_bin_array(dw%u_c,        unit, iostat, iomsg)
        call write_bin_array(dw%v_c,        unit, iostat, iomsg)
    end subroutine write_t_dyn_work

    subroutine read_t_dyn_work(dw, unit)
        class(t_dyn_work), intent(inout) :: dw
        integer,           intent(in)    :: unit
        integer :: iostat
        character(len=1024) :: iomsg
        call read_bin_array(dw%uvnode_rhs, unit, iostat, iomsg)
        call read_bin_array(dw%u_c,        unit, iostat, iomsg)
        call read_bin_array(dw%v_c,        unit, iostat, iomsg)
    end subroutine read_t_dyn_work

    subroutine write_t_dyn(dyn, unit, iostat, iomsg)
        class(t_dyn), intent(in)    :: dyn
        integer,      intent(in)    :: unit
        integer,      intent(out)   :: iostat
        character(*), intent(inout) :: iomsg
        write(unit, iostat=iostat, iomsg=iomsg) dyn%AB_order, dyn%opt_visc, &
            dyn%momadv_opt
        write(unit, iostat=iostat, iomsg=iomsg) dyn%check_opt_visc, &
            dyn%use_ivertvisc, dyn%use_freeslip, dyn%use_wsplit
        write(unit, iostat=iostat, iomsg=iomsg) dyn%visc_gamma0, dyn%visc_gamma1, &
            dyn%visc_gamma2, dyn%visc_gamma0_h, dyn%visc_gamma1_h, dyn%wsplit_maxcfl
        call dyn%solverinfo%write_si(unit)
        call dyn%work%write_dw(unit)
        call write_bin_array(dyn%uv,          unit, iostat, iomsg)
        call write_bin_array(dyn%uv_rhs,      unit, iostat, iomsg)
        call write_bin_array(dyn%uv_rhsAB,    unit, iostat, iomsg)
        call write_bin_array(dyn%uvnode,      unit, iostat, iomsg)
        call write_bin_array(dyn%w,           unit, iostat, iomsg)
        call write_bin_array(dyn%w_e,         unit, iostat, iomsg)
        call write_bin_array(dyn%w_i,         unit, iostat, iomsg)
        call write_bin_array(dyn%w_old,       unit, iostat, iomsg)
        call write_bin_array(dyn%cfl_z,       unit, iostat, iomsg)
        call write_bin_array(dyn%eta_n,       unit, iostat, iomsg)
        call write_bin_array(dyn%d_eta,       unit, iostat, iomsg)
        call write_bin_array(dyn%ssh_rhs,     unit, iostat, iomsg)
        call write_bin_array(dyn%ssh_rhs_old, unit, iostat, iomsg)
    end subroutine write_t_dyn

    subroutine read_t_dyn(dyn, unit, iostat, iomsg)
        class(t_dyn), intent(inout) :: dyn
        integer,      intent(in)    :: unit
        integer,      intent(out)   :: iostat
        character(*), intent(inout) :: iomsg
        read(unit, iostat=iostat, iomsg=iomsg) dyn%AB_order, dyn%opt_visc, &
            dyn%momadv_opt
        read(unit, iostat=iostat, iomsg=iomsg) dyn%check_opt_visc, &
            dyn%use_ivertvisc, dyn%use_freeslip, dyn%use_wsplit
        read(unit, iostat=iostat, iomsg=iomsg) dyn%visc_gamma0, dyn%visc_gamma1, &
            dyn%visc_gamma2, dyn%visc_gamma0_h, dyn%visc_gamma1_h, dyn%wsplit_maxcfl
        call dyn%solverinfo%read_si(unit)
        call dyn%work%read_dw(unit)
        call read_bin_array(dyn%uv,          unit, iostat, iomsg)
        call read_bin_array(dyn%uv_rhs,      unit, iostat, iomsg)
        call read_bin_array(dyn%uv_rhsAB,    unit, iostat, iomsg)
        call read_bin_array(dyn%uvnode,      unit, iostat, iomsg)
        call read_bin_array(dyn%w,           unit, iostat, iomsg)
        call read_bin_array(dyn%w_e,         unit, iostat, iomsg)
        call read_bin_array(dyn%w_i,         unit, iostat, iomsg)
        call read_bin_array(dyn%w_old,       unit, iostat, iomsg)
        call read_bin_array(dyn%cfl_z,       unit, iostat, iomsg)
        call read_bin_array(dyn%eta_n,       unit, iostat, iomsg)
        call read_bin_array(dyn%d_eta,       unit, iostat, iomsg)
        call read_bin_array(dyn%ssh_rhs,     unit, iostat, iomsg)
        call read_bin_array(dyn%ssh_rhs_old, unit, iostat, iomsg)
    end subroutine read_t_dyn

end module mod_dyn
