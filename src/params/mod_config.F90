module mod_config
    ! Write-once run configuration (decision D4): module variables hold the run
    ! settings, populated once by read_config() from namelist.config at init and
    ! treated read-only thereafter. Per-entity config (e.g. per-tracer scheme
    ! strings) rides in the data types, NOT here.
    !
    ! Field names, defaults and namelist groups are transcribed from FESOM2 v2.7.3
    ! src/gen_modules_config.F90 (module g_config). Only the groups needed through
    ! M2 are included; more are added per milestone.
    use mod_precision, only: WP, MAX_PATH
    implicit none
    public
    save

    ! --- &modelname (gen_modules_config.F90:13) ---
    character(10)      :: runid = 'test1'
    namelist /modelname/ runid

    ! --- &timestep (gen_modules_config.F90:18-21) ---
    integer            :: step_per_day     = 72
    integer            :: run_length       = 1
    character          :: run_length_unit  = 'y'      ! y, d, s
    namelist /timestep/ step_per_day, run_length, run_length_unit

    ! --- &paths (gen_modules_config.F90:27-37) ---
    character(MAX_PATH) :: MeshPath        = './mesh/'
    character(MAX_PATH) :: ClimateDataPath = './hydrography/'
    character(MAX_PATH) :: ResultPath      = './result/'
    character(MAX_PATH) :: RestartInPath   = ''
    character(MAX_PATH) :: RestartOutPath  = ''
    character(20)       :: MeshId          = 'NONE'
    namelist /paths/ MeshPath, ClimateDataPath, ResultPath, MeshId, &
                     RestartInPath, RestartOutPath

    ! --- &ale_def (gen_modules_config.F90:80) ---
    character(20)      :: which_ALE          = 'linfs'   ! 'linfs','zlevel','zstar'
    logical            :: use_partial_cell   = .false.
    real(kind=WP)      :: partial_cell_thresh = 0.0_WP
    real(kind=WP)      :: min_hnode          = 0.5_WP
    integer            :: lzstar_lev         = 4
    real(kind=WP)      :: max_ice_loading    = 5.0_WP
    namelist /ale_def/ which_ALE, use_partial_cell, partial_cell_thresh, &
                       min_hnode, lzstar_lev, max_ice_loading

    ! --- &geometry (gen_modules_config.F90:108) ---
    logical            :: cartesian      = .false.
    logical            :: fplane         = .false.
    real(kind=WP)      :: cyclic_length  = 360.0_WP    ! [degree]
    logical            :: rotated_grid   = .true.
    logical            :: force_rotation = .true.
    real(kind=WP)      :: alphaEuler     = 50.0_WP     ! Euler angles [degree]: z, new-x, new-z
    real(kind=WP)      :: betaEuler      = 15.0_WP
    real(kind=WP)      :: gammaEuler     = -90.0_WP
    namelist /geometry/ cartesian, fplane, cyclic_length, rotated_grid, &
                        force_rotation, alphaEuler, betaEuler, gammaEuler

    ! --- &calendar (gen_modules_config.F90:121) ---
    logical            :: include_fleapyear = .false.
    logical            :: use_flpyrcheck    = .true.
    namelist /calendar/ include_fleapyear, use_flpyrcheck

    ! --- &configuration (subset; gen_modules_config.F90 configuration group) ---
    ! use_sw_pene: FESOM2 g_config default is .true., but the FESOM3 byte-gate drivers configure
    ! physics in-code (they do not read namelist.config), so the in-code default is .false.
    ! ("unset = off") — a fail-safe: the M5c sw_3d tracer term (diff_ver_part_impl_ale) and the
    ! KPP bldepth read this flag, and dyn%work%sw_3d is allocated only when sw_pene is enabled, so
    ! any PP/reduced driver that leaves it unset must see .false. The sw_pene drivers (M5c
    ! fesom_lifecycle_native, M5a fesom_pressuredump) set it .true. explicitly.
    logical            :: use_sw_pene = .false.
    logical            :: use_ice     = .false.
    namelist /configuration/ use_sw_pene, use_ice

    ! --- tracer time stepping (FESOM2 o_PARAM, oce_modules.F90:92) ---
    ! Adams-Bashforth(2) off-centring offset used by init_tracers_AB. Kept a
    ! (non-parameter) module variable to mirror FESOM2's runtime `epsilon` so the
    ! compiler cannot fold (1.5_WP+ab_epsilon) to a literal — see oce_tracer_mod.
    real(kind=WP)      :: ab_epsilon = 0.1_WP

    ! --- derived (not from namelist) ---
    real(kind=WP)      :: dt = 1200.0_WP    ! [s], = 86400/step_per_day after read_config

contains

    subroutine read_config(nml_path, ierr)
        ! Read namelist.config once. ierr=0 on success; >0 on parse error or
        ! missing file. Missing optional groups keep their defaults (warn).
        character(len=*), intent(in)  :: nml_path
        integer,          intent(out) :: ierr
        integer :: u, ios

        ierr = 0
        open(newunit=u, file=trim(nml_path), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(a)') 'read_config: cannot open '//trim(nml_path)
            ierr = 1; return
        end if

        call read_group_config(u, 'modelname',     ierr); if (ierr > 0) return
        call read_group_config(u, 'timestep',      ierr); if (ierr > 0) return
        call read_group_config(u, 'paths',         ierr); if (ierr > 0) return
        call read_group_config(u, 'ale_def',       ierr); if (ierr > 0) return
        call read_group_config(u, 'geometry',      ierr); if (ierr > 0) return
        call read_group_config(u, 'calendar',      ierr); if (ierr > 0) return
        call read_group_config(u, 'configuration', ierr); if (ierr > 0) return

        close(u)

        if (step_per_day <= 0) then
            write(*,'(a)') 'read_config: step_per_day must be > 0'
            ierr = 1; return
        end if
        dt = 86400.0_WP / real(step_per_day, WP)
    end subroutine read_config

    subroutine read_group_config(u, group, ierr)
        ! Rewind + read one namelist group. ios>0 (malformed) -> fatal (ierr=2);
        ! ios<0 (group absent) -> keep defaults. Order-independent.
        integer,          intent(in)  :: u
        character(len=*), intent(in)  :: group
        integer,          intent(out) :: ierr
        integer :: ios

        ierr = 0
        rewind(u)
        select case (group)
        case ('modelname');     read(u, nml=modelname,     iostat=ios)
        case ('timestep');      read(u, nml=timestep,      iostat=ios)
        case ('paths');         read(u, nml=paths,         iostat=ios)
        case ('ale_def');       read(u, nml=ale_def,       iostat=ios)
        case ('geometry');      read(u, nml=geometry,      iostat=ios)
        case ('calendar');      read(u, nml=calendar,      iostat=ios)
        case ('configuration'); read(u, nml=configuration, iostat=ios)
        case default;           ios = 0
        end select

        if (ios > 0) then
            write(*,'(a)') 'read_config: malformed namelist group &'//trim(group)
            ierr = 2
        end if
    end subroutine read_group_config

end module mod_config
