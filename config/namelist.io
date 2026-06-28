! ============ FESOM3 output configuration (M9 Zarr output) ====================================
! TEMPLATE / SCHEMA. As of Task 2.5 the output knobs are driven by FESOM3_* environment variables
! (listed per knob below); the &nml_general / &nml_list NAMELIST PARSING lands in Task 2.6, which
! will read this file from the rundir (mirroring the FESOM2 oracle namelists). The FESOM2 analog is
! config/namelist.io in the v2.7.3 tree (&diag_list + &nml_list).
!
! Output is emitted as hand-rolled Zarr v2 stores, one per variable per year:
!   <out_dir>/<name>.fesom.<YYYY>.zarr     (xarray-readable, UGRID mesh in fesom.mesh.diag.zarr)

&nml_general
  n_writers      = 0           ! writer-rank subset (0 => one writer per chunk-block; FESOM3_N_WRITERS)
  chunk_time     = 1           ! time chunk (1 => append-only, no read-modify-write; >1 is Task 2.6)
  chunk_vert     = 0           ! vertical chunk (0 => full depth single chunk; Task 2.6)
  chunk_horiz    = 500000      ! horizontal (node/elem) chunk size            (FESOM3_CHUNK_HORIZ)
  compressor     = 'none'      ! 'none' | 'lz4'                               (Task 2.6)
  filesplit_freq = 'y'         ! 'y' | 'm'  -> per-year / per-month stores    (Task 2.6)
  ! vector frame for velocity/wind pairs (unod/vnod, ...): 'geographic' r2g-rotates to true
  ! east/north (= FESOM2 vec_autorotate=.true., the production default); 'native' writes the raw
  ! rotated-mesh components (= vec_autorotate=.false.).                       (FESOM3_VEC_FRAME)
  vec_frame      = 'geographic'
/

! OUTPUT VARIABLE LIST — rows: '<name>', <freq>, '<unit y|m|d|h|s>', <precision 4|8>, '<mean|snap>'
! (Task 2.6 parses this; the lifecycle currently registers a fixed snapshot set + unod/vnod via env.)
&nml_list
  io_list = 'ssh   ', 1, 'm', 4, 'snap',
            'sst   ', 1, 'm', 4, 'snap',
            'sss   ', 1, 'm', 4, 'snap',
            'temp  ', 1, 'm', 4, 'snap',
            'salt  ', 1, 'm', 4, 'snap',
            'unod  ', 1, 'm', 4, 'mean',     ! zonal velocity at nodes  [m/s] (r2g-rotated; vec pair)
            'vnod  ', 1, 'm', 4, 'mean',     ! meridional velocity at nodes [m/s]
/
