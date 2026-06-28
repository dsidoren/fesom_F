! ============ FESOM3 output configuration (M9 Zarr output) ====================================
! Read from the rundir at runtime (mirroring the FESOM2 oracle namelists); set FESOM3_OUTPUT=<dir>
! to enable output and FESOM3_NAMELIST_IO=<path> to point elsewhere (default ./namelist.io). The
! per-knob FESOM3_* environment variables OVERRIDE the values here. The FESOM2 analog is
! config/namelist.io in the v2.7.3 tree (&nml_general + &nml_list).
!
! Output is emitted as hand-rolled Zarr v2 stores, one per variable per period:
!   <out_dir>/<name>.fesom.<YYYY>.zarr            (filesplit_freq='y')
!   <out_dir>/<name>.fesom.<YYYY>_<MM>.zarr       (filesplit_freq='m')
! xarray-readable, UGRID mesh in fesom.mesh.diag.zarr.

&nml_general
  n_writers      = 0           ! writer-rank subset (0 => one writer per chunk-block; FESOM3_N_WRITERS)
  chunk_time     = 1           ! time chunk (1 => append-only; >1 => read-modify-write) FESOM3_CHUNK_TIME
  chunk_vert     = 0           ! vertical chunk (0 => full depth single chunk)           FESOM3_CHUNK_VERT
  chunk_horiz    = 500000      ! horizontal (node/elem) chunk size                       FESOM3_CHUNK_HORIZ
  compressor     = 'none'      ! 'none' | 'lz4'   (data-array codec)                     FESOM3_COMPRESSOR
  filesplit_freq = 'y'         ! 'y' | 'm'  -> per-year / per-month stores               FESOM3_FILESPLIT
  ! vector frame for velocity/wind pairs (unod/vnod, ...): 'geographic' r2g-rotates to true
  ! east/north (= FESOM2 vec_autorotate=.true., the production default); 'native' writes the raw
  ! rotated-mesh components (= vec_autorotate=.false.).                       (FESOM3_VEC_FRAME)
  vec_frame      = 'geographic'
/

! OUTPUT VARIABLE LIST — rows: '<name>', <freq>, '<unit y|m|d|h|s>', <precision 4|8>, '<mean|snap>'
! Each row registers a stream; the per-field cadence (freq+unit) drives the FESOM2-ported events
! (annual/monthly/daily/hourly/step). Known names: ssh/sst/sss/a_ice/m_ice/m_snow (2-D node),
! temp/salt/w (3-D node), unod/vnod (node velocity vector pair, r2g-rotated per vec_frame).
&nml_list
  io_list = 'ssh   ', 1, 'm', 4, 'snap',
            'sst   ', 1, 'm', 4, 'snap',
            'sss   ', 1, 'm', 4, 'snap',
            'temp  ', 1, 'm', 4, 'snap',
            'salt  ', 1, 'm', 4, 'snap',
            'unod  ', 1, 'm', 4, 'mean',     ! zonal velocity at nodes  [m/s] (r2g-rotated; vec pair)
            'vnod  ', 1, 'm', 4, 'mean',     ! meridional velocity at nodes [m/s]
            'u     ', 1, 'm', 4, 'mean',     ! zonal velocity at elements [m/s] (r2g-rotated; vec pair)
            'v     ', 1, 'm', 4, 'mean',     ! meridional velocity at elements [m/s]
            'Av    ', 1, 'm', 4, 'mean',     ! vertical viscosity at elements (full levels nz) [m2/s]
/
! Element GM bolus (uncomment when Fer_GM is on): 'bolus_u'/'bolus_v' (dyn%fer_uv, vector pair).
