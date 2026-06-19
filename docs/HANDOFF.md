# FESOM3 — Handoff (durable state across sessions)

Single source of truth for "where are we / what's next". Update at the end of every task.

## Where we are

- **Milestone:** M2 (minimal ocean dynamical core) — **IN PROGRESS.** **M2.1 `pressure_bv` +
  M2.2 hydrostatic PGF + M2.3 vel_rhs + M2.4 momentum advection (FULL `UV_rhs`) + M2.4 biharmonic
  viscosity + M2.5 implicit vertical viscosity (TDMA) COMPLETE ✓** (M2.1: EOS split-form
  `density_m_rho0` + top-down `hpressure` + N²/`bvfreq` raw+smoothed; M2.2: `gradient_sca`-contraction
  `pgf_x`/`pgf_y`; M2.3: `coriolis` geometry field + `compute_vel_rhs` Coriolis+AB2+PGF+SSH-gradient;
  M2.4: `momentum_adv_scalar` (w·du/dz + u·du/dx) ADDED into the same `UV_rhsAB(1,1:2,·)` slot → the
  FULL `UV_rhs`, THEN `viscosity_filter(7)` `visc_filt_bidiff` biharmonic viscosity as a SEPARATE
  operator on the post-`compute_vel_rhs` `UV_rhs`; M2.5: `impl_vert_visc_ale` per-element tridiagonal
  Thomas solve [implicit `Av` + vertical advection `w_i` + wind-stress/bottom-drag BCs] OVERWRITING
  `UV_rhs`, a SEPARATE operator after `viscosity_filter`; all `max|Δ|=0` vs FESOM2 on pi 1-rank; one
  gate `tools/run_pressure_gate.sh`, **28 fields**). **Next: M2.6 SSH — stiffness matrix + `ssh_rhs` +
  CG solve.** M1 (tracer advection) — **COMPLETE ✓ at the 1-rank anchor** (tag `m1`).
  M1.1–M1.4 byte-gates `max|Δ|=0` vs FESOM2 on pi 1-rank (geometry + horiz/vert advection +
  FCT limiter + the assembled driver/step). **M1.5 (multi-rank) is FOLDED INTO M2.12**
  (decision 2026-06-19): a multi-rank advection byte-gate needs the local-mesh remap
  (global→local numbering/connectivity/geometry), which `read_mesh` does not build yet
  (`mod_mesh_read.F90:28` errors at `npes/=1`); M2.12 builds it for the whole-model
  multi-rank byte-match, which subsumes advection. The 1-rank anchor IS the D0/D7 bit-identity
  gold standard ("serial == 1-rank MPI"). **Next: M2** (dynamics; reduced-M2 oracle namelist).
  M0 (Foundation) — **COMPLETE ✓** (tag `m0`); 13/13 ctest green Intel+GNU dp, debug clean.
- **M1 detail (all CLOSED ✓ on pi 1-rank).** **Geometry byte-gate** — FESOM3 mesh geometry is
  `max|Δ|=0` vs FESOM2 on pi (1-rank): elem_area, elem_cos, metric_factor, gradient_sca,
  edge_dxdy, edge_cross_dxdy, area/areasvol(+inv), coord_nod2D, elem2D_nodes, edges,
  edge_tri, all level arrays. Run it: `tools/run_geom_gate.sh`.
  ⚠️ **Caveat:** pi has 0/5839 CW swaps, so the `enforce_cw_orientation` vertex-reorder
  path is byte-identical-by-construction to FESOM2 `test_tri` but NOT yet empirically
  gated (soufflet=228/5700, CORE2 will have swaps). Confirm at the M2.11 CORE2 gate or on
  soufflet. See LESSONS L8.
- **M1.1 horizontal tracer advection byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs
  FESOM2 on every field: `del_ttf_advhoriz` (the gate target) AND `adv_flux_hor` for BOTH
  upwind (UPW1) and MUSCL, plus every intermediate (helem, nboundary_lay, edge_up_dn_tri,
  tr_xy, edge_up_dn_grad). Run it: `tools/run_advhor_gate.sh`. Built: `src/oce/`
  oce_tracer_grad (tr_xy), oce_muscl_adv (nboundary_lay/edge_up_dn_tri/edge_up_dn_grad),
  oce_adv_tra_hor (upw1/muscl/mfct), oce_adv_tra_flux (scatter).
- **M1.2 vertical tracer advection byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs
  FESOM2 on `del_ttf_advvert` (gate target) AND `adv_flux_ver` for BOTH upwind (UPW1) and
  QR4C, plus the vertical geometry it consumes (`zbar_3d_n`, `Z_3d_n`, `area`) and the
  prescribed `wvel`. Built `src/oce/oce_adv_tra_ver.F90` (adv_tra_ver_upw1, adv_tra_ver_qr4c).
  Same `tools/run_advhor_gate.sh` (extended). PASSED first run; debug `-check all` clean.
- **M1.3 FCT (Zalesak) limiter byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs FESOM2
  on the whole FCT path: `fct_LO` (low-order solution), the clipped `adv_flux_{hor,ver}_fct`,
  the limiter internals (`fct_ttf_max/min` bounds, `fct_plus/minus` factors, pre-clip
  `*_fct_ho` fluxes), the final `del_ttf_{advhoriz,advvert}_fct`, plus `hnode/hnode_new` and
  the standalone `adv_flux_hor_mfct`. Config = pi namelist (MFCT opth=0.0 / QR4C optv=1.0 /
  FCT). Built `src/oce/oce_adv_tra_fct.F90` (`oce_tra_adv_fct`); the `use_lo` branch of
  `oce_adv_tra_flux.F90` is now LIVE. Same `tools/run_advhor_gate.sh` (extended). PASSED first
  run; debug `-check all` clean. The limiter actively clipped (~41% of nodes, fct_plus/minus
  min=0; ~14% of fluxes changed), so b1/b2/b3 were genuinely exercised. See LESSONS L11.
- **M1.4 assembled tracer-advection step byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs
  FESOM2 on the REAL driver: `valuesAB` (init_tracers_AB AB(2) interpolation), the per-tracer
  `del_ttf_advhoriz/advvert/del_ttf` for the FCT path (do_oce_adv_tra dispatch + the
  adv_tracers_ale `del_ttf += advhoriz+advvert` accumulation) AND a non-FCT config
  (MUSCL/QR4C/NON, ph=pv=0.75 → exercises the `do_zero_flux=.true.` path + per-tracer order
  knobs). Built `src/oce/oce_adv_tra_driver.F90` (`do_oce_adv_tra`), `src/oce/oce_tracer_mod.F90`
  (`init_tracers_AB`, AB offset `ab_epsilon=0.1` in mod_config), `src/oce/oce_ale_tracer.F90`
  (`adv_tracers_ale`/`advect_tracer`). The gate now drives FESOM2's REAL `init_tracers_AB` +
  `do_oce_adv_tra` (oracle shim `fesom_advhor_dump.F90` extended; `libfesom.so` rebuilt) — not
  the inline orchestration — so the AB-interpolation is oracle-gated, not self-checked. Same
  `tools/run_advhor_gate.sh` (41 fields). PASSED first run; debug `-check all` clean. valuesAB
  exactly = 1.6·values−0.6·valuesold (no L7-fold), FCT vs non-FCT del_ttf differ by max 26.
  See LESSONS L12.
  **M1 COMPLETE at the 1-rank anchor → tag `m1`.** M1.5 (multi-rank) folded into M2.12 (above).
- **Done:** M0.1 ✓ build. M0.2 ✓ params/. M0.3 ✓ types/. M0.4 ✓ mod_partitioning
  (par_init/par_ex/set_partition; dist_<NP>/ reader transcribed from oce_mesh.F90;
  1-rank synthesis D7). test_partit passes 1/2/8-rank, Intel+GNU dp.
  M0.5 ✓ mod_halo (generic exchange_nod/exchange_elem, broadcast-only manual pack,
  stale-halo probe; test_halo 1/2/8-rank Intel+GNU dp).
  M0.6 ✓ mod_dump (gid-keyed node+elem dump, byte-format-identical to FESOM2
  fesom_dump_shim.F90) + tools/dump_diff.py (first divergent substep, --selftest);
  test_dump + dump_diff_selftest pass.
  M0.7 ✓ mesh/ (rotate, read, areas, analytic; self-consistency gated).
  M0.8 ✓ step/mod_model (t_model + model_init/step/finalize) + drivers/fesom_analytic;
  runs end-to-end 1+2 ranks; debug build clean.
  M0.7-geom ✓ (closed in M1) mesh geometry byte-matches FESOM2 on pi 1-rank.
  M1.1 ✓ horizontal tracer advection (upw1 + MUSCL): operator-diff `max|Δ|=0` on
  `del_ttf_advhoriz` + `adv_flux_hor` + all intermediates vs FESOM2 on pi 1-rank.
  M1.2 ✓ vertical tracer advection (upw1 + QR4C): operator-diff `max|Δ|=0` on
  `del_ttf_advvert` (oce_ale_tracer.F90:240) + `adv_flux_ver` + vertical geometry
  (zbar_3d_n/Z_3d_n/area) vs FESOM2 on pi 1-rank.
  M1.3 ✓ FCT (Zalesak) limiter (oce_adv_tra_fct.F90 → oce_tra_adv_fct): operator-diff
  `max|Δ|=0` on fct_LO + clipped adv_flux_{hor,ver}_fct + fct_ttf_max/min + fct_plus/minus
  + del_ttf_*_fct + hnode/hnode_new + standalone MFCT vs FESOM2 on pi 1-rank. Activated the
  `use_lo` hnode/hnode_new branch of oce_adv_tra_flux. Limiter clipped ~41% of nodes.
  M1.4 ✓ assembled advection step (do_oce_adv_tra + init_tracers_AB + adv_tracers_ale):
  operator-diff `max|Δ|=0` on valuesAB + del_ttf_{advhoriz,advvert,}_step (FCT) + the
  non-FCT (do_zero_flux) del_ttf vs FESOM2's REAL driver on pi 1-rank. Oracle now runs
  FESOM2's own init_tracers_AB + do_oce_adv_tra (shim extended, libfesom.so rebuilt).
- **M2.1 `pressure_bv` byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs FESOM2 on
  `density_m_rho0` (EOS PGF density anomaly), `hpressure` (top-down hydrostatic integration) AND
  `bvfreq` BOTH raw (pre-smoothing) and smoothed (the horizontal `smooth_nod` mass-matrix sweep,
  N2smth_hidx=1), plus every input (temp, salt, density_ref, zbar_3d_n, Z_3d_n, hnode). Built
  `src/oce/oce_pressure_bv.F90` (`pressure_bv` + `densityJM_components` split EOS + `smooth_nod`);
  added density_m_rho0/density_ref/hpressure/bvfreq to `t_dyn_work` (mod_dyn; recomputed each step,
  not serialized). Gate `tools/run_pressure_gate.sh` (FESOM2 1-rank shim + FESOM3 driver +
  `pressure_diff.py`). PASSED first run; smoother changed 100% of valid entries (non-vacuous);
  Debug `-check all` clean; M1 advhor gate + 13/13 ctest still green. See LESSONS L13.
- **M2.2 hydrostatic PGF byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs FESOM2 on `pgf_x`/`pgf_y`,
  the element pressure-gradient force = `Σ_k gradient_sca(k,elem)·hpressure(nz,elnodes_k)/density_0`
  (the M2.1 `hpressure` contracted with the geometry-gated `gradient_sca`, same shape as M1.1
  `tracer_gradient_elements`). Built `src/oce/oce_pgf.F90` (`pressure_force_4_linfs_fullcell`); added
  `pgf_x`/`pgf_y` to `t_dyn_work` (mod_dyn; recomputed each step, not serialized). Same gate
  `tools/run_pressure_gate.sh` (now 12 fields; extended the FESOM2 shim to call the REAL
  `pressure_force_4_linfs_fullcell` + the FESOM3 driver to run `oce_pgf`). PASSED first run; pgf 66%
  non-zero at ~1e-5 m/s² (non-vacuous); Debug `-check all` clean; M1 advhor gate + 13/13 ctest still
  green. **Also fixed a `configure.sh` footgun** (`--debug` was clobbering the Release `build_intel_dp`;
  now goes to `build_intel_dp_debug`). See LESSONS L14.
- **M2.3 partial `compute_vel_rhs` byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs FESOM2 on the
  Coriolis + AB2 + PGF + SSH-gradient assembly: the new geometry field `coriolis`
  (`2·omega·sin(lat_geo)` via `r2g`, both hemispheres ±1.4e-4), the Coriolis term `uv_rhsAB_cor`
  (`UV·coriolis·elem_area`, ±8.8e6), AND the partial `UV_rhs` for BOTH the first-step Euler path
  (`uv_rhs_eul`, ff=1.0) and the AB2-steady path (`uv_rhs_ab2`, ff=ab2=1.6) — they differ over 66% of
  entries, so the ff branch is genuinely exercised. Built `src/oce/oce_dyn_velrhs.F90`
  (`compute_vel_rhs`, the non-advection part) + `compute_coriolis` in `src/mesh/mod_mesh_areas.F90`.
  Same gate `tools/run_pressure_gate.sh` (now **19 fields**): the FESOM2 shim drives the REAL
  `compute_vel_rhs` TWICE (lfirst Euler → AB2) with `momadv_opt=0` (skip `momentum_adv_scalar`),
  `ldiag_ke=.false.` (pi default is `.true.`!), a minimal fake `ice` (use_pice=0 on linfs → never
  dereferenced), `dt`/`r_restart` pinned. Momentum advection is **M2.4**. PASSED first run; Debug
  `-check all` caught a real `elnodes` shape-mismatch (`elem2D_nodes` is `MAX_NV=4` → slice `(1:3)`),
  fixed; M1 advhor gate + 13/13 ctest still green. See LESSONS L15.
- **M2.4 momentum advection / FULL `UV_rhs` byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs FESOM2
  on the new `uvnode_rhs` (the momentum-advection nodal intermediate, `w·du/dz` + `u·du/dx` post-
  normalize: 68.9% non-zero, both signs → both edge-scatter branches exercised) AND the re-gated
  `uv_rhsAB_cor`/`uv_rhs_eul`/`uv_rhs_ab2` (now Coriolis + momadv = the FULL `UV_rhs`), plus the
  prescribed input `w_e`. Ported `momentum_adv_scalar` (`oce_ale_vel_rhs.F90:335-589`) into
  `src/oce/oce_dyn_velrhs.F90` (3 passes: vertical `w·du/dz` on scalar CVs averaging elemental `UV`
  to prism faces ×`elem_area`×`w_e`, horizontal `u·du/dx` over edges via `edge_cross_dxdy`, then
  `×areasvol_inv` + 1-rank `exchange_nod` no-op + vertice→element `/3` ADD into `UV_rhsAB(1,1:2,·)`)
  and wired it into `compute_vel_rhs` at the `momadv_opt==2` site (FESOM2 :271-273, AFTER the
  Coriolis/PGF elem loop, BEFORE the AB blend). Same gate `tools/run_pressure_gate.sh` (now **21
  fields**; the shim flips `momadv_opt` 0→2 + prescribes `dynamics%w_e`, the FESOM3 driver allocates
  `dyn%w_e`/`dyn%work%uvnode_rhs` + prescribes the identical analytic `w_e`). PASSED first run (L9
  transitive-gate: every operand — `UV`/`elem_area`/`hnode`/`areasvol_inv`/`edge_cross_dxdy` +
  `nod_in_elem2D` & `edges` accumulation order — already byte-pinned). Debug `-check all` clean (the
  L15 trap avoided: `elem2D_nodes(1:3,el)` accessed per-component, no MAX_NV shape mismatch); M1
  advhor gate + 13/13 ctest still green. See LESSONS L16.
- **M2.4 biharmonic viscosity (`opt_visc=7`) byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs FESOM2
  on the first-stage Laplacian intermediate `visc_u_c`/`visc_v_c` (the element field `visc_filt_bidiff`
  pass 1 builds) AND the post-viscosity `uv_rhs_visc` (the gate target). Built `src/oce/oce_dyn_visc.F90`
  (`viscosity_filter` dispatcher + `visc_filt_bidiff` — a biharmonic = edge-based ∇² applied TWICE over
  INTERIOR edges only, free slip on the boundary; the FESOM2 `oce_dyn.F90:591-744` non-subcycl branch),
  a SEPARATE operator run AFTER `compute_vel_rhs` (FESOM2 `oce_ale.F90:3822`, NOT inside it like momadv).
  Added `visc_gamma0_h`/`visc_gamma1_h` to `t_dyn` (default 0 → pure biharmonic). **NO new geometry**
  (the `gradient_vec` worry was unfounded for opt_visc=7 — it reads only the geom/area-gated `edge_tri`/
  `elem_area`/`ulevels`/`nlevels`/`edge2D_in`). First gated kernel to use `edge2D_in` (interior-edge
  filter, byte-pinned transitively). Same gate `tools/run_pressure_gate.sh` (now **24 fields**); the shim
  forces `opt_visc=7` + pi gammas (`visc_gamma0=0.003` OVERRIDES the type default 0.03; `gamma_h=0`) and
  calls the REAL `visc_filt_bidiff`. **Bumped the shared prescribed `UV` to 2.0/1.5 m/s** so `|du|` spans
  all three flow-aware branches `max(γ0,γ1,γ2)` (selected on **20.4/78.6/1.0%** of edge-levels — a
  driver-side diagnostic verifies; γ2 is a near-dead production path, exercised synthetically). M2.3/M2.4
  re-gated `max|Δ|=0` with the new UV. PASSED first run (L9 transitive). Debug `-check all` clean (no L15
  trap — `edge_tri(:,ed)` is a size-2 slot, `elem2D_nodes` untouched); M1 advhor gate + 13/13 ctest still
  green. See LESSONS L17.
- **M2.5 implicit vertical viscosity (TDMA) byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs FESOM2 on the
  post-solve `uv_rhs_ivv` (the gate target) AND the prescribed inputs (`Av`/`stress_surf`/`w_i`). Built
  `src/oce/oce_dyn_ivertvisc.F90` (`impl_vert_visc_ale` — a per-element tridiagonal Thomas solve: implicit
  vertical viscosity `Av` + vertical advection [`w_i` upwind] + wind-stress top BC + quadratic bottom drag,
  OVERWRITING `UV_rhs` with the solution; `UV` is read-only). A SEPARATE operator run AFTER
  `viscosity_filter` (FESOM2 `oce_ale.F90:3874`, the `use_ssh_se_subcycl=.false.` branch). **The TDMA is a
  strictly SEQUENTIAL recurrence (forward+backward) → no summation/scatter order ambiguity → byte-matches by
  pure L9 transitivity** (every operand already pinned: `UV`/post-visc `UV_rhs`/`helem`/`zbar_e_bot`/levels +
  prescribed `w_i`/`Av`/`stress_surf`). Added `zbar_e_bot` to `t_mesh` (+ serialization; `helem` also now
  built in the driver, full cells). `Av`/`stress_surf` are passed as **explicit kernel args** (PP mixing M2.8 /
  forcing M2.10 not yet ported → prescribed analytically for this gate; when they land the caller sources
  them, kernel unchanged). Same gate `tools/run_pressure_gate.sh` (now **28 fields**); the shim pins
  `C_d=0.0025`, prescribes `Av`/`stress_surf`/`w_i`, and calls the REAL `impl_vert_visc_ale` (via an explicit
  interface block — no auto-gen `*_interface` module exists for it). PASSED first run; non-vacuous
  (`max|d(uv_rhs)|=0.985`, `w_i>0/<0` both branches fire). Debug `-check all` clean — pi has NO single-layer
  columns (min elem `nlevels=5`) so the FESOM2 single-layer `Z_n(0)`/`UV(:,0)` benign-OOB never triggers
  (deferred to M2.11 CORE2 shelf columns). M1 advhor gate + 13/13 ctest still green. See LESSONS L18.
- **Current task:** M2.6 SSH — stiffness matrix + `ssh_rhs` + CG solve (`oce_ale_ssh_splitexpl_*` / the
  implicit-SSH solver path). First M2 kernel that runs PAST forcing → assemble the reduced-M2 namelist
  (PP/no-GM/no-Redi/linfs/opt_visc=7) from `work_pi/namelist.*`. M1's multi-rank advection gate still rides M2.12.

## Geometry byte-gate (CLOSED ✓) — the 1-rank FESOM2 oracle recipe

The first true byte-gate vs the live oracle. Reusable for ALL M1+ kernel gates.
- **FESOM2 must run 1-rank** (so per-node `area`/`nod_in_elem2D` accumulation order ==
  FESOM3 global order; multi-rank reorders → ULP diffs). METIS can't make `dist_1`, so
  it is HAND-CRAFTED in `port2/fesom2/tests/data/MESHES/pi/dist_1/` (mirrors
  `save_dist_mesh` for np=1: `rpart.out` = npes/counts/identity node→contiguous map;
  `my_list00000.out` = identity lists; `com_info00000.out` = empty com_structs with
  BLANK lines for the zero-size halo arrays, which the reader's zero-trip `read(*,*)`
  skips). FESOM2's reader accepts it (validated).
- **FESOM2 geom dump:** `src/fesom_geom_dump.F90` (NEW), wired at end of `mesh_setup`
  (oce_mesh.F90), env-gated `FESOM_GEOM_DUMP`, npes==1 only, writes full arrays and
  STOPS before forcing (1-rank forcing init hangs on login node — irrelevant, geometry
  is complete at mesh_setup). Rebuilt: `build/bin/fesom.x` + `build/lib64/libfesom.so`
  (the proven oracle `bin/fesom.x`+`lib64/` is UNTOUCHED; backups `.proven` exist).
- **FESOM3 side:** `src/infra/mod_geom_dump.F90` + driver `src/drivers/fesom_geomdump.F90`
  (same binary format) + `tools/geom_diff.py` (field-by-field max|Δ|). Rotation matches
  the pi namelist (alpha/beta/gamma=50/15/-90, cyclic 360, force_rotation).
- **Run:** `tools/run_geom_gate.sh` → PASS. Individually: `tools/run_geomdump_pi.sh`
  (FESOM2) then `fesom_geomdump` (FESOM3) then `geom_diff.py`.

## M1.1–M1.4 advection byte-gate (CLOSED ✓) — the kernel-gate recipe (reusable for M1.5+)

Same 1-rank oracle as geometry, but the FESOM2 dump fires LATER (end of `ocean_setup`,
after `init_thickness_ale` builds `helem` + `muscl_adv_init` builds nboundary_lay/
edge_up_dn_tri), still BEFORE forcing (the 1-rank hang). The dump PRESCRIBES analytic
inputs and drives the REAL FESOM2 kernels, so it gates actual FESOM2 code:
- **FESOM2 side:** `src/fesom_advhor_dump.F90` (NEW), wired at end of `ocean_setup`
  (`oce_setup_step.F90`), env-gated `FESOM_ADVHOR_DUMP`, npes==1, STOPS after dumping.
  It sets `ttf`(nodes)/`vel`(elements) from an analytic formula of the rotated coords,
  then calls FESOM2's own `tracer_gradient_elements`→`fill_up_dn_grad`→`adv_tra_hor_upw1`/
  `_muscl`→`oce_tra_adv_flux2dtracer`. Built into `build/` (proven oracle `bin/*.proven`
  UNTOUCHED). dt=1800, num_ord=0.75 are pinned constants shared with FESOM3.
- **FESOM3 side:** `src/drivers/fesom_advhordump.F90` + `src/infra/mod_advhor_dump.F90`
  (FADVHDMP format, adds a 3-D-array writer). Same analytic prescription (byte-identical
  coords ⇒ byte-identical `ttf`/`vel`), same transcribed kernels.
- **Inputs that had to match:** `helem` (linfs/zstar agree at init since hbar=eta=0:
  `helem(nz,e)=zbar(nz)-zbar(nz+1)`), `areasvol` (geom-proven), and — crucially —
  `nod_in_elem2D` ORDERING, which the area gate already proved transitively (see L9).
- **Run:** `tools/run_advhor_gate.sh` → PASS (now gates BOTH horizontal M1.1 and vertical
  M1.2 fields). Individually: `tools/run_advhordump_pi.sh` (FESOM2) then `fesom_advhordump`
  (FESOM3) then `advhor_diff.py`.
- **M1.2 extension (DONE):** both the shim and the FESOM3 driver prescribe an analytic
  `wvel(nz,n)=1e-4·sin(2·lon)·cos(lat)·cos(0.3·nz)` (varies sign in space+depth → both
  upwind branches + QR4C interior exercised), call `adv_tra_ver_upw1`/`_qr4c` → scatter
  with a ZERO horizontal flux (so only the vertical branch fills `del_ttf_advvert`), and
  dump `wvel/zbar_3d_n/Z_3d_n/area/adv_flux_ver_{upw1,qr4c}/del_ttf_advvert_{upw1,qr4c}`.
  The FESOM2 shim uses the live `mesh%Z_3d_n/zbar_3d_n` (built by `init_ale` earlier in
  ocean_setup); the FESOM3 driver builds them from the init_ale formula (LESSONS L10).
- **M1.3 extension (DONE):** both the shim and the FESOM3 driver prescribe a second tracer
  `ttfAB` (sharp/large oscillatory, to force clipping), build `hnode/hnode_new`, and run the
  FCT branch inline — LO upw1(ttf) horiz+vert → `fct_LO`; HO MFCT(ttfAB)@opth + QR4C(ttfAB)@optv
  → antidiffusive `adv_flux_*`; `oce_tra_adv_fct` clips; `oce_tra_adv_flux2dtracer(use_lo=.TRUE.)`
  scatters. The shim drives the REAL FESOM2 `oce_tra_adv_fct` (passing a SEPARATE `aux_scratch`,
  not `edge_up_dn_grad`, so the latter survives for its dump). 15 FCT fields gated `max|Δ|=0`.
  Key facts in LESSONS L11: pi config = MFCT(opth=0)/QR4C(optv=1)/FCT; `edge_up_dn_grad`=grad(ttf)
  not grad(ttfAB); the a2 bignumber bottom-layer fill; the AUX cavity caveat.
- **M1.4 extension (DONE):** unlike M1.1–M1.3 (which inlined the orchestration in the shim), the
  M1.4 section drives FESOM2's REAL `init_tracers_AB` + `do_oce_adv_tra`. It runs AFTER all the
  M1.1/1.2/1.3 records are written (those driver calls overwrite tracers%work/tracers%data(1)),
  prescribing `values`=ttf (smooth) + `valuesold(1)`=ttfAB (sharp) on tracer 1 so the AB(2)
  interpolation `valuesAB = -(0.5+ε)·valuesold + (1.5+ε)·values` (ε=0.1) is exercised + gated
  (record `valuesAB`). Two configs on tracer 1: FCT (MFCT/QR4C/FCT, opth=0/optv=1) and non-FCT
  (MUSCL/QR4C/NON, ph=pv=0.75 → `do_zero_flux` + order knobs). pi runs `use_wsplit=.true.` in
  production but the shim FORCES `dynamics%use_wsplit=.false.` (= the FESOM3 dyn default) so the
  gate tests the matched explicit path (w==w_e; the implicit adv_tra_vert_impl is M2). del_ttf is
  accumulated (= advhoriz+advvert) as adv_tracers_ale does. 7 records: valuesAB,
  del_ttf_{advhoriz,advvert,}_step, del_ttf_{advhoriz,advvert}_stepnon, del_ttf_step_non. See L12.

## M1 entry notes (read before starting)

- M1 = first byte-exact physics slice (tracer advection, FCT, prescribed velocity),
  gate `max|Δ|=0` vs FESOM2 on pi. **This is the first task whose gate needs the
  instrumented FESOM2 oracle running** — also closes the M0.7 geometry byte-gate.
- Oracle prep required: build instrumented FESOM2 v2.7.3 (`port2/fesom2`), add the
  element dump shims (HANDOFF TODO above), run pi with prescribed UV/Wvel +
  controlled-input replay, produce reference dumps; compare with `tools/dump_diff.py`.
- Transcribe math from FESOM2 `oce_adv_tra_{hor,ver,fct}.F90` + `oce_muscl_adv.F90`;
  structure from dwarf `oce_adv_tra_*`. Geometry needed (elem_area, gradient_sca,
  edge_cross_dxdy, areas) is in `mod_mesh_areas` — verify it against the oracle here.
- **TODO (oracle-side, when running M1+ byte-gate):** add ELEMENT dump shims to
  FESOM2 `port2/fesom2/src/fesom_dump_shim.F90` (node shim exists; element fields
  uv/UV_rhsAB/Av/pgf_x/pgf_y need new shims, gid set [1000,2000,3000,4000,5000]).
- **pi mesh:** `/home/a/a270088/port2/fesom2/tests/data/MESHES/pi` (nod2D=3140,
  elem2D=5839, edge2D=8986; dist_2, dist_8). Nodes unique-owned; elem/edge myDim
  boundary-redundant (sum > global). com_info/my_list are ASCII free-field, 1-based.
- **Plan:** `docs/plans/2026-06-18-fesom3-architecture.md` (M0–M2 detailed, M3–M6+ roadmap).
- **Decisions:** project memory `project-brainstorm-decisions.md` (D0–D9).

## Build

```bash
./configure.sh --compiler intel --precision dp --clean --build   # anchor
./configure.sh --compiler gnu   --precision dp --clean --build   # portability
cd build_intel_dp && ctest --output-on-failure                   # self-tests
```
Anchor = Intel + DP + FESOM2-v2.7.3-exact flags (see docs/LESSONS.md L1). Build dirs:
`build_<compiler>_<precision>/` (Release); `--debug` builds into `build_<compiler>_<precision>_debug/`
(kept SEPARATE so a Debug binary never clobbers the Release anchor and silently breaks the byte-gates —
the L14 footgun, now fixed in configure.sh). Login-node runs of 1–8 ranks are fine for self-tests.
Debug is for `-check all` OOB/FPE only, NOT byte-comparable to the Release oracle (L10).

## Oracle — PROVEN RUNNABLE (2026-06-19) ✅

The whole byte-gate pipeline is validated end-to-end:
- `bin/fesom.x` (FESOM2 v2.7.3, git SHA 9271ae92) is prebuilt. The node dump shim is
  already wired into `oce_ale.F90`.
- **`tools/run_oracle_pi.sh [run_dir] [nsteps] [nranks]`** runs pi in ~0.25 s and writes
  `<prefix>.<rank>` dumps (density/pressure/bvfreq/sw_alpha-beta/Kv/ssh_rhs/d_eta/hbar/
  eta_n/w/T/S/hnode, substeps 1–16). Quirks it handles: fresh `fesom.clock` = `0 1 1948`
  twice; `run_length_unit='s'` = steps; absolute ClimateDataPath.
- The dump byte-format is **identical** to FESOM3 `mod_dump`; `tools/dump_diff.py` parses
  real FESOM2 dumps (150 records) and self-compares `max|Δ|=0`. A fresh run byte-matches
  the committed fixture `test/refdata/pi_oracle_default/dump.*` → FESOM2 is deterministic.
- M1 (DONE) used the end-of-`ocean_setup` shim `fesom_advhor_dump.F90` (prescribe inputs, call
  the real kernels, dump, stop before forcing). **M2.1 reuses that shim pattern** for EOS/pressure
  (see "Next task"); later M2 kernels that run past forcing need the reduced-M2 namelist (below).

## Canonical references

- **Algorithm oracle + byte-gate target:** FESOM2 v2.7.3 `/home/a/a270088/port2/fesom2/src/`
  (tag `fesom2.7.3-cport-instr`). Transcribe FROM here, gate `max|Δ|=0` AGAINST here.
  Run via `tools/run_oracle_pi.sh` (proven).
- **Structural template (structure only, NOT math/flags):**
  `/home/a/a270088/fesom3/design_refs/tracer_dwarf/lib/`.
- **Reduced M2 oracle namelist:** `mix_scheme='PP'`, `Fer_GM=.false.`, `Redi=.false.`,
  `which_ale='linfs'`, `opt_visc=7` (NOT the shipped KPP/GM CORE2 namelist).

## Environment (Levante)

- Toolchains via modules (see `env/levante.dkrz.de/shell.{intel,gnu}`): intel-oneapi
  2022.0.1 + openmpi 4.1.2-intel; gcc 11.2.0 + openmpi 4.1.2-gcc. netCDF loaded (only
  needed M2.10+). Login `gfortran` is 8.5.0 with no MPI — always build via `configure.sh`.

## M2.1 pressure/EOS/N² + M2.2 PGF + M2.3 vel_rhs + M2.4 momadv + M2.4 viscosity + M2.5 ivertvisc byte-gate (CLOSED ✓) — the dynamics-kernel-gate recipe (reusable for M2.6+)

Same end-of-`ocean_setup` 1-rank shim pattern as the M1 advection gate, the FIRST dynamics
kernel. `pressure_bv` runs DURING the timestep but 1-rank forcing hangs on the login node (L8),
so the shim fires at the end of `ocean_setup` (after `init_ale` built Z_3d_n/zbar_3d_n,
`init_thickness_ale` built hnode, `arrays_init` allocated+set density_ref=density_0), BEFORE
forcing. It PRESCRIBES analytic T/S and drives FESOM2's REAL `pressure_bv`, so it gates actual
FESOM2 code:
- **FESOM2 side:** `src/fesom_pressure_dump.F90` (NEW), wired at end of `ocean_setup`
  (`oce_setup_step.F90`, AFTER `advhor_dump_write`), env-gated `FESOM_PRESSURE_DUMP`, npes==1,
  STOPS. Sets `T/S = f(rotated coords, nz)`, FORCES the gate knobs (state_equation=1,
  which_ale='linfs', N2smth_v=.false./N2smth_hidx=1, ldiag_dMOC=.false., mix_scheme_nmb=-1 to
  skip dbsfc), then calls `pressure_bv` TWICE — once `N2smth_h=.false.` (raw bvfreq), once
  `.true.` (smoothed). Library rebuilt: `make -C port2/fesom2/build fesom.x` rebuilds
  `build/lib64/libfesom.so` (fesom.x loads it at runtime; the exe itself is not relinked, fine).
- **FESOM3 side:** `src/oce/oce_pressure_bv.F90` (`pressure_bv` + `densityJM_components` split EOS
  + `smooth_nod`), driver `src/drivers/fesom_pressuredump.F90` (builds Z_3d_n/zbar_3d_n/hnode like
  fesom_advhordump, same analytic T/S, same two-call raw/smoothed toggle), `tools/pressure_diff.py`
  + `tools/run_pressuredump_pi.sh` + `tools/run_pressure_gate.sh` — reuse the `mod_advhor_dump`
  FADVHDMP binary format. Added density_m_rho0/density_ref/hpressure/bvfreq to `t_dyn_work`.
- **Run:** `tools/run_pressure_gate.sh` → PASS (10 fields `max|Δ|=0`). Key facts in LESSONS L13:
  density_ref==density_0 (use_density_ref=.false.); the smoother byte-matches by faithful
  transcription (elem_area + nod_in_elem2D order are geom/area-proven, L9); the EOS split form +
  the two-call raw/smoothed dump; caller pre-zeros the outputs (pressure_bv leaves below-bottom
  entries as-is); cavity/use_density_ref branches transcribed but ungated on pi (M2.11).
- **M2.2 extension (DONE):** the SAME gate now also covers the hydrostatic PGF. After the two
  pressure_bv calls (hpressure unchanged by smoothing), both the shim and the FESOM3 driver pre-zero
  `pgf_x`/`pgf_y`, call the full-cell PGF (`pressure_force_4_linfs_fullcell` — shim calls the REAL
  FESOM2 routine, FESOM3 runs `oce_pgf`), and dump `pgf_x`/`pgf_y`. `pressure_diff.py` picks the two
  new records up automatically (12 fields now). PGF is element-based (nl-1, elem2D), the
  `gradient_sca`·`hpressure`/density_0 contraction — same shape as M1.1's tracer_gradient_elements,
  both operands already gated, so `max|Δ|=0` first run. See LESSONS L14 (+ the configure.sh
  `--debug` Release-clobber footgun fixed there).
- **M2.3 extension (DONE):** the SAME gate now also covers the partial vel_rhs. The geometry step grew
  `compute_coriolis` (`coriolis`/`coriolis_node` = `2·omega·sin(lat_geo)` via `r2g` on the rotated
  centroid/node; byte-matches because `r2g` reuses the g2r-proven rotation matrix — L9 transitive).
  Both the shim and the FESOM3 driver prescribe analytic `uv`(elements)/`eta_n`(nodes)/`uv_rhsAB`(prev),
  copy in the live `pgf_x`/`pgf_y`, and run `compute_vel_rhs` TWICE (lfirst Euler ff=1.0 → AB2 ff=ab2);
  the shim forces `momadv_opt=0`/`ldiag_ke=.false.`/`use_ssh_se_subcycl=.false.`, pins `dt`/`r_restart`,
  and passes a minimal fake `ice` (`use_pice=0` on linfs → `m_ice`/`m_snow` associated but never read).
  7 new records (coriolis, eta_n, uv_in, uv_rhsAB_prev, uv_rhsAB_cor, uv_rhs_eul, uv_rhs_ab2);
  `pressure_diff.py` picks them up automatically (19 fields). `max|Δ|=0` first run. See LESSONS L15.
- **M2.4 extension (DONE):** the SAME gate now also covers momentum advection → the FULL `UV_rhs`.
  `momentum_adv_scalar` lives in `src/oce/oce_dyn_velrhs.F90` (a private routine called by
  `compute_vel_rhs` when `momadv_opt==2`, mirroring FESOM2's own file layout). Both the shim and the
  FESOM3 driver prescribe an analytic `w_e`(nodes) `=1e-4·sin(2·lon)·cos(lat)·cos(0.3·nz)` (sign varies
  in space AND depth → non-trivial `w·du/dz`; the shim flips `momadv_opt` 0→2). 2 new records: `w_e`
  (input) + `uvnode_rhs` (the post-normalize momadv nodal intermediate); `uv_rhsAB_cor`/`uv_rhs_eul`/
  `uv_rhs_ab2` now carry momadv. `pressure_diff.py` picks them up automatically (**21 fields**). The
  momadv operator is non-vacuous (uvnode_rhs 68.9% non-zero, both signs); since the gate runs the REAL
  FESOM2 vertical+horizontal passes, FESOM3's transcription byte-matches BOTH by construction. Same L9
  transitive reason as all prior — `max|Δ|=0` first run. See LESSONS L16.
- **M2.4-visc extension (DONE):** the SAME gate now also covers biharmonic viscosity (`opt_visc=7`), a
  SEPARATE operator run AFTER the two `compute_vel_rhs` calls (FESOM2 `oce_ale.F90:3822`). FESOM3
  `src/oce/oce_dyn_visc.F90` (`viscosity_filter`→`visc_filt_bidiff`); the shim calls the REAL FESOM2
  `visc_filt_bidiff` (the leaf; the `viscosity_filter` dispatcher is pure branching, like M2.2). Both
  force `opt_visc=7` + pi gammas (`visc_gamma0=0.003`, `gamma_h=0`) and the shared `UV` is bumped to
  2.0/1.5 m/s (so `|du|` spans all three `max(γ0,γ1,γ2)` branches — selected 20.4/78.6/1.0%; a driver
  diagnostic prints it). 3 new records: `visc_u_c`/`visc_v_c` (the pass-1 Laplacian intermediate) +
  `uv_rhs_visc` (the gate target). `pressure_diff.py` picks them up automatically (**24 fields**). NO new
  geometry; `edge2D_in` (interior-edge filter) is byte-pinned transitively. `max|Δ|=0` first run; Debug
  `-check all` clean. See LESSONS L17.
- **M2.5-ivertvisc extension (DONE):** the SAME gate now also covers the implicit vertical viscosity TDMA, a
  SEPARATE operator run AFTER `visc_filt_bidiff` (FESOM2 `oce_ale.F90:3874`, the `use_ssh_se_subcycl=.false.`
  branch → `impl_vert_visc_ale`). FESOM3 `src/oce/oce_dyn_ivertvisc.F90` (`impl_vert_visc_ale`); the shim calls
  the REAL FESOM2 `impl_vert_visc_ale` via an **explicit interface block** (no auto-gen `*_interface` module
  exists for it — only `_vtransp` does). Both prescribe `Av` (strictly-positive vertical viscosity, all
  `nz=1..nl`), `stress_surf` (sign-varying wind stress) and `dynamics%w_i` (implicit vertical velocity, a
  DISTINCT sign-varying formula from `w_e`), and pin `C_d=0.0025`. `Av`/`stress_surf` are FESOM3 **explicit
  kernel args** (M2.8/M2.10 not ported); the driver builds `helem`+`zbar_e_bot` (full cells). 4 new records:
  `Av`/`stress_surf`/`w_i` (inputs) + `uv_rhs_ivv` (the post-solve gate target). `pressure_diff.py` picks them
  up automatically (**28 fields**). The TDMA is a strictly sequential Thomas recurrence (no order ambiguity),
  so `max|Δ|=0` first run by pure L9 transitivity; non-vacuous (`max|d(uv_rhs)|=0.985`, `w_i>0/<0` both fire).
  Debug `-check all` clean (pi min `nlevels=5` → the single-layer `Z_n(0)` benign-OOB never triggers). See
  LESSONS L18.

## Next task

M2.6 SSH — stiffness matrix + `ssh_rhs` + CG solve (plan Task M2.6). FESOM2 `oce_ale.F90` SSH path:
build `ssh_stiff` (CSR; negative factor `-g·dt·α·hbar`), assemble `ssh_rhs`, then the CG solve
(`soltol=1e-5`, `maxiter`). **Approach:**
- **This is the FIRST M2 kernel that runs PAST forcing** (the prior M2.1–M2.5 kernels depend only on
  prescribed T/S/UV/eta/w_e/w_i/Av/stress_surf + mesh, so the end-of-`ocean_setup` shim — which STOPS
  before forcing — sufficed). The CG solve needs the assembled stiffness matrix + `ssh_rhs`; decide
  whether the same prescribe-inputs-at-end-of-`ocean_setup` shim still works (prescribe `eta`/`d_eta`/
  `UV`/`hbar` and drive the REAL `ssh` assembly + CG) OR whether it needs the **reduced-M2 namelist**
  (PP/no-GM/no-Redi/linfs/opt_visc=7) assembled from `work_pi/namelist.*` to run the timestep further.
  Prefer the prescribe-and-stop shim if the SSH assembly + CG can be driven standalone.
- **Gate:** operator-diff `max|Δ|=0` on `ssh_rhs` / the CG solution / `d_eta`; plus the
  **`Σ ssh_rhs` over owned nodes telescopes to ~1e-13** consistency check (plan M2.6). The CG iteration
  is a sequence of dot-products + SpMVs — like the TDMA a deterministic recurrence given byte-identical
  operands, BUT the reduction order in the dot-products is the new bit-identity risk (1-rank: a single
  serial sum, so order-stable; the multi-rank reduction order is an M2.12 concern).
- **Watch:** the CSR stiffness matrix construction (sparsity pattern + values) is new structure; the
  `pcg` preconditioner choice; `g`/`alpha`/`hbar` operands. The L15 `elem2D_nodes(1:3,·)` MAX_NV-slice
  trap on any new element-indexed code.

- **Reduced-M2 namelist** (PP/no-GM/no-Redi/linfs/opt_visc=7) was NOT needed for M2.1–M2.5 (EOS + PGF +
  vel_rhs + momadv + viscosity + ivertvisc depend only on prescribed T/S/UV/eta/w_e/w_i/Av/stress_surf +
  mesh). Assemble it from `work_pi/namelist.*` if/when a kernel that runs PAST forcing (PP mixing, SSH CG)
  needs the full timestep (M2.6+).

**M1 multi-rank gate (folded into M2.12):** M2.12 builds the local-mesh remap (global→local
numbering/connectivity/`nod_in_elem2D` order/geometry, com-structs) needed for ANY multi-rank
byte-match, then gates the WHOLE model (incl. advection) on 1/8/32-rank. Re-confirm there that
the M1.1–M1.4 kernels' dropped halo exchanges (tr_xy/edge_up_dn_grad/fct_LO/fct_plus_minus/
del_ttf) + loop bounds (myDim vs myDim+eDim) are correctly lifted; the gate target is the
post-exchange OWNED values on the SAME partition (L8 accumulation-order caveat). Also deferred
to a richer mesh (M2.11 CORE2): the FCT `AUX`/`edge_up_dn_grad`-scratch cavity caveat (L11);
`enforce_cw_orientation`'s swap path (L8, pi has 0 swaps).

## Open notes / risks

- The FESOM2 oracle is built + PROVEN (M0.7 geometry + M1.1–M1.4 advection + M2.1–M2.5
  pressure/EOS/N²/PGF/vel_rhs/momadv/viscosity/ivertvisc gates all `max|Δ|=0`, 28 fields). Self-tests
  run standalone (13/13 ctest). M2+ kernel gates extend the proven 1-rank end-of-`ocean_setup` shim
  pattern (`fesom_pressure_dump.F90` now drives `pressure_bv` + `pressure_force_4_linfs_fullcell` + the
  REAL `compute_vel_rhs` with `momadv_opt=2` → `momentum_adv_scalar` + the REAL `visc_filt_bidiff` with
  `opt_visc=7` + the REAL `impl_vert_visc_ale` on the post-`viscosity_filter` `UV_rhs`); the next gate
  (M2.6 SSH stiffness/`ssh_rhs`/CG) is the FIRST that may need to run PAST forcing — see "Next task".
- **M2.12 is now heavy** (the local-mesh remap + the folded M1 advection multi-rank gate + the
  whole-model multi-rank byte-match + the deferred cavity/CW-swap caveats). Consider splitting the
  local-mesh remap into its own early-M2 task once a dynamics kernel first needs halos at multi-rank.
- The FESOM2 oracle shim edits (`fesom_advhor_dump.F90`, `fesom_pressure_dump.F90`, `oce_setup_step.F90`,
  geom/ale shims) live as UNCOMMITTED working-tree instrumentation in `port2/fesom2` (not committed
  there, by design); `libfesom.so` must be rebuilt (`make -C port2/fesom2/build fesom.x`, then re-run
  `cmake .` first if a NEW shim file was added so the GLOB picks it up) after editing a shim.
