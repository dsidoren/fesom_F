# FESOM3 — Handoff ARCHIVE (pre-M2.12 detail)

Archived from `docs/HANDOFF.md` on 2026-06-21 to keep the live handoff lean. This file
holds the full per-milestone completion writeups (M1 → M2.11c) and the detailed per-gate
recipes (geometry, M1 advection, M2.1–M2.9a dynamics kernels, M2.10 forcing, M2.11b IC).
Everything here is CLOSED, `max|Δ|=0` history. The live `docs/HANDOFF.md` keeps the
current status, a compact milestone summary, the next task, and the active entry notes.

---

## Completed-milestone detail (M1 → M2.11c) — verbatim from the old "Where we are"

- **M1 detail (all CLOSED ✓ on pi 1-rank).** **Geometry byte-gate** — FESOM3 mesh geometry is
  `max|Δ|=0` vs FESOM2 on pi (1-rank): elem_area, elem_cos, metric_factor, gradient_sca,
  edge_dxdy, edge_cross_dxdy, area/areasvol(+inv), coord_nod2D, elem2D_nodes, edges,
  edge_tri, all level arrays. Run it: `tools/run_geom_gate.sh`.
  ✅ **Caveat CLOSED (M2.11a, 2026-06-20):** the CW-swap path is now EMPIRICALLY gated — CORE2 needs
  **244654/244659** `enforce_cw_orientation` swaps (≈100%, mesh stored CCW; pi had 0), and the post-swap
  `elem2D_nodes` + centroid + `elem_area`/`gradient_sca`/`area` are all `max|Δ|=0` vs FESOM2's runtime
  `test_tri` on CORE2 1-rank (`tools/run_geom_gate_core2.sh`, 19 fields). See LESSONS L8 + **L26**.
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
- **M2.6 SSH (stiffness + `ssh_rhs` + CG solve) byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs FESOM2 on
  `ssh_stiff_diag` (the stiffness diagonal = mass + self-stiffness), `ssh_Aeta` (the FULL matvec `A·eta_n` →
  exercises every CSR non-zero), `ssh_rhs` (the edge-divergence rhs) AND `d_eta` (the CG solution — the gate
  target). Built `src/oce/oce_ssh_rhs.F90` (`init_stiff_mat_ale` CSR stiffness [`factor=g·dt·α·θ` ×
  `(zbar_e_bot−zbar_e_srf)`·gradient·edge_cross + mass `areasvol/dt`] + `compute_ssh_rhs_ale` edge-scatter
  divergence of `α(UV+UV_rhs)`) + `src/oce/oce_ssh_solve.F90` (`ssh_solve_preconditioner` MITgcm Jacobi-
  symmetrised M⁻¹ + `ssh_solve_cg` + `solve_ssh_ale`) — the **FIRST iterative solver** in the port.
  linfs builds the matrix ONCE (`update_stiff_mat_ale` skipped, `oce_ale.F90:3921`); the prescribe-and-stop
  shim STILL works (`init_stiff_mat_ale` runs at `ocean_setup:140`, before the shim → reduced-M2 namelist NOT
  needed). The CG converged in **37 iters** (byte-identical → the whole recurrence matches by L9). Same gate
  `tools/run_pressure_gate.sh` (now **32 fields**); the shim drives the REAL `compute_ssh_rhs_ale` +
  `solve_ssh_ale` via explicit interfaces (no auto-gen `*_interface`). 1-rank drops the global remap +
  `exchange_nod`/`MPI_Allreduce` (CG on local CSR `colind_loc`/`rowptr_loc`); `α=θ=1` ⇒ `(1−α)·ssh_rhs_old=0`;
  stiffness dt = pi namelist 86400/36=2400 (NOT the shim's 1800). `Σ ssh_rhs ≈ −3.7e-5` = ~1e-13 relative
  (telescoping). **CLEAN-rebuild footgun:** the first gate FAILED at uniform ~few-ULP on EVERY field because
  the incremental build (CMake reconfigure on the 2 NEW files → mixed `.mod` interfaces) ULP-drifted the
  codegen; `./configure.sh --clean --build` restored `max|Δ|=0`. Debug `-check all` clean; M1 advhor +
  13/13 ctest still green. See LESSONS L19.
- **M2.7 ALE (linfs) velocity/SSH/thickness-W update byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs FESOM2
  on `uv_upd` (post-`update_vel` velocity), `ssh_rhs_old` (the `compute_hbar_ale` edge-divergence intermediate),
  `hbar` + `dhe` (the new elevation + element thickness change), `eta_n_upd` (the α-blend), `w` (the
  `vert_vel_ale` vertical velocity — the big target), `hnode_new` (=`hnode`, the linfs branch leaves it
  untouched), `cfl_z` (the `compute_CFLz` intermediate), AND `w_split_e`/`w_split_i` (the `compute_Wvel_split`
  explicit/implicit split) — plus the prescribed `hbar_in`. Built `src/oce/oce_ale.F90` (`update_vel` [FESOM2
  `oce_dyn.F90:88-173`] + `compute_hbar_ale` + `update_eta_n` + `vert_vel_ale` + private `compute_CFLz`/
  `compute_Wvel_split` [FESOM2 `oce_ale.F90:2165/2323/3126/3217`]). The whole post-CG chain is `max|Δ|=0` by
  **L9 transitivity** (every operand pre-gated — `d_eta` M2.6, post-TDMA `UV_rhs` M2.5, `gradient_sca`/edge
  order/`helem`/`area`/`areasvol` geometry+M2.6; `hbar` prescribed identically; `dt=1800`/`θ=1`/`α=1` pinned).
  Same gate `tools/run_pressure_gate.sh` (now **43 fields**; the shim drives the REAL `update_vel`/
  `compute_hbar_ale`/`vert_vel_ale` via explicit interface blocks — no auto-gen `*_interface`, like
  `impl_vert_visc_ale`). linfs: the `zlevel`/`zstar` thickness redistribution + the `compute_hbar_ale`
  water-flux term vanish (`hnode_new` stays = `hnode`); Fer_GM/ldiag_ke branches dropped. **3 prescribed inputs
  (`eta_n`/`w_e`/`w_i`) are overwritten IN PLACE** by the M2.7 kernels → saved copies before the chain so the
  M2.3/M2.4/M2.5 input records still echo the prescription (L20). Non-vacuous: the Wvel split fired on 13253
  (nz,node) (`CFL_z>wsplit_maxcfl`, `use_wsplit=.true.`), `max|w|=0.042` m/s, `max|uv_upd|=2.5`. Debug
  `-check all` clean (compute clean — L15 `elem2D_nodes(1:3,·)` slices correct; the I/O dump writer needs
  `ulimit -s unlimited` for the big array temporary, a harness concern not a kernel bug). M1 advhor + 13/13
  ctest still green. See LESSONS L20.
- **M2.8 PP vertical mixing byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs FESOM2 on the prescribed
  nodal-velocity shear `uvnode` (input) AND the two gate targets `pp_Kv` (vertical diffusivity, nodes) /
  `pp_Av` (vertical viscosity, elements). Built `src/oce/oce_ale_mixing_pp.F90` (`oce_mixing_pp` — three
  SEQUENTIAL passes: pass 1 the inverse-Richardson factor `f=shear/(shear+5·max(N²,0)+1e-14)` stored in `Kv`,
  pass 2 `Av=mix_coeff_PP·mean₃(f²)+A_ver` on elements, pass 3 `Kv=mix_coeff_PP·f³+K_ver` on nodes — the
  ordering is load-bearing; + `Kv0_background_qiang`/`Kv0_background` transcribed but ungated, the
  `Kv0_const=.false.` path). Added `Kv`/`Av` to `t_dyn_work` (recomputed each step, NOT serialized) and
  `Kv0_const`/`use_instabmix`/`instabmix_kv`/`use_momix`/`momix_*`/`use_windmix`/`windmix_*` to
  `mod_param_phys`. A STANDALONE operator appended after M2.7 (the M2.9 step wires it in FRONT of
  `compute_vel_rhs`, and sources `impl_vert_visc_ale`'s `Av` from `dyn%work%Av`). **Prescribe `dyn%uvnode`
  DIRECTLY** (strong vertical shear → Ri factor spans **[5.8e-7, 0.93]**, 62.4% > 0.1, `max|Kv|`=8e-3 [800×
  K_ver], `max|Av|`=8.7e-3) rather than running `compute_vel_nodes` from the existing `UV` — isolates the add
  (no M2.3-M2.7 re-gate) and controls non-vacuity (the M2.5 prescribe-the-unsourced-input precedent;
  `compute_vel_nodes` is gated with the step at M2.9). Force the pi knobs (`A_ver=1e-4` — the pi NAMELIST
  override of the 1e-3 module default; `mix_coeff_PP=0.01`/`K_ver=1e-5`/`Kv0_const=.true.`). Same gate
  `tools/run_pressure_gate.sh` (now **46 fields**); the shim saves the M2.5 prescribed `Av` before PP
  overwrites the shared o_ARRAYS `Av`, dumps `uvnode`/`pp_Kv`/`pp_Av`. `max|Δ|=0` first run (L9 transitive —
  sequential passes, no scatter/order ambiguity). Debug `-check all` clean (the L15 `elem2D_nodes(1:3)` slice;
  `target` on the `dyn` dummy for the pointer aliases — the only compile fix). M1 advhor + 13/13 ctest still
  green. See LESSONS L21.
- **M2.8b `mo_convect` convective adjustment byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs FESOM2 on
  `moc_Kv` (post-adjustment nodal diffusivity) / `moc_Av` (post-adjustment elemental viscosity) AND the full
  M2.1-M2.8 cascade (the unstable T/S re-verifies every prior record). Built `src/oce/oce_mo_conv.F90`
  (`mo_convect` — the static-instability adjustment: `Kv(nz,node)=max(Kv,instabmix_kv)` where `bvfreq<0`;
  `Av(nz,elem)=max(Av,instabmix_kv)` where `any(bvfreq(nz,elnodes)<0)`; `use_windmix` transcribed guarded-off;
  `use_momix` TB04 OMITTED — needs forcing/ice not ported, deferred to M2.10). Run AFTER PP (the M2.9 step
  order). **Non-vacuity needed `bvfreq<0`** — added a localized unstable T band (`+8·max(0,cos lat·cos lon)·
  min(nz,8)/8`, a warm subsurface lens) → `bvfreq<0` on **5165 node-levels / 9936 elem-levels**, the floor
  fires (`max|ΔKv|`=0.090, `max|ΔAv|`=0.096, ~0.01→0.1). The T/S change CASCADES (T→density→…→every field) but
  all M2.1-M2.8 records re-verify `max|Δ|=0` automatically. Same gate (now **48 fields**); the shim saves the
  PP output (`pp_Kv_save`/`pp_Av_save`) before `mo_convect` overwrites `Kv`/`Av`, extends `ice_dummy` with
  `uice`/`vice`/`data(1)` (FESOM2 `mo_convect` pointer-assigns them unconditionally), and calls the REAL
  `mo_convect` via an explicit interface block. `max|Δ|=0` first run (L9 transitive — `max()` is order-free).
  Debug `-check all` clean; M1 advhor + 13/13 ctest still green. See LESSONS L22.
- **M2.9a tracer-solve assembly byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs FESOM2 on the per-tracer
  post-diffusion `del_ttf` (`tsol_del_ttf_T`/`tsol_del_ttf_S`) AND the solved T/S (`tsol_T`/`tsol_S`), plus the
  prescribed inputs (`Ki`, `heat_flux`, `virtual_salt`, `relax_salt`, the advection tendency `del_ttf_in`). Built
  `diff_tracers_ale` + `diff_part_hor_redi` (horizontal diffusion, Redi-off `Kh·tr_xy` edge scatter) +
  `diff_ver_part_impl_ale` (implicit vertical-diffusion Thomas/TDMA) + `bc_surface` (surface heat/virtual-salt/
  relaxation flux) into `src/oce/oce_ale_tracer.F90` (mirrors FESOM2's file). **The TDMA is the FIRST consumer of
  a PP output** — it reads `dyn%work%Kv` sourced LIVE post-`mo_convect` (no longer prescribed; the M2.8/M2.5
  prescribe-the-input precedent ends here). On pi it reduces to pure vertical diffusion + the surface BC (Redi off
  → `isredi=0`; FCT → `do_wimpl=.false.`; `i_vert_diff=.true.`; PP → no KPP nonlocal; `smooth_bh_tra=.false.`; linfs
  `hnode_new==hnode` → the ALE-reconstruct `T·(hnode−hnode_new)` term vanishes). Same gate
  `tools/run_pressure_gate.sh` (now **57 fields**): the shim prescribes `Ki`/surface-fluxes/`del_ttf_in`, sets the
  FCT/`i_vert_diff`/Redi-off config, computes `tr_xy` via the REAL `o_tracers::tracer_gradient_elements`, and drives
  the REAL FESOM2 `diff_tracers_ale` per tracer (T,S). Non-vacuous: Kv consumed=[0, 0.1] (live), `max|dT|`=1.26 °C,
  `max|dS|`=0.11, `max|del_ttf_T|`=0.50. PASSED first gate run. Debug `-check all` clean; M1 advhor + 13/13 ctest
  green. Two oracle traps fixed: `slope_tapered`/`tr_z` must be ZEROED in the shim (Redi-off reads them ×`isredi=0`
  before the multiply → uninitialised `NaN·0=NaN` would mismatch FESOM3 which omits the terms); and
  `tracer_gradient_elements` is a `o_tracers` MODULE procedure → `use` it (an explicit external interface gives
  `undefined symbol`). See LESSONS L23.
- **M2.9b `step_oce` ASSEMBLY byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs the REAL FESOM2
  `oce_timestep_ale` on **13 NODE substeps × 5 probes = 65 records** (density/pressure/bvfreq / Kv / ssh_rhs /
  d_eta / hbar / eta_n / hnode_new / w / T / S / hnode). Built `src/step/mod_step_oce.F90` (`step_oce` — the
  faithful sequence with LIVE data flow), added `compute_vel_nodes` + `update_thickness_ale` to
  `src/oce/oce_ale.F90`, `solve_tracers_ale` (the full per-tracer wrapper: `advect_tracer` → `tracer_gradient_elements`
  → `diff_tracers_ale` → salinity clamp) to `src/oce/oce_ale_tracer.F90`, and wired `model_ocean_step` into
  `src/step/mod_model.F90`. The gate drives the REAL `compute_vel_nodes` + `oce_timestep_ale` (its OWN built-in
  `dump_shim_record_node` — no new oracle dump code) vs FESOM3's `step_oce` (mirrored `mod_dump` dumps) on identical
  prescribed state. **This is the FIRST gate of the ASSEMBLY (data flow), not isolated kernels:** uvnode← `compute_vel_nodes`,
  Av← PP into `impl_vert_visc_ale` (the M2.5-deferred wiring), Kv← PP+convection into the tracer TDMA, d_eta← CG, T/S←
  advection+diffusion. `max|Δ|=0` first gate run (L9 transitive — every operand pre-pinned; M2.9b only adds wiring).
  Built `src/drivers/fesom_stepdump.F90`, FESOM2 shim `port2/fesom2/src/fesom_step_dump.F90` (forces the reduced-M2
  DISPATCH over the pi namelist, so all GM/Redi/KPP arrays stay allocated; SW_AB dumped but `--ignore-substep=2`'d),
  `tools/run_step{dump_pi,_gate}.sh` + `dump_diff.py --ignore-substep`. Non-vacuous (bvfreq<0 → convective Kv=0.1,
  d_eta∈[-4.7,2.2], T/S evolved); Debug `-check all` clean; M1 advhor + the M2.1-M2.9a pressure gate (57 fields) +
  13/13 ctest still green. **Scope:** `use_wsplit=.false.` (the FCT implicit vertical-advection `adv_tra_vert_impl`,
  do_oce_adv_tra's `use_wsplit=.true.` path, is a distinct unported kernel — M1.4 precedent; `impl_vert_visc_ale` still
  runs on the prescribed `w_i`). **The prescribe-and-stop shim STILL sufficed** (the "reduced-M2 namelist past forcing"
  worry was unfounded — `oce_timestep_ale` reads forcing ARRAYS, not files). Physical UV 0.50/0.40 m/s keeps eta_n in
  the ±10 `check_blowup` guard so the step completes cleanly. See LESSONS L24.
- **M2.10a forcing READ byte-gate CLOSED ✓ (on pi 1-rank) — the FIRST netCDF I/O.** `max|Δ|=0` vs the REAL FESOM2
  `sbc_do` on all 8 atmospheric fields (`u_wind`/`v_wind` [g2r-rotated] / `Tair` / `shum` / `shortwave` / `longwave` /
  `prec_rain` / `prec_snow`), read from the CORE2 NCAR stubs (`test/input/global/{u_10,v_10,q_10,ncar_rad,t_10,
  ncar_precip}.1948.nc`). Built `src/io/mod_io_netcdf.F90` (thin `use netcdf` wrapper) + `src/forcing/mod_forcing_read.F90`
  (julday [noleap=365·yyyy] + binarysearch + time-axis transform + periodic-lon halo + spatial bilinear + the two-stage
  time-interp `atmdata=rdate·coef_a+coef_b` ~710820-scale cancellation + g2r wind-coef rotation) + `vector_g2r` into
  `mod_mesh_rotate` + driver `src/drivers/fesom_forcingdump.F90`. **Added netCDF to FESOM3's CMake** (via `nf-config` +
  `-Wl,-rpath` → self-contained binaries; same spack `netcdf-fortran-4.5.3` the oracle links). **Fixed the L8 1-rank
  forcing hang** = a REAL infinite recursion in FESOM2 `next_io_rank` at npes==1 (async IO-rank selector) — patched
  `port2/fesom2/src/io_netcdf_workaround_module.F90` to short-circuit to sequential I/O on rank 0 (value-identical).
  Oracle shim `port2/fesom2/src/fesom_forcing_dump.F90` (after `forcing_setup`, drives REAL `sbc_do`, dumps, stops).
  Gate `tools/run_forcing_gate.sh` (8 fields). `max|Δ|=0` first gate run + after CLEAN rebuild; ctest 13/13 + the M2.1-M2.9a
  pressure gate (57 fields) still green (no regression). See LESSONS L25.
- **M2.10b bulk transfer coeffs + wind stress byte-gate CLOSED ✓ (on pi 1-rank, 17 fields).** `max|Δ|=0` vs FESOM2 on
  `cd_atm_oce`/`ch_atm_oce`/`ce_atm_oce` (the REAL `ncar_ocean_fluxes_mode` — Large&Yeager 2004 + Large-2009 drag, a
  5-iteration Monin-Obukhov stability solve) + `stress_atmoce_x`/`stress_atmoce_y` (`Cd·ρ_air·|Δu|·Δu`) + `stress_surf`
  (the 2×elem2D node→elem `/3` average, `a_ice=0`) + the prescribed `sst`/`srfoce_u`/`srfoce_v`. Built
  `src/forcing/mod_forcing_bulk.F90` (`forcing_bulk_ncar` + `forcing_wind_stress` + `forcing_stress_surf`). The bulk
  reads the LIVE M2.10a `u_wind`/`v_wind`/`Tair`/`shum`; SST + surface velocity are PRESCRIBED (a dummy ice supplies the
  thermo type-defaults `inv_rhoair=1./1.3`/`tmelt=273.15`/`rhoair=1.3`; the M2.5/M2.8 prescribe-the-unsourced-input
  pattern). Byte traps (LESSONS L25): `inc_ratio=1.0e-4`/`inv_rhoair=1./1.3`/`tmelt=273.15`/`rhoair=1.3` transcribed
  VERBATIM (un-suffixed default-real literals); `(ustar*ustar)` not `**2`; `atan(1.0_WP)` runtime; `elem2D_nodes(1:3)`
  (L15). The shim extends `fesom_forcing_dump.F90` (REAL bulk + inline stress); same gate `tools/run_forcing_gate.sh`
  (now **17 fields**). Non-vacuous (`Cd∈[5.5e-5,1.7e-2]`, both stability branches via the wide Tair−SST range); `max|Δ|=0`
  first gate run; Debug `-check all` clean; ctest 13/13 + the pressure gate (57 fields) still green.
- **M2.10c shortwave penetration byte-gate CLOSED ✓ (on pi 1-rank, 20 fields total).** `max|Δ|=0` vs FESOM2 on the
  floored `chl`, the visible-removed `heat_flux_sw`, and `sw_3d` (the 48×nod2D depth profile) — the REAL
  `cal_shortwave_rad` (Morel&Antoine/Sweeney 2005: `swsurf=(1-albw)·shortwave·0.54`, `chl` floor 0.02, the v1/v2/sc1/sc2
  polynomial in `log10(chl)`, the two-exponential `sw_3d` profile over `zbar_3d_n`, `/vcpw`). Built
  `src/oce/oce_shortwave_pene.F90`. Consumes the LIVE M2.10a `shortwave`; `chl`/`heat_flux`/`a_ice=0` prescribed;
  `albw=0.066_WP` forced both sides (avoids the un-suffixed default-real ambiguity); the live `mesh%zbar_3d_n`
  (init_ale, pressure-gate-proven) drives the depth profile. The shim extends `fesom_forcing_dump.F90` (a dummy
  `ice%data(1)`=a_ice + the REAL `cal_shortwave_rad` via an explicit interface). Non-vacuous (chl floor fires on 476
  polar nodes, `sw_3d` decays over depth). `max|Δ|=0` first gate run; Debug `-check all` clean; ctest 13/13 + pressure
  gate (57 fields) still green. **M2.10 forcing COMPLETE.**
- **M2.11b initial conditions (`do_ic3d`) byte-gate CLOSED ✓ (on CORE2 1-rank).** `max|Δ|=0` vs the REAL FESOM2
  `oce_initial_state` on `ic_temp` (potential T, post-`insitu2pot`), `ic_salt` AND `Z_3d_n` (the init_ale depth the
  interpolation consumes, dumped as an input check). Built `src/oce/oce_initial_state.F90` (`do_ic3d` + `nc_readGrid` +
  `nc_ic3d_ini` + `getcoeffld` + `extrap_nod3D`) reading **phc3.0_winter.nc** (360×180×33): the FIRST 3D netCDF read
  (`nc_get_var3d_dp` added to `mod_io_netcdf`, real(8)) + NaN→dummy mask (phc3.0 land = NaN, no `_FillValue`) + periodic-lon
  halo + spatial bilinear + vertical LINEAR interp onto `Z_3d_n` + `extrap_nod3D` (Gauss-Seidel neighbour-average sweep +
  downward fill — the partition-order step, deterministic at 1-rank) + Kelvin guard + **`insitu2pot`** (Bryden-1973 RK4
  `ptheta`/`atg`, added to `oce_pressure_bv.F90`; uses 1-D `Z(nz)` for the pressure proxy). `idlist=2,1` ⇒ salt (data(2))
  read FIRST, temp (data(1)) SECOND, `t_insitu=.true.`. Reused `forcing_binarysearch` (identical bisection) + `mod_io_netcdf`.
  Gate `tools/run_ic_gate_core2.sh`: the FESOM2 shim `port2/fesom2/src/fesom_ic_dump.F90` (env `FESOM_IC_DUMP`, npes==1) simply
  DUMPS the live `tracers%data(1)/(2)%values` the REAL `do_ic3d` already produced (no prescribe — `oce_initial_state` runs at
  `oce_setup_step.F90:253`, before the end-of-ocean_setup dump) + stops; the FESOM3 driver `src/drivers/fesom_icdump.F90` reads
  the CORE2 mesh, builds linfs `Z_3d_n`, runs the transcribed `do_ic3d`, dumps. Oracle namelist reduced-M2 override
  (`which_ALE→linfs` so its `Z_3d_n` = FESOM3's linfs build; CORE2 ocean_setup completes 1-rank in ~19 s, no hang — the dump
  stops before forcing). `max|Δ|=0` FIRST gate run; Debug `-check all` clean; pressure gate (57) + ctest (13) still green. The
  T/S `do_ic3d` "global min" print loops WET levels only (5.628) vs the FESOM3 driver's full-array `minval` (0.0, bottom-zeroed)
  — a print-mask red herring, NOT a data mismatch (the full-array byte-gate is 0). See LESSONS L27.
- **M2.11c (full multi-step lifecycle on CORE2) — ✅ DONE; the WHOLE dynamical core (incl. the CG `d_eta`) byte-matches `max|Δ|=0` over multiple steps (the "CG floor" was SOLVED — L29).** Built the FIRST real time-stepping run (not a prescribe-and-stop shim):
  FESOM3 `src/drivers/fesom_lifecycle.F90` (cold-start CORE2 mesh + `do_ic3d` phc3.0 IC + N-step runloop calling
  `mod_step_oce::step_oce`, multi-step AB2 evolving in place) vs the REAL FESOM2 multi-step lifecycle
  (`tools/run_lifecycle_core2.sh`, the built-in per-substep `dump_shim` over N steps, `use_ice=.false.` UNFORCED).
  **Result:** the **entire ported dynamical core is byte-identical (`max|Δ|=0`) on the 40× CORE2 mesh** — proven by the
  NEW `tools/run_pressure_gate_core2.sh` (the M2.1-M2.9 per-kernel pressure gate, run on CORE2): all 27 fields
  `max|Δ|=0` incl. the FULL SSH stiffness matrix (`ssh_stiff_diag` + `ssh_Aeta` matvec) + `ssh_rhs` + PP mixing +
  the tracer solve. **ONE field initially diverged — the preconditioned-CG `d_eta` (`~4e-14` on CORE2) — first read
  as an "iterative-solver reproducibility floor" (old L28), then SOLVED 2026-06-21 (L29): it was a BUG, an
  auto-vectorised preconditioner divide.** Root cause: `ssh_solve_preconditioner`'s off-diagonal divide compiled to
  packed `divpd` in FESOM3 but scalar `divsd` in the oracle (FESOM3 writes the component `ssh_stiff%pr_values`,
  provably non-aliasing → vectorised; the oracle writes a local pointer → scalar); packed vs scalar divide differ
  ~1 ULP under `-no-prec-div -fimf-use-svml`, so the **un-gated** `pr_values` drifted in 1299/870146 entries, seeding
  `z=M⁻¹r`, sub-ULP through iters 1–5, surfacing in the CG residual at ~iter 6 (CORE2 136 iters; pi's 37 stayed below
  the bit — why pi was `max|Δ|=0`). The earlier "ruled out: preconditioner" was WRONG (`pr_values` was never dumped).
  **Fix: one `!DIR$ NOVECTOR`** → scalar `divsd`, `pr_values` `max|Δ|=0`. Now `run_pressure_gate_core2.sh` PASSes all
  28 fields incl. `d_eta`, and the multi-step `run_lifecycle_gate_core2.sh` MATCHes (195 records, worst |Δ|=0). The
  whole dynamical core is byte-exact across multiple steps on CORE2. Localisation method + generalisable lessons in
  **LESSONS L29**. Debug `-check all` clean; FESOM3 lifecycle step-1 diagnostics match the oracle (eta_n=0.34864,
  uv=0.23185). M2.11c-2 (forced: prescribe the oracle's per-step `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt`/
  `stress_surf` via a flux dump shim, `use_ice=.true.`) is ✅ DONE — it inherited the same byte-exact CG; the forced
  gate `tools/run_lifecycle_forced_gate_core2.sh` MATCHes (195 records, worst |Δ|=0, re-verified post-L29 2026-06-21).
  **M2.11a geometry ✅ + M2.11b initial conditions ✅ + M2.11c dynamical-core-on-CORE2 ✅** (`run_geom_gate_core2.sh` 19,
  `run_ic_gate_core2.sh` 3, `run_pressure_gate_core2.sh` 27, all CORE2 1-rank `max|Δ|=0`). **SCOPE carried (L25):**
  `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt` air-sea budget = M3; M1 multi-rank gate rides M2.12.


---

## Gate recipes (geometry + M1 advection + M1 entry notes)

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


---

## M2.1–M2.9a dynamics-kernel byte-gate (CLOSED ✓, `run_pressure_gate.sh`, 57 fields) + M2.9b whole-step ASSEMBLY gate (CLOSED ✓, `run_step_gate.sh`, 65 records) — the dynamics-gate recipes (reusable for M2.10+)

**M2.9b step-ASSEMBLY gate (a DISTINCT vehicle from the per-kernel FADVHDMP gate below).** Instead of dumping
each kernel's output in the FADVHDMP format (`pressure_diff.py`), the step gate reuses FESOM2's OWN built-in
per-substep `dump_shim_record_node` instrumentation inside `oce_timestep_ale` (the `run_oracle_pi` / `mod_dump`
format), which FESOM3's `mod_dump` is byte-identical to. So the FESOM2 shim `fesom_step_dump.F90` only prescribes
the clean state + forcing arrays + reduced-M2 config and calls the REAL `compute_vel_nodes` + `oce_timestep_ale`
ONCE (no new oracle dump code); FESOM3's `step_oce` (`mod_step_oce.F90`) emits the same records inline; `dump_diff.py
--glob --ignore-substep=2` compares the 13 NODE substeps (the SW_AB substep is FESOM2-only, dead in M2). Run it:
`tools/run_step_gate.sh`. This gates the ASSEMBLY (LIVE data flow), the per-kernel gate below gates each leaf.
The two coexist (separate scratch files, separate comparators). See LESSONS L24.

### the per-kernel FADVHDMP recipe (M2.1–M2.9a, reused since M1):

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
- **M2.6-SSH extension (DONE):** the SAME gate now also covers the free-surface solve — `init_stiff_mat_ale`
  (CSR stiffness, FESOM2 `oce_ale.F90:1584`) + `compute_ssh_rhs_ale` (`:2012`) + `solve_ssh_ale` (`:3272`)
  + the CG/preconditioner (`solver.F90`), run AFTER `impl_vert_visc_ale` (FESOM2 `oce_ale.F90:3920-3930`).
  FESOM3 `src/oce/oce_ssh_rhs.F90` + `oce_ssh_solve.F90` (the FIRST iterative solver). The shim drives the
  REAL `compute_ssh_rhs_ale` + `solve_ssh_ale` via **explicit interface blocks** (no auto-gen `*_interface`,
  like `impl_vert_visc_ale`); the stiffness was ALREADY built by `init_stiff_mat_ale` at `ocean_setup:140`
  (BEFORE the shim, with the pi namelist dt=2400 — NOT the shim's 1800), and linfs never updates it. Both
  force `α=θ=1.0` (FESOM2 default ⇒ `(1−α)·ssh_rhs_old=0`), zero `d_eta` (CG x0) + `ssh_rhs_old`, chain the
  rhs off the post-TDMA `uv_rhs`. FESOM3 builds its stiffness with `dt_ssh=86400._WP/real(36,WP)` (FESOM2's
  formula, NOT a 2400.0 literal — the L7 trap). 4 new records: `ssh_stiff_diag` (matrix diagonal), `ssh_Aeta`
  (full matvec `A·eta_n` → every CSR non-zero), `ssh_rhs`, `d_eta` (CG solution). `pressure_diff.py` picks
  them up automatically (**32 fields**). CG converged in **37 iters** (byte-identical ⇒ the whole recurrence
  matches by L9); `Σ ssh_rhs≈−3.7e-5` (~1e-13 relative). 1-rank drops the global remap + `exchange_nod`/
  `MPI_Allreduce` (local CSR `colind_loc`/`rowptr_loc`); serial DO-loop dot-products (OpenMP-off form) +
  `sum()` matvecs. **The first gate FAILED at uniform ~few-ULP on EVERY field — the incremental build's
  CMake-reconfigure-on-new-files left mixed `.mod` codegen; a CLEAN `./configure.sh --clean --build` fixed
  it (always clean-rebuild after adding NEW files).** Debug `-check all` clean. See LESSONS L19.
- **M2.7-ALE extension (DONE):** the SAME gate now also covers the post-CG velocity/SSH/thickness-W update —
  `update_vel` (FESOM2 `oce_dyn.F90:88`) + `compute_hbar_ale` (`oce_ale.F90:2165`) + the inline `eta_n` blend
  (`:3973`) + `vert_vel_ale` (`:2323` → `compute_CFLz` `:3126` + `compute_Wvel_split` `:3217`), run AFTER
  `solve_ssh_ale` (FESOM2 `oce_ale.F90:3946-4081`). FESOM3 `src/oce/oce_ale.F90`. The shim drives the REAL
  routines via **explicit interface blocks** (no auto-gen `*_interface`, like `impl_vert_visc_ale`); the
  `eta_n` blend is inlined in both (FESOM2 keeps it inline in the step). Both prescribe an analytic `hbar`
  (the previous-step elevation), force `use_wsplit=.true.`/`wsplit_maxcfl=1.0` (pi production), and chain off
  the post-CG `d_eta`/post-TDMA `UV_rhs`. **3 prescribed inputs (`eta_n`/`w_e`/`w_i`) are overwritten IN PLACE**
  (eta_n by the blend, w_e/w_i by `compute_Wvel_split`) → saved copies before the chain so the M2.3/M2.4/M2.5
  input records still echo the prescription (dump `eta_n_upd`/`w_split_e`/`w_split_i` as NEW records). 11 new
  records: `hbar_in`, `uv_upd`, `ssh_rhs_old`, `hbar`, `dhe`, `eta_n_upd`, `w`, `hnode_new`, `cfl_z`,
  `w_split_e`, `w_split_i`; `pressure_diff.py` picks them up automatically (**43 fields**). `max|Δ|=0` first
  run (L9 transitive — every operand pre-gated). Non-vacuous: the Wvel split fired on 13253 (nz,node), `max|w|=
  0.042` m/s. Debug `-check all` clean (`ulimit -s unlimited` for the I/O writer's big array temp — L20). See
  LESSONS L20.
- **M2.8-PP extension (DONE):** the SAME gate now also covers PP Richardson-number mixing (FESOM2
  `oce_ale.F90:3728` → `oce_mixing_PP`, the first dynamics operator, here appended STANDALONE after M2.7).
  FESOM3 `src/oce/oce_ale_mixing_pp.F90` (`oce_mixing_pp` — 3 sequential passes: Ri factor in `Kv`, then
  `Av=mix·mean₃(f²)+A_ver` elements, then `Kv=mix·f³+K_ver` nodes; `Kv0_background_qiang`/`Kv0_background`
  transcribed, ungated). The shim calls the REAL `oce_mixing_PP` via an **explicit interface block** (free
  subroutine, no auto-gen `*_interface`, like `impl_vert_visc_ale`). **Both PRESCRIBE `dynamics%uvnode`
  directly** (strong vertical shear, identical formula — the Ri factor spans [5.8e-7, 0.93], non-vacuous)
  rather than running `compute_vel_nodes` from `UV` (isolates the add; `compute_vel_nodes` gated at M2.9).
  Force `A_ver=1e-4` (pi NAMELIST, not the 1e-3 module default), `mix_coeff_PP=0.01`/`K_ver=1e-5`/
  `Kv0_const=.true.`. The shim saves the M2.5 prescribed o_ARRAYS `Av` (`Av_in`) before PP overwrites it; both
  pre-zero `Kv`/`Av`. 3 new records: `uvnode` (input), `pp_Kv`/`pp_Av` (gate targets); `pressure_diff.py`
  picks them up automatically (**46 fields**). `max|Δ|=0` first run (L9 transitive — sequential passes, no
  scatter/order ambiguity). Debug `-check all` clean. See LESSONS L21.
- **M2.8b-mo_convect extension (DONE):** the SAME gate now also covers the convective adjustment (FESOM2
  `oce_ale.F90:3729` → `mo_convect`, run AFTER PP). FESOM3 `src/oce/oce_mo_conv.F90` (`mo_convect`:
  `Kv`/`Av`=max(·, `instabmix_kv`=0.1) where `bvfreq<0`; `use_windmix` transcribed guarded-off; `use_momix`
  TB04 OMITTED → forcing/ice deferred to M2.10). The shim calls the REAL `mo_convect` via an **explicit
  interface block**; its `(ice, partit, mesh)` signature pointer-assigns `ice%uice`/`vice`/`data(1)`
  unconditionally → extend `ice_dummy` to allocate them (even though `use_momix=.false.` never reads them).
  **Non-vacuity: a deliberate unstable T band** (`+8·max(0,cos lat·cos lon)·min(nz,8)/8`, a warm subsurface
  lens) makes `bvfreq<0` on **5165 node-levels / 9936 elem-levels** so the floor fires (`max|ΔKv|`=0.090). The
  FIRST gate to perturb the shared T/S — it CASCADES through every downstream field, but all M2.1-M2.8 records
  re-verify `max|Δ|=0` (both sides identical T). Save the PP output (`pp_Kv_save`/`pp_Av_save`) before
  `mo_convect` overwrites `Kv`/`Av` (the L20 pattern); 2 new records `moc_Kv`/`moc_Av` (**48 fields**). `max|Δ|=0`
  first run. Debug `-check all` clean. See LESSONS L22.
- **M2.9a-tracer-solve extension (DONE):** the SAME gate now also covers the per-tracer diffusion solve (FESOM2
  `oce_ale_tracer.F90:335` → `diff_tracers_ale`, run per tracer inside `solve_tracers_ale` after advection). FESOM3
  `src/oce/oce_ale_tracer.F90` grew `diff_tracers_ale` + `diff_part_hor_redi` + `diff_ver_part_impl_ale` +
  `bc_surface`. The shim drives the REAL FESOM2 `diff_tracers_ale` (via `diff_tracers_ale_interface`) per tracer
  (T,S), ISOLATED on prescribed inputs (the M2.8 pattern): `del_ttf` enters with a prescribed advection tendency
  (`del_ttf_in`), the horizontal diffusivity `Ki` + surface fluxes (`heat_flux`/`water_flux`/`virtual_salt`/
  `relax_salt`) are prescribed, `tr_xy` is set by the REAL `o_tracers::tracer_gradient_elements`, and the implicit
  vertical-diffusion TDMA CONSUMES the LIVE `o_ARRAYS Kv` (post-`mo_convect`, the **FIRST consumer of a PP
  output**). Force the gate config: `Redi=.false.` (isredi=0), `tra_adv_lim='FCT'` (`do_wimpl=.false.`),
  `i_vert_diff=.true.`, `use_sw_pene`/`use_kpp_nonlclflx`/`use_icebergs=.false.`, `smooth_bh_tra=.false.`,
  `is_nonlinfs=0`. **CRITICAL oracle trap:** zero `slope_tapered`/`tr_z` in the shim — Redi-off reads them
  ×`isredi=0` BEFORE the multiply, and the shipped `Fer_GM/Redi=.true.` namelist allocated-but-didn't-fill them →
  uninitialised `NaN·0=NaN` would poison the result (FESOM3 omits those terms). 9 new records: `tsol_Ki`/
  `tsol_heat_flux`/`tsol_virtual_salt`/`tsol_relax_salt`/`tsol_del_ttf_in` (inputs) + `tsol_del_ttf_T`/
  `tsol_del_ttf_S`/`tsol_T`/`tsol_S` (gate targets); `pressure_diff.py` picks them up automatically (**57 fields**).
  Non-vacuous (Kv consumed=[0,0.1] live, `max|dT|`=1.26 °C, `max|dS|`=0.11). `max|Δ|=0` first run; Debug `-check
  all` clean. The 3D `relax_to_clim` (`clim_relax=0` on pi) + the explicit `diff_ver_part_expl_ale` (`i_vert_diff`)
  + biharmonic `diff_part_bh` (`smooth_bh_tra=.false.`) are transcribed-deferred (not on the pi path). See LESSONS L23.

## M2.10 forcing byte-gate (CLOSED ✓, `run_forcing_gate.sh`, 20 fields: read + bulk + SW penetration) — the FIRST netCDF I/O, a NEW gate vehicle

A DISTINCT gate from the per-kernel FADVHDMP pressure gate (it reads REAL files). The FESOM2 oracle shim
`port2/fesom2/src/fesom_forcing_dump.F90` is wired AFTER `forcing_setup` (`fesom_module.F90:316`, NOT
end-of-`ocean_setup` — `sbc_ini` only runs inside `forcing_setup`), env-gated `FESOM_FORCING_DUMP`, npes==1: it pins
the model time (yearold:=yearnew to skip the year-change reload; `timenew=43200`), drives the REAL `sbc_do`, maps
`atmdata→u_wind/.../prec_snow` (`update_atm_forcing:681-694`), dumps (FADVHDMP), stops. **Prereq: the np=1
`next_io_rank` fix** in `port2/fesom2/src/io_netcdf_workaround_module.F90` (else the async IO-rank selector
infinite-recurses — the L8 hang). FESOM3 side: `src/io/mod_io_netcdf.F90` (`use netcdf` wrapper) + transcribed
`src/forcing/mod_forcing_read.F90` + `vector_g2r` (in `mod_mesh_rotate`) + driver `src/drivers/fesom_forcingdump.F90`,
hardcoding the pi forcing config (8 CORE2 files + var names + iyear=1948/freq=1). **CMake now links netCDF** (`nf-config`
+ rpath). Run it: `tools/run_forcing_gate.sh` (8 fields `max|Δ|=0`); the diff reuses `pressure_diff.py` (its trailing
"PRESSURE/EOS…" summary line is cosmetic — the per-field PASS table is the result). The dump drivers don't call
`set_partition` → loop over `mesh%nod2D` (`frc%nnod`), not `partit%myDim_nod2D`. See LESSONS L25.

## M2.11b initial-conditions byte-gate (CLOSED ✓, `run_ic_gate_core2.sh`, 3 fields) — the FIRST 3D netCDF read, on CORE2

A DISTINCT gate vehicle from the per-kernel pressure gate and the forcing gate: it runs on the **CORE2 mesh** (1-rank) and
reads the **phc3.0_winter.nc** climatology, but the oracle shim PRESCRIBES NOTHING — the REAL `oce_initial_state`→`do_ic3d`
runs DURING `ocean_setup` (`oce_setup_step.F90:253`), so by the end-of-ocean_setup dump point the live
`tracers%data(1)/(2)%values` already hold the final potential-T/S. So `port2/fesom2/src/fesom_ic_dump.F90` (env `FESOM_IC_DUMP`,
npes==1, wired FIRST before `advhor_dump_write`) just DUMPS `data(1)/(2)%values` + `mesh%Z_3d_n` (FADVHDMP) and stops. FESOM3
side: `src/oce/oce_initial_state.F90` (`do_ic3d` + helpers) + `insitu2pot`/`ptheta`/`atg` in `oce_pressure_bv.F90` +
`nc_get_var3d_dp` in `mod_io_netcdf.F90` + driver `src/drivers/fesom_icdump.F90` (reads CORE2 mesh, builds linfs `Z_3d_n`,
runs `do_ic3d`, dumps). The oracle namelist gets the reduced-M2 override (`which_ALE 'zlevel'→'linfs'` so its init `Z_3d_n`
= FESOM3's linfs full-cell build — gated as an input; `mix_scheme→PP`, `Fer_GM/Redi→.false.`). Run it:
`tools/run_ic_gate_core2.sh` (3 fields `max|Δ|=0`, ~19 s; reuses `pressure_diff.py`). **CORE2 ocean_setup completes at
1-rank** (no hang — the dump stops before forcing). Reuses `forcing_binarysearch` (identical bisection). See LESSONS L27.

