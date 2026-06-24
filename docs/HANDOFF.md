# FESOM3 — Handoff (durable state across sessions)

Single source of truth for "where are we / what's next". Update at the end of every task.
Full pre-M2.12 milestone detail + the per-gate recipes live in [`HANDOFF-archive.md`](HANDOFF-archive.md).

## Where we are

- **Milestone:** M2 (minimal ocean dynamical core) — **multi-rank MVP byte-match ACHIEVED (tag `m2-mvp`).**
  **→ M3 (sea ice EVP) ✅ COMPLETE 2026-06-23 (all M3a–M3f DONE, `max|Δ|=0` 1-rank AND multi-rank) — ✅ COMMITTED
  (`b591153`, tag `m3`).**
  **→ M4 (GM/Redi) ✅ COMPLETE 2026-06-24 (all M4a–M4f DONE, `max|Δ|=0` 1-rank AND multi-rank) — ✅ COMMITTED (`7a088de`, tag `m4`).**
  M4a (producers) + M4b (GM diffusivity/streamfunction/bolus-velocity) + M4c (GM bolus into advection) + M4d (Redi
  isopycnal diffusion) + M4e (GM+Redi in the FORCED/fully-native lifecycle) + M4f (multi-rank) all byte-gated
  `max|Δ|=0` vs FESOM2 (CORE2 1-rank AND dist_2/dist_8, both whichEVP). Plan: `docs/plans/2026-06-23-m4-gm-redi.md`.
  **→ M5 (KPP vertical mixing + sw_pene) IN PROGRESS 2026-06-24.** Scoped + decomposed M5a–M5d
  (plan: `docs/plans/2026-06-24-m5-kpp.md`). **M5a (isolated KPP module gate, 95 pressure-gate
  fields) ✅ DONE; M5b (wire KPP into `step_oce` + UNFORCED CORE2 lifecycle) ✅ DONE 2026-06-24:
  `max|Δ|=0` KPP-alone 195 (3-step) + 325 (5-step) AND KPP+GM+Redi 195** (`tools/run_lifecycle_kpp_gate_core2.sh`).
  M5a-1 (init + wscale) + M5a-2 (ri_iwmix) + M5a-3 (prestep dVsq/ustar/Bo + dbsfc + bldepth) +
  M5a-4 (blmix_kpp + enhance + combine + viscAE average) all byte-exact on BOTH pi AND CORE2, all
  FOUR first try; M5b assembled the real `oce_mixing_kpp_driver` + the `mix_scheme_nmb==1` dispatch
  (one new physics step beyond M5a: `smooth_blmc` — see the M5 block + LESSONS **L45**).
  **RESUME at M5c (forced/fully-native lifecycle + sw_pene: the `ghats` nonlocal + `sw_3d` tracer
  TDMA terms + wire `cal_shortwave_rad`).** See "Milestone ordering" + the "↳ RESUME HERE" section below.
  (M3 scoped 2026-06-22 into M3a–M3f.) M3a (ice foundation + cold-start IC + FCT mass
  matrix) ✅ DONE 2026-06-22 (`max|Δ|=0`, 4 fields, CORE2 1-rank, `tools/run_ice_gate_core2.sh`). M3b (ocean2ice + EVP
  dynamics) ✅ DONE 2026-06-22 (`max|Δ|=0`, 7 fields × BOTH whichEVP=0 standard-EVP AND whichEVP=1 mEVP, CORE2 1-rank,
  `tools/run_evp_gate_core2.sh`). M3c (ice FCT advection) ✅ DONE 2026-06-22 (`max|Δ|=0`, 6 fields = `rhs_a/m/ms` +
  post-advection `a_ice/m_ice/m_snow`, × BOTH whichEVP=0 AND whichEVP=1, CORE2 1-rank, `tools/run_icefct_gate_core2.sh`;
  first try). M3d (thermodynamics: `cut_off` + `thermodynamics` + `therm_ice`/`budget`/`obudget`/`flooding`) ✅ DONE
  2026-06-23 (`max|Δ|=0`, 6 fields = post-thermo `a_ice/m_ice/m_snow` + `flx_h`/`flx_fw` + `t_skin`, × BOTH whichEVP=0
  AND whichEVP=1, CORE2 1-rank, `tools/run_icethermo_gate_core2.sh`; first try, L38). M3e (`oce_fluxes` coupling-out —
  THE PAYOFF) ✅ DONE 2026-06-23 (`max|Δ|=0`, 5 fields = `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt` +
  `stress_surf`, × BOTH whichEVP=0 AND whichEVP=1, CORE2 1-rank, `tools/run_iceflux_gate_core2.sh`; first try, L39 —
  the prescribed M2.11c-2 flux dump is now NATIVE). **M3f-1 (`ice_timestep` assembly) + M3f-2 (native-flux coupled
  multi-step forced lifecycle) ✅ DONE 2026-06-23 (first try, L40): `max|Δ|=0` on the 195-record (3-step) AND
  325-record (5-step) NODE substep set, BOTH whichEVP=0 AND whichEVP=1 — the native sea-ice + air-sea budget
  (`ocean2ice → ice_timestep → oce_fluxes_mom → oce_fluxes`) drives the multi-step coupled CORE2 ocean step
  byte-exactly; `tools/run_lifecycle_native_gate_core2.sh`. M3f-3a (native CORE2 NCAR forcing READ) + M3f-3c (native
  NCAR bulk Ch/Ce + stress_atmoce + wind-on-ice stress_atmice) ✅ DONE 2026-06-23 (L41) — 14/16 atmospheric arrays
  native, `max|Δ|=0` self-checked in the coupled lifecycle, both EVP variants. **M3f-3b (runoff + SSS monthly
  climatology — the last 2/16) ✅ DONE 2026-06-23 (first try, L42) → M3f-3 COMPLETE: the FULLY-NATIVE CORE2 lifecycle
  (NO prescribed atmosphere — the whole air-sea forcing computed in-driver) is `max|Δ|=0` vs FESOM2, 195 records (3-step)
  AND 325 records (5-step), BOTH whichEVP=0 and whichEVP=1 (`tools/run_lifecycle_fullynative_gate_core2.sh`). **M3f-4
  (multi-rank) ✅ DONE 2026-06-23 (first try, L43) → M3 COMPLETE: the FULLY-NATIVE lifecycle byte-matches FESOM2
  `max|Δ|=0` at MULTI-RANK (CORE2 dist_2 + dist_8, BOTH whichEVP, 195 + 325 records) —
  `tools/run_lifecycle_fullynative_gate_multirank.sh`. The M3a–e kernels were already MR-ready (M2.12 optional-partit);
  the only new work was `mesh%ocean_area` allreduce + the `mod_forcing_bulk` nNodL loop bound + the
  `fesom_lifecycle_native_mr` driver. Tag `m3`.**
  **CURRENT STATUS (2026-06-22):** the whole M2 dynamical core is BYTE-COMPLETE — `max|Δ|=0` vs FESOM2 on pi 1-rank
  (per-kernel + whole-step) AND across multi-step CORE2 (unforced + forced lifecycle, incl. the free-surface CG) AND
  **the WHOLE ASSEMBLED ocean step now byte-matches FESOM2 at MULTI-RANK** (pi dist_2 + dist_8). **M2.12 (multi-rank,
  MVP exit) is ✅ DONE: M2.12a (geometry) + M2.12b (tracer advection) + M2.12c-1 (pre-SSH DYNAMICS → ssh_rhs) + c-2
  (SSH stiffness + free-surface CG → d_eta) + c-3 (post-SSH ALE update + tracer SOLVE + the WHOLE step via step_oce,
  all 65 substep records `max|Δ|=0`).** The whole-model multi-rank byte-match — the architectural MVP — holds at dist_2
  AND dist_8 with no 1-rank regression. **MVP CLOSED 2026-06-22: baseline re-confirmed after a clean rebuild (all 8
  byte-gates `max|Δ|=0`/green) → M2.12b–c committed → tag `m2-mvp` created.** **PRODUCTION VALIDATION DONE 2026-06-22:
  the multi-rank FREE-RUNNING CORE2 lifecycle (dt=1800, reduced-M2) is now BYTE-EXACT vs FESOM2 (195 records
  `max|Δ|=0`); it exposed two latent multi-rank bugs — a UV halo over-exchange (FIXED, eDim) and a `tr_xy` halo "drift"
  that turned out to be an OpenMPI `vader` KNEM single-copy corruption of large (>131 KB) messages, fixed by
  `env.sh` exporting `OMPI_MCA_btl_vader_single_copy_mechanism=none` (LESSONS L35; not a FESOM3 bug).** Next milestone:
  M3 (sea ice EVP + `oce_fluxes` air-sea budget). See "Next task".

A compact list; full per-milestone detail + the gate recipes are archived in
[`HANDOFF-archive.md`](HANDOFF-archive.md).

### Completed milestones (all `max|Δ|=0` vs FESOM2 at the stated anchor)
- **M0** Foundation — tag `m0`. CMake + params/types/partitioning/halo/dump/mesh/step; 13/13 ctest (Intel+GNU dp).
- **M1** Tracer advection — tag `m1`, pi 1-rank. Geometry + horiz/vert advection (upw1/MUSCL/QR4C) + FCT (Zalesak) limiter + assembled step (`tools/run_advhor_gate.sh`).
- **M2.1–M2.9** Minimal dynamical core — pi 1-rank. EOS/`hpressure`/N² + PGF + Coriolis/AB2 `vel_rhs` + momentum advection + biharmonic viscosity + implicit vertical viscosity (TDMA) + SSH (CSR stiffness + preconditioned CG) + ALE update + PP mixing + `mo_convect` + tracer solve, AND the whole `step_oce` assembly (`tools/run_pressure_gate.sh` 57 fields + `tools/run_step_gate.sh` 65 records).
- **M2.10** Forcing — pi 1-rank, the FIRST netCDF I/O. Read + NCAR bulk coeffs/wind stress + SW penetration (`tools/run_forcing_gate.sh` 20 fields).
- **M2.11a** CORE2 geometry — CORE2 1-rank, 19 fields (`tools/run_geom_gate_core2.sh`); closed the CW-swap caveat (≈244654/244659 swaps).
- **M2.11b** Initial conditions (`do_ic3d`) — CORE2 1-rank, the FIRST 3D netCDF read (`tools/run_ic_gate_core2.sh`).
- **M2.11c** Multi-step CORE2 lifecycle — the whole dynamical core byte-exact over multiple steps, UNFORCED + FORCED (`tools/run_lifecycle_gate_core2.sh` + `…_forced_…`, 195 records each). The CORE2 free-surface-CG "reproducibility floor" was SOLVED — a vectorised preconditioner divide, one `!DIR$ NOVECTOR` (LESSONS **L29**).
- **M2.12a** Multi-rank local-mesh remap + geometry gate — pi dist_2 + dist_8, 19 fields every rank (`tools/run_geom_gate_multirank.sh`); the FIRST multi-rank byte-gate (LESSONS **L30**). Detail in "Next task".
- **M2.12b** Multi-rank tracer advection gate — pi dist_2 + dist_8, 13 fields every rank (`tools/run_advhor_gate_multirank.sh`): the `find_neighbors` halo dance + `exchange_elem_full(tr_xy/elem_area)` + `exchange_nod(fct_LO/fct_plus/minus)` + owned/halo loop bounds, all byte-exact; `del_ttf` (FCT + MUSCL) `max|Δ|=0` (LESSONS **L31**). Detail in "Next task".
- **M2.12c-1** Multi-rank pre-SSH DYNAMICS chain — pi dist_2 + dist_8, `max|Δ|=0` (25 records: density/pressure/bvfreq/Kv/ssh_rhs) (`tools/run_stepdyn_gate_multirank.sh`). Lifted the whole chain `compute_vel_nodes → pressure_bv(+smooth_nod) → pressure_force_4_linfs → oce_mixing_pp → mo_convect → compute_vel_rhs(+momentum_adv_scalar) → viscosity_filter(visc_filt_bidiff) → impl_vert_visc_ale → compute_ssh_rhs_ale` to multi-rank via the M2.12b OPTIONAL-`partit` pattern: `mod_part_bounds.local_dims` (owned/halo/full bounds) + the FESOM2 exchanges (`exchange_nod` UVnode/bvfreq-per-smooth-sweep/UVnode_rhs/ssh_rhs, `exchange_elem` U_c/V_c) + the visc interior-edge test on the GLOBAL id (`myList_edge2D(ed)>edge2D_in`). New `exchange_nod` rank-3 node-block variant in `mod_halo`. **NO 1-rank regression** (step 65 + advhor 40 + lifecycle 195 + 13/13 ctest all `max|Δ|=0`). See LESSONS **L32**. Detail in "Next task".
- **M2.12c-2** Multi-rank SSH stiffness + free-surface CG — pi dist_2 + dist_8, `max|Δ|=0` on `d_eta` (30 records = the 25 c-1 + 5 `d_eta`; `tools/run_stepdyn_gate_multirank.sh`, substep 9 now compared). The FIRST multi-rank ITERATIVE solver (38 CG iters on 8 ranks, byte-exact). Lifted `init_stiff_mat_ale` + `solve_ssh_ale` (`ssh_solve_preconditioner`/`ssh_solve_cg`) via the same OPTIONAL-`partit` pattern: owned loops (`nNodO`), arrays `nNodL`, `exchange_nod(diag_values)` in the precond + `exchange_nod(pp/rr/x)` in the CG + new `allreduce_sum` (in `mod_halo`) for the dot-products; `rtol`/convergence denominators stay GLOBAL (`mesh%nod2D`). **KEY FINDING: NO mesh-infra extension was needed** — the HANDOFF/L32 hypothesis (extend `elem2D_nodes`/`gradient_sca` to the halo + `enforce_cw` on owned+eDim) was WRONG; the stiffness owned rows assemble FULLY LOCALLY (verified partition invariant: for owned edges both triangles are owned — FESOM2's own `elem2D_nodes`/`gradient_sca` are owned-only too). **NO 1-rank regression** (step 65 + pressure 57 + advhor MR + 13/13 ctest all `max|Δ|=0`). See LESSONS **L33**. Detail in "Next task".
- **M2.12c-3** Multi-rank post-SSH ALE update + tracer SOLVE + the WHOLE assembled step — pi dist_2 + dist_8, `max|Δ|=0` on ALL **65 substep records** (density/pressure/bvfreq/Kv/ssh_rhs/d_eta/hbar/eta_n/hnode_new/w/T/S/hnode; `tools/run_step_gate_multirank.sh`). **The whole multi-rank ocean step byte-matches FESOM2 — the architectural MVP.** Lifted `update_vel`/`compute_hbar_ale`/`update_eta_n`/`vert_vel_ale`(+`compute_CFLz`/`compute_Wvel_split`)/`update_thickness_ale` (`oce_ale.F90`) + `solve_tracers_ale`/`diff_tracers_ale`/`diff_part_hor_redi`/`diff_ver_part_impl_ale` (`oce_ale_tracer.F90`; the FCT advection was M2.12b) via the OPTIONAL-`partit` pattern (`owned_bounds`/`local_dims` + the FESOM2 exchanges), and threaded the optional `partit` through `mod_step_oce::step_oce` (absent ⇒ 1-rank verbatim; the 1-rank callers `fesom_stepdump`/`fesom_lifecycle` unchanged). New whole-step MR driver `fesom_stepfull_mr` (= `step_oce` through `partit`, gid-keyed per-rank dumps); the FESOM2 oracle is UNCHANGED (`fesom_step_dump` npes>1 already drives the whole `oce_timestep_ale` + dumps all substeps). pi ships only dist_1/2/8 (no dist_32; dist_2+dist_8 = the c-1/c-2 coverage). **NO 1-rank regression** (step 65 + pressure 57 + advhor MR dist_2/8 + stepdyn MR dist_2/8 + 13/13 ctest all `max|Δ|=0`). See LESSONS **L34**. Detail in "Next task".

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


## Gate recipes (archived)

The detailed per-gate recipes — the 1-rank FESOM2 oracle setup (hand-crafted `dist_1`, the
end-of-`ocean_setup` prescribe-and-stop shim, the FADVHDMP dump format), the M1 advection
kernel-gate recipe, the M2.1–M2.9a dynamics-kernel recipe (+ every per-kernel extension), the
M2.10 forcing gate, and the M2.11b IC gate — are in [`HANDOFF-archive.md`](HANDOFF-archive.md).
The live "Oracle" section above + the "M2.12 entry notes" below cover what the current work needs.

## Next task

**M2.10 forcing DONE + M2.11a geometry DONE + M2.11b initial conditions DONE + M2.11c lifecycle DONE (unforced AND
forced) — the multi-step CORE2 lifecycle is `max|Δ|=0` (the "CG floor" was SOLVED 2026-06-21; the forced re-verify
CONFIRMED 2026-06-21: 195 records, worst |Δ|=0).** **All of M2 is now byte-exact end-to-end on CORE2.**

**→ M2.12 (multi-rank, MVP exit) IN PROGRESS.** Decomposed into M2.12a/b/c in the plan
(`docs/plans/2026-06-18-fesom3-architecture.md` Task M2.12) — see "M2.12 entry notes" below. **Key finding: the
multi-rank foundation ALREADY EXISTS** (`mod_halo.F90` real MPI exchange_nod/elem + `mod_partitioning` myList/com,
tested 1/2/8-rank), so M2.12 is faithful transcription of FESOM2's `read_mesh`/`find_neighbors`/`mesh_areas` (with
their halo exchanges) gated PER-RANK vs same-partition FESOM2 — not a from-scratch parallel build.
- **M2.12a — local-mesh remap + per-rank GEOMETRY byte-gate — ✅ DONE (2026-06-21).** `max|Δ|=0` on ALL 19 geometry
  fields, every rank, on pi **dist_2 + dist_8** (`tools/run_geom_gate_multirank.sh`). `read_mesh_local` (global→local
  scatter + owned `nod_in_elem2D` + `enforce_cw`(owned)) + partition-aware `compute_geometry` (owned centers
  `exchange_elem`'d for owned-edge `edge_cross_dxdy`; owned-node areas local). See LESSONS L30.
- **M2.12b — multi-rank TRACER ADVECTION byte-gate — ✅ DONE (2026-06-21).** `max|Δ|=0` on ALL 13 advection fields,
  every rank, pi **dist_2 + dist_8** vs same-partition FESOM2 (`tools/run_advhor_gate_multirank.sh`). Lifted the
  advection subtree (init_tracers_AB → tracer_gradient_elements → muscl_adv_init/find_up_downwind_triangles/
  fill_up_dn_grad → adv_tra_hor/ver → oce_tra_adv_fct → do_oce_adv_tra → advect_tracer) to multi-rank via an
  **OPTIONAL `partit`** (absent ⇒ the proven 1-rank path VERBATIM; present+npes>1 ⇒ owned/halo bounds + the FESOM2
  exchanges). The explicit-shape dummies + per-element arithmetic are UNCHANGED so codegen — and the byte-match incl.
  any vectorised divide (L29) — is preserved; only loop bounds change + guarded exchanges are added. New machinery:
  `src/infra/mod_part_bounds.F90` (`owned_bounds`/`is_multirank`); the `find_neighbors` halo dance
  (`complete_nod_in_elem_halo` in `read_mesh_local` — fills halo-node `nod_in_elem2D`, re-localised through eXDim);
  `exchange_elem_full` (2D-r/2D-i/3D-r on `com_elem2D_full`) + `core_blk_r` in `mod_halo`; `compute_geometry` halo
  exchanges (`elem_area` full + `area`/`areasvol`/`*_inv`); `exchange_elem(tr_xy)` (full) + `exchange_nod(fct_LO)` +
  `exchange_nod(fct_plus/minus)`; `find_up_downwind_triangles` `coord_elem`/`e_nodes` (full-halo, global-id match);
  `nboundary_lay` owned+halo, NO exchange (partition-local — matches FESOM2). Gate: `tools/run_advhordump_multirank.sh`
  (oracle `fesom_advhor_dump.F90` npes>1 branch) + `fesom_advhordump_mr` (F3) → `advhor_diff.py` per rank.
  ⚠️ `edge_up_dn_grad` is captured PRE-`do_oce_adv_tra` (FESOM2 REUSES it as the FCT `AUX` scratch → bignumber=1e3).
  **No regression:** 1-rank geom(19)/advhor(40)/pressure(57)/step(65) + CORE2 lifecycle(195) + 13/13 ctest all
  `max|Δ|=0`/green. See LESSONS L31.
- **M2.12c-1 — multi-rank pre-SSH DYNAMICS chain — ✅ DONE (2026-06-22).** `max|Δ|=0` on pi dist_2 + dist_8 (25
  records density/pressure/bvfreq/Kv/ssh_rhs, `tools/run_stepdyn_gate_multirank.sh`). Lifted EVERY dynamics kernel
  up to `compute_ssh_rhs_ale` with the M2.12b optional-`partit` pattern (`mod_part_bounds.local_dims` for the
  owned/halo/full bounds + the exact FESOM2 exchanges, arithmetic UNCHANGED). New machinery: `mod_halo` rank-3
  `exchange_nod` node-block variant (UVnode/UVnode_rhs `(2,nl-1,nod)`); `local_dims` (adds `nEdgeL`/`nElemL` to
  `owned_bounds`). The exchanges added at the FESOM2 sites: `exchange_nod(UVnode)` (compute_vel_nodes),
  `smooth_nod` per-sweep `exchange_nod(bvfreq)`, `exchange_nod(UVnode_rhs)` (momentum_adv_scalar),
  `exchange_elem(U_c/V_c)` (visc_filt_bidiff, between its two Laplacian sweeps), `exchange_nod(ssh_rhs)`
  (compute_ssh_rhs_ale). visc interior-edge test lifted to the GLOBAL id (`myList_edge2D(ed)>edge2D_in`). Gate:
  FESOM3 `fesom_stepdump_mr` (dynamics chain only, no step_oce — the un-lifted CG is never reached) + the FESOM2
  `fesom_step_dump` shim extended to npes>1 (owned-prescribe + `exchange_elem(UV)`, runs the REAL `oce_timestep_ale`
  whose built-in per-rank gid-keyed `dump_shim_record_node` emit the probes; the post-ssh_rhs substeps 9/11/12/13/15/16
  are `--ignore-substep`'d since FESOM3 stops at ssh_rhs). Key finding: c-1 needs NO mesh-infra change — every kernel
  reads `elem2D_nodes`/`gradient_sca` only at OWNED elements (the rest read UV/UV_rhs/helem/elem_area at the M2.12b-
  exchanged full halo). See LESSONS **L32**.
- **M2.12c-2 — multi-rank SSH stiffness + free-surface CG — ✅ DONE (2026-06-22).** `max|Δ|=0` on `d_eta`, pi dist_2
  + dist_8 (30 records, `tools/run_stepdyn_gate_multirank.sh` with substep 9 un-ignored; 38 CG iters on 8 ranks). The
  FIRST multi-rank ITERATIVE solver. **KEY FINDING — the scoped "mesh-infra extension" (this bullet's prior text +
  L32) was WRONG and NOT needed.** A `dist_N`-file check (for each owned edge, is `edge_tri ≤ myDim_elem2D`? → 0 halo
  on every rank) proved the stiffness owned rows assemble FULLY LOCALLY. Invariants (pi dist_2/8; universal in FESOM2
  since its own `elem2D_nodes`/`gradient_sca` are owned-only, `oce_mesh.F90:497`/`:2466`): (i) every edge incident to
  an owned node is owned; (ii) both triangles of an owned edge are owned; (iii) owned-element nodes ≤ `nNodL`. So
  `init_stiff_mat_ale` lifts with ONLY bounds (`nod2D→nNodO`, `edge2D→nEdgeO`; `n_num(nNodL)`, `n_pos(12,nNodO)`;
  `ssh_stiff%dim=mesh%nod2D` global) — NO `elem2D_nodes`/`gradient_sca` halo, NO `enforce_cw` on halo, NO exchange in
  the assembly. The CG (`oce_ssh_solve`): arrays `nNodL`, owned loops `nNodO`, `exchange_nod(diag_values)` in the
  precond + `exchange_nod(pp/rr/x)` per-iter + NEW `mod_halo::allreduce_sum` (scalar+vec, `MPI_IN_PLACE`/`MPI_SUM`/
  `MPI_COMM_FESOM` = the oracle `solver.F90`) for the owned-partial-sum dot-products; `rtol`/convergence denominators
  GLOBAL (`mesh%nod2D`). Optional `partit` absent ⇒ 1-rank VERBATIM. The cross-rank `MPI_Allreduce` byte-matches
  (same OpenMPI+comm-size+op+8-byte-type → deterministic tree, L6); the L29 NOVECTOR precond divide carries over. Gate:
  `fesom_stepdump_mr` extended (`init_stiff_mat_ale`+`solve_ssh_ale`, dump `d_eta` at `DUMP_SUBSTEP_SSH_SOLVE`=9), the
  oracle already dumps `d_eta` at substep 9 (NO oracle change). **NO 1-rank regression** (step 65 + pressure 57 + advhor
  MR dist_2 + 13/13 ctest all `max|Δ|=0`). See LESSONS **L33**.
- **M2.12c-3 — multi-rank post-SSH ALE update + tracer SOLVE + the WHOLE step — ✅ DONE (2026-06-22).** `max|Δ|=0` on
  ALL **65 substep records** (density/pressure/bvfreq/Kv/ssh_rhs/d_eta/hbar/eta_n/hnode_new/w/T/S/hnode), pi dist_2 +
  dist_8 (`tools/run_step_gate_multirank.sh`, only SW_AB id=2 ignored = the 1-rank `run_step_gate.sh` set). **The whole
  multi-rank ocean step byte-matches FESOM2 — the architectural MVP.** Lifted EVERY post-SSH kernel via the M2.12b/c
  optional-`partit` pattern (bounds from `owned_bounds`/`local_dims`, arithmetic UNCHANGED): `update_vel` (owned
  elements + `exchange_elem_full(UV)` — rank-3, the eDim+eXDim superset of FESOM2's `exchange_elem`); `compute_hbar_ale`
  (ssh_rhs_old zero owned+halo, edge-div over owned edges, hbar update owned, then `exchange_nod(hbar)` ALWAYS — outside
  the linfs guard — so the owned-element `dhe` reads valid halo hbar); `update_eta_n` (owned+halo, no exchange);
  `vert_vel_ale` (Wvel zero owned+halo, edge-div owned, cumsum/divide owned, `exchange_nod(w)` + `exchange_nod(hnode_new)`)
  + `compute_CFLz`/`compute_Wvel_split` (owned+halo); `update_thickness_ale` (linfs no-op, helem unchanged + already
  full-halo-valid ⇒ FESOM2's `exchange_elem(helem)` is value-neutral and skipped, L33 "no machinery you don't need").
  Tracer SOLVE (`oce_ale_tracer.F90`): `solve_tracers_ale` (per-tracer `exchange_nod(values)` after the diffusion solve +
  the salinity clamp over owned+halo; `tr_xy` sized to the LOCAL `nElemF`); `diff_tracers_ale` ALE-reconstruct owned;
  `diff_part_hor_redi` owned edges (owned-node del_ttf complete, invariant i; `Ki`/`areasvol` read at the edge's halo
  node are prescribed/M2.12b-exchanged); `diff_ver_part_impl_ale` owned-node TDMA. Threaded the optional `partit` through
  `mod_step_oce::step_oce` (absent ⇒ 1-rank VERBATIM; the 1-rank callers `fesom_stepdump`/`fesom_lifecycle` unchanged).
  New whole-step MR driver `src/drivers/fesom_stepfull_mr.F90` (= `step_oce` through `partit`, local-sized state + tracer
  machinery + forcing inputs, gid-keyed per-rank dumps); the FESOM2 oracle is UNCHANGED (`fesom_step_dump` npes>1 already
  drives the whole `oce_timestep_ale` + its built-in dump_shim emits all substeps). The byte-gate on the rich analytic
  state (non-trivial T/S/UV/SSH at 2 AND 8 ranks) subsumes the suggested rest-at-rest / gravity-wave / `stale_halo_max_*`
  sanity probes — matching FESOM2's exact owned loop bounds + exchanges makes any stale halo unobservable in the owned
  dumps. **NO 1-rank regression** (step 65 + pressure 57 + advhor MR dist_2/8 + stepdyn MR dist_2/8 + 13/13 ctest all
  `max|Δ|=0`). See LESSONS **L34**.
- **MVP CLOSED 2026-06-22** — the whole-model multi-rank byte-match (M2.12) is committed and tagged `m2-mvp` (baseline
  re-confirmed after a clean rebuild: all 8 byte-gates `max|Δ|=0`/green, see the gate list in "RESUME HERE").
- **PRODUCTION VALIDATION (multi-rank FREE-RUNNING lifecycle) — ✅ DONE, BYTE-EXACT.** Built the multi-rank analog
  of the 1-rank CORE2 lifecycle (`fesom_lifecycle_mr` + `tools/run_lifecycle_gate_multirank.sh`, CORE2 dist_2, dt=1800,
  reduced-M2) — the FIRST test of *multi-rank × free-running × many-steps* (the m2-mvp gates are all single-step
  prescribe-and-stop, which MASK live-multi-step halo bugs). It found two latent bugs, **both now fixed**:
  1. ✅ **UV halo blow-up FIXED** — `update_vel` exchanged UV over `com_elem2D_full` (eDim+eXDim); FESOM2 uses
     `com_elem2D` (eDim). Wrong UV halo → viscosity → blow-up. Fix (`oce_ale.F90`): eDim `exchange_elem(UV)`. Model now
     STABLE + matches FESOM2's physical values to all printed digits.
  2. ✅ **`tr_xy` halo "drift" SOLVED — it was an OpenMPI `vader` KNEM single-copy bug, NOT a FESOM3 bug.** The
     exchange logic is correct; KNEM (the `vader` shared-memory single-copy path; CMA is blocked by
     `kernel.yama.ptrace_scope=3`) CORRUPTS messages >~131 KB on levante. A cross-rank buffer trace localized it to the
     MPI transport (`sbuf` packed correct → `rbuf` corrupt from byte 134656 → unpack faithful). The CORE2 `tr_xy`
     block exchange (478 KB) is the FIRST live F3 message above the threshold (every pi gate stays under it; FESOM2's
     `MPI_TYPE_INDEXED` path dodges KNEM) — which is why it only surfaced in the big-mesh free-running lifecycle.
     **Fix: `env.sh` exports `OMPI_MCA_btl_vader_single_copy_mechanism=none`** (forces byte-faithful copy-in/copy-out).
     Verified: CORE2 dist_2 lifecycle gate **195 records `max|Δ|=0`**; no 1-rank/pi regression. Full writeup +
     meta-lessons: [`HANDOFF-multirank-lifecycle-trxy.md`](HANDOFF-multirank-lifecycle-trxy.md) + LESSONS **L35**.
- **ACTIVE NEXT: M3 — sea ice (EVP).** The prescribed `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt`/`stress_surf`
  become the real air-sea budget (sea ice EVP + advection + thermo → `oce_fluxes`). **Scoped 2026-06-22** (oracle map +
  call-sequence pass — see "M3 entry notes" below). M3 mirrors M2: byte-gate each ice kernel against FESOM2 v2.7.3,
  `max|Δ|=0`. Standard EVP path only (`whichEVP=0`; NO icepack/meltponds/cavity/oasis-yac). Decomposed:
  - **M3a — ice foundation + cold-start IC — ✅ DONE (2026-06-22).** `max|Δ|=0` on ALL 4 fields (`ice_a_ice`,
    `ice_m_ice`, `ice_m_snow`, `ice_massmatrix` [870146 CSR entries]), CORE2 1-rank (`tools/run_ice_gate_core2.sh`).
    Rewrote `src/types/mod_ice.F90` (full `t_ice`: `data(1:3)`=a_ice/m_ice/m_snow + work EVP stress/strain + thermo +
    rheology params, ALL initialisers VERBATIM from `MOD_ICE.F90` incl. un-suffixed singles); new
    `src/ice/mod_ice_setup.F90` (`ice_setup`/`ice_allocate`/`ice_mass_matrix_fill`/`ice_initial_state`, optional-`partit`
    pattern). **`ice_mass_matrix_fill` rides the M2.6 `ssh_stiff` CSR** — instead of FESOM2's `nn_num`/`nn_pos` it scans
    `ssh_stiff%colind_loc(rowptr_loc(row):..)` (the SAME ordered neighbour list `init_stiff_mat_ale` fills), and it
    byte-matched → confirms the FESOM3 CSR column order == FESOM2 `nn_pos`. Oracle: NEW `port2/fesom2/src/fesom_ice_dump.F90`
    (env `FESOM_ICE_DUMP`, dumps after `ice_setup`, stops) wired in `fesom_module.F90` after `ice_setup`; runner
    `tools/run_icedump_core2.sh` (the forced-lifecycle init path, use_ice=.true.). NEW driver `src/drivers/fesom_icedump.F90`.
    Restart binary I/O for ice DEFERRED (M3 is cold-start; gates use mod_dump). **No regression** (13/13 ctest +
    step-65 `max|Δ|=0`). ⚠️ CORE2 `namelist.ice` (the oracle's effective config) DIFFERS from the `MOD_ICE.F90` type
    defaults: `ice_gamma_fct=0.5` (not 0.25), `ice_diff=0.0` (NO artificial ice diffusion, not 10.0), `albw=0.1` (not
    0.066); and the namelist parses literals as DOUBLES (so `cd_oce_ice=0.0055` etc. = `_WP` doubles, NOT single→WP) —
    bake the namelist doubles into the M3b driver config (these don't affect the SST-sign IC / geometry mass matrix).
  - **M3b — ocean2ice + EVP dynamics — ✅ DONE (2026-06-22).** `max|Δ|=0` on ALL 7 fields (`ice_srfoce_u`/`srfoce_v`
    [ocean2ice], `ice_uice`/`vice`, `ice_sigma11`/`12`/`22`), CORE2 1-rank, for **BOTH `whichEVP=0` (standard EVP) AND
    `whichEVP=1` (modified EVP / mEVP)** (`tools/run_evp_gate_core2.sh`; mEVP added per user request). New
    `src/ice/mod_ice_dyn.F90`: `ocean2ice` (srfoce_u/v/temp/salt/ssh from the ocean surface; node-avg of `UV(:,1,:)` +
    `exchange_nod`) + std-EVP `EVPdynamics` (`stress_tensor` strain→stress + `stress2rhs` ∇·σ + 120-subcycle velocity
    update [ocean drag + Coriolis + ∇·σ + wind-on-ice stress] + edge BC + per-subcycle `exchange_nod(u_ice/v_ice)`) +
    `EVPdynamics_m` (mEVP, Bouillon-2013 α/β pseudo-time iteration on `uice_aux`/`vice_aux`, inlined
    `stress_tensor_m`+`stress2rhs_m`+`ssh2rhs`, the `bc_index_nod2D` node BC) + `EVPdynamics_solve` dispatcher
    (whichEVP SELECT CASE; aEVP=2 deferred). Optional-`partit` pattern throughout (1-rank verbatim; exchanges guarded by
    `is_multirank`). Gate driver `src/drivers/fesom_evpdump.F90` (`FESOM3_WHICHEVP` env) prescribes the analytic ocean
    forcing (surface UV / hbar / stress_atmice — identical formulas to the oracle, from byte-identical rotated coords),
    runs ocean2ice + the EVP solve, dumps. Oracle: `port2/fesom2/src/fesom_ice_dump.F90::evp_dump_write` (env
    `FESOM_EVP_DUMP`, dispatches on namelist `whichEVP`, runs the REAL `ocean2ice`+`EVPdynamics`/`EVPdynamics_m`, dumps,
    stops; wired in `fesom_module.F90` after `ice_dump_write`) + `tools/run_evpdump_core2.sh [run] [out] [whichEVP]`.
    **M3a-deferred items built:** `bc_index_nod2D` (the mEVP node boundary mask) + `uice_aux`/`vice_aux` allocation in
    `mod_ice_setup` (`ice_allocate`/`build_bc_index_nod2D`). ⚠️ **`bc_index_nod2D` lives in `t_ice`, NOT `t_mesh`** —
    adding any component to `t_mesh` (which has type-bound `write(unformatted)`) trips an ifort `#7976` generic-
    serialization cascade across every `write_bin_array` call; `t_ice` has no type-bound I/O so it's safe. Value is
    byte-identical (0/1 mask). ⚠️ `cd_oce_ice`/`delta_min` set to the namelist DOUBLES (`0.0055_WP`/`1.0e-11_WP`) in the
    driver, NOT the `t_ice` single→WP defaults (the M3a precision note). mEVP uses `rdt=ice_dt` (NOT `ice_dt/steps`),
    `alpha_evp`/`beta_evp`=250 (exact). **No regression** (step 65 + M3a ice 4 + 13/13 ctest all `max|Δ|=0`).
  - **M3c — ice FCT advection — ✅ DONE (2026-06-22; first try).** `max|Δ|=0` on ALL 6 fields (`ice_rhs_a`/`rhs_m`/
    `rhs_ms` from `ice_TG_rhs` + post-advection `ice_a_ice`/`m_ice`/`m_snow` from `ice_fct_solve`), CORE2 1-rank, for
    **BOTH `whichEVP=0` AND `whichEVP=1`** feeding the advection (`tools/run_icefct_gate_core2.sh`). New
    `src/ice/mod_ice_fct.F90` (`ice_TG_rhs` Taylor–Galerkin rhs + `ice_fct_solve` = `ice_solve_high_order` 3-iter Jacobi
    → `ice_solve_low_order` → `ice_fem_fct` Zalesak limiter ×3), optional-`partit` pattern. **KEY: the neighbour-list
    basis is the ssh_stiff CSR (`rowptr_loc`/`colind_loc`), NOT FESOM3's own `nn_pos`** — M3a proved `fct_massmatrix`
    byte-matches FESOM2's `mass_matrix` position-by-position, so `colind_loc(clo:clo2)` IS FESOM2's `nn_pos(1:nn_num,row)`
    as an ordered list (incl. self), and `sum(fct_massmatrix(clo:clo2)*field(colind_loc(clo:clo2)))` reproduces
    `sum(mass_matrix*field(nn_pos))` byte-for-byte (same operands + reduction order). Driver
    `src/drivers/fesom_icefctdump.F90` (extends `fesom_evpdump`: ocean2ice + EVP → `uice/vice`, then `ice_TG_rhs` +
    `ice_fct_solve`); oracle `fesom_ice_dump.F90::ice_fct_dump_write` (env `FESOM_FCT_DUMP`); runners
    `run_icefctdump_core2.sh`/`run_icefct_gate_core2.sh`. ⚠️ `ice_diff=0.0`/`ice_gamma_fct=0.5` set by the driver to the
    CORE2 namelist values (NOT the `t_ice` 10.0/0.25 defaults; both exactly representable, no ULP subtlety). Oracle
    `ENABLE_OPENMP=OFF` ⇒ the `ice_fem_fct` `!$OMP ORDERED`/`ATOMIC` flux scatters run serially in element order = the
    FESOM3 serial path (accumulation byte-matches). `ice_fct_solve` is called as an external (not in
    `ice_fct_interfaces`); `ice_TG_rhs` IS. **No regression** (step 65 + M3a ice 4 + M3b evp 7×2 + M3c fct 6×2 + 13/13
    ctest all `max|Δ|=0`; the byte-exact `rhs` independently re-confirms M3b — the `uice/vice` feeding `ice_TG_rhs` were
    byte-identical). See LESSONS **L37**.
  - **M3d — ice thermodynamics — ✅ DONE (2026-06-23; first try, L38).** `max|Δ|=0` on ALL 6 fields (post-thermo
    `ice_a_ice`/`m_ice`/`m_snow` + `ice_flx_h`/`flx_fw` + `ice_t_skin`), CORE2 1-rank, for **BOTH `whichEVP=0` AND
    `whichEVP=1`** feeding the advection upstream (`tools/run_icethermo_gate_core2.sh`). New `src/ice/mod_ice_thermo.F90`
    (`cut_off` hmin/Armin clamp + `thermodynamics` per-node driver + `therm_ice` 0-layer Semtner/Hibler-1984 ice-class
    growth + `budget` 5-iter Newton-Raphson surface temp + `obudget` open-ocean growth+evaporation + `flooding` snow→ice
    + `TFrez` Millero-1978 + solar-zenith/albw helpers [ported for linkage, dead at `open_water_albedo=0`]),
    optional-`partit` (`exchange_nod(ustar)` guarded). **Atmospheric forcing + thermo flux diagnostics bundled into
    `t_atmflux`** (FESOM3 explicit-dataflow vs FESOM2's `g_forcing_arrays`/`g_forcing_param`/`o_arrays` globals); the gated
    outputs live in `t_ice`. Added `h_cutoff`/`hpdf` to `t_ice_thermo` (read unconditionally by `therm_ice`, dead at
    `new_iclasses=.false.`). Driver `src/drivers/fesom_icethermodump.F90` (extends `fesom_icefctdump`: + prescribed
    analytic atmosphere + `cut_off` + `thermodynamics`); oracle `fesom_ice_dump.F90::ice_thermo_dump_write` (env
    `FESOM_THERMO_DUMP`, wired in `fesom_module.F90` after `ice_fct_dump_write`). **KEY (L38): the `&ice_therm` params are
    read by the oracle as namelist DOUBLES (con=2.1656/consn/hmin/Armin=0.01/emiss=0.97/albedos/albw=0.1 [not 0.066]/
    Sice/h0/c_melt) — NOT the `t_ice_thermo` single→WP defaults; the driver re-sets them as `_WP` doubles. Likewise
    `Ch_atm_ice=Ce_atm_ice=0.00175_WP` (namelist.forcing doubles), `ref_sss_local=.true.` (⇒ rsss=S_oc),
    `use_virt_salt=.true.` (linfs), `l_snow=.true.`. cc=rhowat·4190 / cl=rhoice·3.34e5 are exactly representable
    (ice_init recompute = the single→WP default). All thermo arithmetic is SCALAR per-node ⇒ no L29 vectorised-divide
    trap.** No regression (step 65 + M3c fct 6×2 + 13/13 ctest all `max|Δ|=0`; M3c re-run since the shared oracle
    `fesom_ice_dump.F90` was edited).
  - **M3e — `oce_fluxes` coupling-out (THE PAYOFF) — ✅ DONE (2026-06-23; first try, L39).** `max|Δ|=0` on ALL 5 fields
    (`heat_flux`/`water_flux`/`virtual_salt`/`relax_salt` [nod2D] + `stress_surf` [2,elem2D]), CORE2 1-rank, for **BOTH
    `whichEVP=0` AND `whichEVP=1`** (`tools/run_iceflux_gate_core2.sh`). **These 5 fields ARE the proven M2.11c-2
    `fesom_flux_dump` set — the `fesom_lifecycle` PRESCRIBED M3-gap is now produced NATIVELY.** New
    `src/ice/mod_ice_oce_coupling.F90`: `oce_fluxes_mom` (ice-ocean drag blended with atm-ocean stress by ice
    concentration on nodes → averaged to `stress_surf` on elements) + `oce_fluxes` (the reduced air-sea budget:
    `heat_flux=-flx_h`, `water_flux=-flx_fw` then globally balanced, `virtual_salt`/`relax_salt` each globally balanced
    to zero net via `integrate_nod`/`ocean_area`) + a private `integrate_nod_2D` (the EXACT FESOM2 `gen_support.F90:318`
    sequential reduction, NOT `sum()`). **KEY (L39): `mesh%ocean_area` — the balancing divisor — was an un-gated landmine:**
    computed since M2.11a as `sum(area(1,1:nNodO))` but never consumed, so never gated; `oce_fluxes` is its first consumer,
    and `sum()` ≠ FESOM2's sequential `areasvol` loop is a 1-ULP trap that would drift every balanced field. Fixed
    `mod_mesh_areas.F90` to the faithful FESOM2 `oce_mesh.F90:2385` loop (geometry-neutral — `ocean_area`/`areawithcav`
    are the only changed scalars, no gate consumes them but M3e). `dens_flux` (MOC diagnostic, needs `sw_alpha`/`sw_beta`
    EOS) DEFERRED (out of M3 scope, not a surface BC). Driver `src/drivers/fesom_icefluxdump.F90` (extends
    `fesom_icethermodump`: + prescribed analytic `stress_atmoce_x/y` + `Ssurf`); oracle
    `fesom_ice_dump.F90::ice_flux_dump_write` (env `FESOM_OCEFLUX_DUMP`, after `ice_thermo_dump_write`; + `wr_r2` for the
    2-D `stress_surf`); runners `run_icefluxdump_core2.sh`/`run_iceflux_gate_core2.sh`. `t_atmflux` extended with the
    oce_fluxes I/O (`heat_flux`/`water_flux`/`virtual_salt`/`relax_salt`/`heat_flux_in`/`stress_node_surf`/
    `stress_atmoce_x/y`/`Ssurf`/`surf_relax_S`=1.929e-06). **No regression** (M3e 5×2 + M3d 6×2 + step 65 1-rank + step 65
    multirank dist_2 + pressure 57 + 13/13 ctest all `max|Δ|=0`).
  - **M3f — whole ice step + forced lifecycle (NATIVE fluxes) + multi-rank.** IN PROGRESS — M3f-1 + M3f-2 ✅ DONE
    (2026-06-23, first try, L40); M3f-3 (native CORE2 forcing read) + M3f-4 (multi-rank) remain.
    - **M3f-1 — `ice_timestep` assembly — ✅ DONE.** New `src/ice/mod_ice_step.F90::ice_timestep` = the faithful
      FESOM2 `ice_setup_step.F90:96` chain `EVPdynamics_solve → ice_TG_rhs → ice_fct_solve → cut_off → thermodynamics`
      (CMIP6 `dyngr*` / `h_ice` diagnostics + cavity cleans omitted — unconsumed in the reduced config). Refactored
      `fesom_icefluxdump` to call it; M3e gate re-run `max|Δ|=0` (byte-neutral, both whichEVP).
    - **M3f-2 — native-flux coupled multi-step forced lifecycle — ✅ DONE (first try).** New driver
      `src/drivers/fesom_lifecycle_native.F90` (= the proven M2.11c-2 ocean init + `step_oce`, but the prescribed
      `FESOM3_FLUX_FILE` is REPLACED by the native chain `ocean2ice → ice_timestep → oce_fluxes_mom → oce_fluxes`,
      runloop order ocean2ice→[atm]→ice_timestep→oce_fluxes_mom→oce_fluxes→step_oce). The per-step POST-bulk
      atmospheric forcing (16 nod2D arrays: shortwave/longwave/Tair/shum/prec_rain/prec_snow/runoff/u_wind/v_wind/
      Ch-Ce_atm_oce/stress_atmoce_x-y/stress_atmice_x-y/Ssurf) is PRESCRIBED from a NEW oracle dump
      (`port2/fesom2/src/fesom_atmflux_dump.F90`, env `FESOM_ATMFLUX_DUMP`, wired in `fesom_module.F90` next to
      `flux_dump_record`) so the native sea-ice + air-sea coupling is isolated from the CORE2 forcing read (M3f-3).
      **Result: `max|Δ|=0` on the 13 NODE substeps × 5 probes — 195 records (3 steps) AND 325 records (5 steps) — for
      BOTH whichEVP=0 (std EVP) AND whichEVP=1 (mEVP).** The driver's built-in per-step native-vs-oracle flux
      self-check (`FESOM3_FLUX_FILE`) reads `max|Δ(hf,wf,vs,rs,ss)|=0` at EVERY step. This is the FIRST test of the
      MULTI-STEP ice evolution (sigma elastic memory + t_skin carry-over + `values_old` — the M3a-e gates were all
      single-step from cold start) AND the FIRST native ice↔ocean coupled run. Gate
      `tools/run_lifecycle_native_gate_core2.sh [run] [nsteps] [whichEVP]` (oracle `run_lifecycle_forced_core2.sh`
      extended with `FESOM_ATMFLUX_DUMP` + a `whichEVP` namelist patch). **No regression** (step 65 1-rank + ctest 13/13
      + M2.11c-2 prescribed-flux forced lifecycle 195 + M3e flux 5×2 all `max|Δ|=0`). See LESSONS **L40**.
    - **M3f-3 — native CORE2 forcing read + bulk. M3f-3a + M3f-3c ✅ DONE 2026-06-23 (first try, L41); M3f-3b remains.**
      **14/16 atmospheric arrays are now NATIVE** (`max|Δ|=0`); only runoff + Ssurf (the monthly climatology) stay
      prescribed.
      - **M3f-3a — native CORE2 NCAR forcing READ — ✅ DONE.** The 8 NCAR fields (shortwave/longwave/Tair/shum/
        prec_rain/prec_snow/u_wind/v_wind from u_10/v_10/q_10/t_10/ncar_rad/ncar_precip) read + bilinear-interp +
        g2r-rotated via the M2.10a `mod_forcing_read` (already byte-proven on these EXACT files, on pi), now on the
        CORE2 mesh + the per-step `rdate` advance (clock `0 1 1948`, dt=1800, `sbc_do:1527` half-step;
        `rdate(n)=julday(1948,1,1,noleap)+(2n-1)·900/86400`; cold-start getcoeffld + per-step timeinterp, NO
        re-trigger within forcing day 1). Driver `src/drivers/fesom_forcing_core2.F90` self-checks the 8 fields vs the
        M3f-2 `atmflux_f2` dump: `max|Δ|=0`, 3 steps.
      - **M3f-3c — native NCAR bulk + wind stress — ✅ DONE.** `Ch_atm_oce_arr`/`Ce_atm_oce_arr` (`forcing_bulk_ncar`,
        M2.10b, fed the LIVE `ice%srfoce_temp/u/v` from ocean2ice — same inputs the oracle's
        `ncar_ocean_fluxes_mode` reads), `stress_atmoce_x/y` (`forcing_wind_stress`, M2.10b), and the NEW
        wind-on-ice `stress_atmice_x/y` (`forcing_ice_stress` in `mod_forcing_bulk`:
        `Cd_atm_ice·rhoair·|uw-uice|·(uw-uice)`, Cd_atm_ice=0.0012 const, fed the previous-step `ice%uice/vice`).
        Validated IN the lifecycle via a per-step self-check in `fesom_lifecycle_native` (env `FESOM3_FORCING_DIR`):
        the native NCAR read + bulk + 2 stresses are recomputed each step and compared to the prescribed atm —
        `max|Δ|=0` on all 14 fields, BOTH whichEVP, while the step itself still USES the prescribed values (so the
        195-record gate cannot regress). No regression: M2.10 forcing gate (pi, `mod_forcing_bulk` changed) PASS,
        ctest 13/13, the 195-record lifecycle gate MATCH.
      - **M3f-3b — runoff + SSS monthly climatology read (the last 2/16) — ✅ DONE (2026-06-23; first try, L42).**
        New `src/forcing/mod_forcing_other.F90` (`read_other_NetCDF` + `interp_2d_field`, faithful from FESOM2
        `gen_modules_read_NetCDF.F90:6` / `gen_interpolation.F90:145`) — a DISTINCT reader from `mod_forcing_read`
        (one 2D slice, raw-grid dummy fill, direct per-node bilinear, NO time-interp). runoff: `read_other_NetCDF
        ('Foxx_o_roff', rec 1, check_dummy=.false.)` → missing→0, then `/1000` (kg/s/m²→m/s); Ssurf: `('SALT', rec
        month=1, check_dummy=.true.)` → 30-neighbour fill. Both vertices (`do_onvert=.true.`; centroid path
        `error stop`-guarded). +2 real64 netCDF wrappers (`nc_get_slice_dp`/`nc_get_att_dp` in `mod_io_netcdf`). Wired
        into `fesom_lifecycle_native` (read ONCE at setup — both constant for a Jan run); self-check `max|d clim
        (runoff,Ssurf)|=0`. **THEN dropped the prescribed atm entirely → M3f-3 COMPLETE: the FULLY-NATIVE lifecycle**
        (`FESOM3_ATMFLUX_FILE` absent ⇒ compute all 16 in-driver via `apply_native_forcing`) is `max|Δ|=0` vs FESOM2,
        195 records (3-step) + 325 records (5-step), BOTH whichEVP, `tools/run_lifecycle_fullynative_gate_core2.sh`.
        ⚠️ TRAP (L42): the standalone 1-rank drivers do NOT populate `partit%myDim_nod2D` (it stays 0) — passing
        `partit` to an optional-`partit` routine whose loop bound is `partit%myDim_*` gives `num=0` → SILENT all-zero
        output (no crash). FIX: omit `partit` at 1-rank (count from `mesh%nod2D`), exactly like `ocean2ice`/
        `ice_timestep`. M3f-4 will pass partit (dims set by `read_dist_partition`).
    - **M3f-4 — multi-rank — ✅ DONE (2026-06-23; first try, L43) → M3 COMPLETE, tag `m3`.** `max|Δ|=0` on the per-rank
      gid-keyed NODE substeps — **195 records (3-step) AND 325 records (5-step), BOTH whichEVP=0 (std EVP) AND
      whichEVP=1 (mEVP), on CORE2 dist_2 AND dist_8** (`tools/run_lifecycle_fullynative_gate_multirank.sh [np] [nsteps]
      [whichEVP]`). The MULTI-RANK fully-native lifecycle (`src/drivers/fesom_lifecycle_native_mr.F90` = the
      `fesom_lifecycle_mr` MR mesh/state/dump scaffold ⊕ the `fesom_lifecycle_native` ice/atm/forcing setup + runloop,
      at LOCAL sizes, `, partit` on every kernel call) drives the whole coupled sea-ice + air-sea budget + atmosphere
      (16/16 native arrays) byte-exactly at multi-rank with NO prescribed input. **KEY (L43): the M3a–e kernels were
      ALREADY MR-ready** (each transcribed with the M2.12 optional-`partit` pattern; `ice_setup` sizes via `local_dims`,
      `ocean2ice`/EVP/FCT exchange their fields, `cut_off`+`thermodynamics` loop owned+halo so outputs stay halo-valid,
      `oce_fluxes` does the cross-rank `integrate_nod_2D`). The EVP needs **NO sigma element-halo** — `stress2rhs`
      scatters owned-elements→owned-nodes and the owned-node-completeness invariant (M2.12a) makes the owned rhs
      complete (FESOM2 `ice_EVP.F90` has only `exchange_nod(U_ice,V_ice)`). Only TWO things were genuinely missing:
      (1) `mesh%ocean_area`/`ocean_areawithcav` were a LOCAL-only owned sum → added `allreduce_sum` at npes>1
      (`mod_mesh_areas.F90`, FESOM2 `oce_mesh.F90:2389` `MPI_AllREDUCE`) — the M3e flux-balance divisor; (2) the 3
      `mod_forcing_bulk` routines looped `mesh%nod2D` (GLOBAL) → optional `partit` + `owned_bounds` → loop `nNodL`
      (FESOM2 computes the bulk over `myDim+eDim`, no exchange). The NCAR read needed only `frc%nnod=nNodL`;
      `read_other_NetCDF` gets `partit` at MR (1-rank omits it, the L42 trap). Oracle: `run_lifecycle_forced_core2.sh`
      now takes `[np]` (FESOM2 auto-reads dist_<np>; the 1-rank-only flux/atmflux dumps skip at npes/=1, the gid-keyed
      `dump_shim` is the gate). **No regression:** ctest 13/13 + 1-rank fully-native 195 + iceflux 5×2 + forcing pi +
      step-65 1-rank + step-65 MR dist_2 all `max|Δ|=0`. ⚠️ multi-rank levante needs the `env.sh` KNEM flag (L35).

**↳ RESUME HERE (next session): M5c — forced/fully-native CORE2 lifecycle + sw_pene.** Add the two
NEW tracer-TDMA terms to `diff_ver_part_impl_ale` (`oce_ale_tracer.F90`): the `ghats` nonlocal
counter-gradient flux (mix_scheme==1, uses `dyn%work%ghats` + `blmc` + heat_flux/water_flux,
oracle `:900-945`) + the `sw_3d` shortwave heating (use_sw_pene, tr ID==1, oracle `:991-996`); wire
`cal_shortwave_rad` (ALREADY byte-proven, M2.10c) into the native runloop to fill `dyn%work%sw_3d`
+ flip `use_sw_pene=.true.`. Gate the FORCED / fully-native lifecycle `max|Δ|=0`, BOTH whichEVP (the
production gate — exercises ghats + sw_3d heating + the full surface coupling). 195 + 325 records.
Then M5d (multi-rank — the kernels are optional-`partit` from the start; expect pure wiring + the
`smooth_nod` exchanges already in place, the M4f lesson). **M5a + M5b ✅ DONE byte-exact (see the M5
block below): M5b assembled `oce_mixing_kpp_driver` + dispatch and the UNFORCED lifecycle is
`max|Δ|=0` (KPP-alone AND KPP+GM+Redi). The KPP arrays live in `dyn%work` (allocated only when KPP);
`Av`/viscAE stays element-based so `impl_vert_visc` is UNCHANGED; `Kv=Kv_double(:,:,1)`.**

### M5 progress (KPP vertical mixing + sw_pene) — scoped + started 2026-06-24. Plan: `docs/plans/2026-06-24-m5-kpp.md`.
- **Decomposition:** M5a (producers, isolated gate; sub-stepped a1–a4 mirroring the C-port K1–K8) →
  M5b (wire into step + UNFORCED lifecycle) → M5c (forced/native + sw_pene + ghats) → M5d (multi-rank). Tag `m5`.
- **KEY STRUCTURAL FINDINGS (from reading the actual `oce_ale_mixing_kpp.F90` + the dispatch + the validated C port):**
  - `oce_ale.F90:3713` dispatch: `oce_mixing_KPP(Av, Kv_double, ...)` — **`Av` (viscAE) stays element-based** (KPP
    averages node→elem internally, `minmix=3e-3` floor) ⇒ `impl_vert_visc` UNCHANGED; then `Kv(:,n)=Kv_double(:,n,1)`
    (T-channel) + `mo_convect` ⇒ the tracer vertical-diffusion TDMA is **UNCHANGED** (same single `Kv` as PP). The
    S-channel `Kv_double(:,:,2)` is computed-but-dead in CORE2 (only a gated output). So KPP integration = the KPP
    module itself + `dbsfc` (in pressure_bv) + sw_pene wiring + TWO new tracer terms (ghats nonlocal + sw_3d).
  - `dbsfc` is filled in `oce_ale_pressure_bv.F90:332-339` (a NEW output of the already-ported pressure_bv) — in
    FESOM3 a `dyn%work%dbsfc` field; `ghats`/`blmc` are READ by the tracer solve ⇒ live in `dyn%work` (not KPP-locals).
  - `cal_shortwave_rad` (sw_3d producer) is ALREADY ported + byte-proven (M2.10c); M5 wires it into the lifecycle +
    adds the TDMA consumer (`oce_ale_tracer.F90:991-996` sw_3d term, `:900-945` ghats term). `ddmix` DEFERRED
    (`double_diffusion=.false.`). Config DOUBLES: `Ricr=0.3`, `concv=1.6`, `visc_sh_limit`/`diff_sh_limit=5e-3`,
    `A_ver=1e-4`, `K_ver=1e-5`, `Kv0_const=.true.`.
- **HARNESS (no oracle rebuild needed):** the oracle `build/lib/libfesom.so` (Jun-23 M4 build) ALREADY carries the
  `FESOM_KPP_DUMP_DIR` instrumentation (`kpp_dump_*` symbols; env-gated off by default → harmless). `bin/fesom.x`
  (Apr-23, stale) is a thin stub that loads that `.so` at runtime — the gates use `build/bin/fesom.x`. **NEW
  `tools/run_kppinit_oracle.sh`** runs CORE2 1-rank with `mix_scheme=KPP` (proven forced-runner setup, no PP
  downgrade) → emits `kpp_init_rank0.txt` + `kpp_wscale_rank0.txt` **AND all the per-node K3/K5/K6/K7/K8 dumps**
  (`kpp_dump_s1_{ri_*,prestep,dVsq,dbsfc,blmc_*,diffK*,ghats,viscA,viscAE,dkm1,bldepth}_rank0.txt`) — the M5a-2/3/4
  reference data is ALREADY GENERATED at `/scratch/a/a270088/kppinit_oracle/kpp_oracle/`. KPP+GM+Redi+linfs runs
  cleanly in the oracle at 1-rank.
- **M5a-1 ✅ DONE (byte-exact, first try):** NEW `src/oce/oce_mixing_kpp.F90` (`oce_mixing_kpp_init` builds Vtc/cg/
  deltaz/deltau + the wmt/wst lookup tables; `wscale`) + NEW driver `src/drivers/fesom_kppdump.F90` (emits the two
  txt files in the oracle's exact es24.16 format, Ricr=0.3/concv=1.6). FESOM3 vs oracle `kpp_init_rank0.txt`
  (429949 lines) + `kpp_wscale_rank0.txt` (20302 lines) are **IDENTICAL** (`diff` clean = `max|Δ|=0`; es24.16 =
  exact double round-trip). Validates the constant arithmetic, the `**1/3`/`1/4`/`1/2` table build (same Intel libm),
  + the wscale INT()/clamp/interp. NEW-FILE clean rebuild OK; no existing source touched (no regression possible).
- **M5a-2 ✅ DONE (byte-exact, first try):** `ri_iwmix` (interior Ri-mixing → `viscA`/`diffK` T&S) appended to the
  CORE2 pressure gate (option (b), the M2-M4 house method — RESOLVED the harness question). It reuses the pressure
  shim's already-prescribed strong-shear `uvnode` (FESOM3:556 == oracle:476, "Ri factor spans [0,~0.7]") + smoothed
  `bvfreq`. **`ri_viscA`/`ri_diffKt`/`ri_diffKs` all `max|Δ|=0`; all 70 pre-existing pressure fields still PASS;
  ctest 13/13** (`tools/run_pressure_gate_core2.sh` now 73 fields). Non-trivial: `ri_viscA` max=5.10e-3 (=visc_sh_limit
  +A_ver ⇒ frit=1 reached), `ri_diffKt/s` max=5.01e-3 — the full shear-instability shape range exercised. KEY: the L29
  divide is kept scalar via local POINTERS (`vA/dK` like `oce_mixing_pp`), and `AMAX1`/`AMIN1` kept VERBATIM (REAL
  intrinsics — Fortran↔Fortran is exact; the C port's 1e-9 drift was from mapping them to fmax). `ri_iwmix` loops
  OWNED nodes (`nNodO`, the oracle's `myDim_nod2D`); optional-`partit` from the start (M5d-ready). Changed: FESOM3
  `oce_mixing_kpp.F90` (+ri_iwmix), `mod_param_phys.F90` (+`diff_sh_limit`), `fesom_pressuredump.F90` (call+dump);
  oracle `oce_ale_mixing_kpp.F90` (`public ri_iwmix`), `fesom_pressure_dump.F90` (call+dump). Both rebuilt; the
  pressure gate's 70-field no-regression confirms the oracle rebuild didn't shift the dynamical core.
- **M5a-3 ✅ DONE (byte-exact, first try, 2026-06-24):** the KPP driver prestep (dVsq/ustar/Bo) + `dbsfc`(pressure_bv)
  + `bldepth` (the C-port HIGHEST-RISK routine: bulk-Ri accumulation + sw interp + ekman/monob limit) appended to BOTH
  pressure gates. **All 11 new fields `max|Δ|=0` on pi (3140) AND CORE2 (126858): `kpp_stress`/`kpp_sw3d` (inputs) +
  `kpp_dVsq`/`kpp_dbsfc`/`kpp_ustar`/`kpp_Bo` (prestep+dbsfc) + `kpp_hbl`/`kpp_kbl`/`kpp_bfsfc`/`kpp_stable`/`kpp_caseA`
  (bldepth outputs).** NON-VACUOUS: hbl spans 0.13→1180 m, BOTH forcing branches fire (57485 stable / 69373 unstable
  nodes on CORE2). Implementation:
  - **`bldepth`** ported into `src/oce/oce_mixing_kpp.F90` (explicit-dataflow args, optional-`partit`, owned-node loops;
    `smooth_hbl=.false.` ⇒ no exchange). All divides are scalar (inner nz-loop EXITs at the Rib crossing + calls
    `wscale`) ⇒ NO L29 SIMD trap; `SIGN`/`AMIN1` kept VERBATIM. Calls the M5a-1 `wscale`/tables.
  - **`dbsfc`** added as an OPTIONAL output of `pressure_bv` (`oce_pressure_bv.F90`): a separate `!DIR$ NOVECTOR` loop
    (scalar divide, matching the oracle's db_max-reduction density loop), guarded by `present(dbsfc)` so the
    lifecycle/step callers (dbsfc absent) are byte-neutral. The gate driver passes it on the smoothed `pressure_bv`
    call. Faithful to FESOM2 `oce_ale_pressure_bv.F90:326-339` (`-g*(rho_surf-rho_full)/rho_full`).
  - **Harness (the M5a-2 method, option (b)):** BOTH shims prescribe surface fluxes (`heat_flux`=200·.../`water_flux`/
    `stress_node_surf`) + a physical decaying `sw_3d` (exp(zbar/20), surface ~2.4e-5 K m/s) ANALYTICALLY (byte-identical
    formula, M3b style; NOT `cal_shortwave_rad` — that wiring is M5c) so Bo spans both signs. Compute dVsq/ustar/Bo
    inline (the oce_mixing_KPP prestep, verbatim both sides), then call `bldepth`. Oracle: `bldepth`+`dVsq`/`ustar`/
    `bfsfc`/`stable`/`caseA`/`kbl` made PUBLIC; `dbsfc` filled by enabling `mix_scheme_nmb=1` + `oce_mixing_kpp_init`
    (allocates the KPP arrays) around the smoothed `pressure_bv`; `sw_3d` allocated in-shim (forcing-init is later).
  - ⚠️ **META-LESSON (uninit-memory gate fragility):** the oracle shim allocated `fer_tapfac` (M4a) WITHOUT init →
    its below-bottom region was uninitialised garbage that COINCIDENTALLY matched FESOM3's 0.0. My new
    `oce_mixing_kpp_init` (allocated ~11 arrays before the `fer_tapfac` alloc) shifted the heap → `fer_tapfac`
    below-bottom became 1.0 on pi (CORE2 still 0.0). FIX: `fer_tapfac = 0.0_WP` after the alloc (matches FESOM3,
    deterministic). The consumed fields (`slope_tapered`/`neutral_slope`/`sigma_xy`) all stayed `max|Δ|=0` — it was
    a dump-only below-bottom artifact, not a physics bug. **Lesson: any shim array that's dumped+gated must be
    explicitly initialised below-bottom (don't rely on fresh-alloc memory); a later allocation can shift the heap.**
  - **No regression:** pi pressure gate (now incl. fer_tapfac fix) + CORE2 pressure gate (84 fields) + 1-rank step-65
    + CORE2 lifecycle 195 + ctest 13/13 all `max|Δ|=0`/green. FESOM3 changed: `oce_mixing_kpp.F90` (+bldepth),
    `oce_pressure_bv.F90` (+dbsfc), `fesom_pressuredump.F90` (M5a-3 block). Oracle (uncommitted): `oce_ale_mixing_kpp.F90`
    (public bldepth+arrays), `fesom_pressure_dump.F90` (M5a-3 block + fer_tapfac init) + `libfesom.so` rebuilt.
- **M5a-4 ✅ DONE (byte-exact, first try, 2026-06-24) → M5a COMPLETE (isolated KPP module gate closed).** `blmix_kpp`
  (BL mixing coeffs `blmc(3)` + `dkm1(3)` + nonlocal `ghats`, eqn 10/11/20: T uses diffK(1)→blmc(2), S diffK(2)→blmc(3),
  mom diffK→blmc(1)) + `enhance` (kbl-1 interface blend) + the driver tail (combine: within-BL `max(.,blmc)`, outside
  `ghats=0`; node→elem viscAE average + `minmix=3e-3` surface floor) appended to BOTH pressure gates. **All 11 new
  fields `max|Δ|=0` on pi (3140) AND CORE2 (126858): `kpp_blmc1/2/3` + `kpp_dkm1m/t/s` + `kpp_ghats` + final
  `kpp_viscA`/`kpp_diffKt`/`kpp_diffKs` + `kpp_viscAE`** (the pressure gates now total **95 fields** each). Implementation:
  - **`blmix_kpp` + `enhance`** ported into `src/oce/oce_mixing_kpp.F90` (explicit-dataflow args — `blmc`/`ghats`/`dkm1`
    passed in, vs the oracle's module arrays — optional-`partit`, owned-node loops; `blmc` pre-zeroed over owned+halo
    like the oracle). They read the M5a-2 `viscA_ri`/`diffK_ri` (interior) + the M5a-3 `hbl`/`kbl`/`bfsfc`/`stable`/
    `caseA`/`ustar` — all byte-proven, so blmix consumes them directly. **NO L29 NOVECTOR needed** (first try): every
    blmix divide (`/dthick`, `/(hbl+epsln)`, `/(wm+epsln)`, `/(ws+epsln)`) is SCALAR per-node — the inner nz-loop EXITs
    at kbl + calls `wscale`, and enhance's `delta` is one scalar divide; `AMIN1`/`MIN`/`INT`/`ABS`/`SIGN` kept VERBATIM.
  - **Combine + node→elem viscAE average** done INLINE in BOTH shims (the M5a-3 prestep pattern), operating on SEPARATE
    copies `viscA_fin`/`diffK_fin` so the M5a-2 `ri_viscA`/`ri_diffKt`/`ri_diffKs` dumps still echo the interior (the
    combine overwrites viscA/diffK in place in the real driver; the copies isolate the two gate points). `viscAE` is the
    `pp_Av`-style element field (`SUM(viscA_fin(nz,elnodes))/3`, gid-order byte-exact at 1-rank); `kpp_ghats` is the
    post-combine nonlocal flux (== post-enhance here: blmix only writes nz<kbl, combine only zeros nz≥kbl which were
    already 0). Oracle: `blmix_kpp`/`enhance`/`dkm1` made PUBLIC (like `bldepth`); the shim calls them on `viscA_ri`/
    `diffK_ri` (`blmc`/`ghats`/`dkm1` already allocated+zeroed by the M5a-3 `oce_mixing_kpp_init`).
  - **No regression:** pi pressure gate (95) + CORE2 pressure gate (95) + ctest 13/13 + step-65 1-rank (65) all
    `max|Δ|=0`/green; the heap-shift-sensitive `tsol_*` + `fer_tapfac` (L:M5a-3) stayed byte-exact on BOTH meshes.
    FESOM3 changed: `oce_mixing_kpp.F90` (+blmix_kpp/+enhance), `fesom_pressuredump.F90` (M5a-4 block + 11 dumps).
    Oracle (uncommitted): `oce_ale_mixing_kpp.F90` (public blmix_kpp/enhance/dkm1), `fesom_pressure_dump.F90` (M5a-4
    block + 11 dumps) + `libfesom.so` rebuilt.
- **M5b ✅ DONE (2026-06-24) — KPP wired into `step_oce`; UNFORCED CORE2 lifecycle byte-exact.** `max|Δ|=0` KPP-alone
  195 (3-step) + 325 (5-step) AND KPP+GM+Redi 195 (`tools/run_lifecycle_kpp_gate_core2.sh [run] [nsteps]`, env
  `FER_GM=1`/`REDI=1`). NEW `oce_mixing_kpp_driver` (renamed from `oce_mixing_KPP` to dodge the case-insensitive clash
  with the module name) + KPP `dyn%work` arrays (allocated only when KPP, not serialized) + `mix_scheme_nmb==1` dispatch
  in `mod_step_oce` (pressure_bv `dbsfc=`; `sw_alpha_beta` fires for `Fer_GM.or.Redi.or.is_kpp`; `Kv=Kv_double(:,:,1)`
  over `nNodL` then `mo_convect`; `Av`/viscAE element-based → `impl_vert_visc` UNCHANGED). `step_oce` gained an OPTIONAL
  `stress_node_surf(2,nod)` (KPP ustar; unforced⇒0, passed unallocated for PP⇒seen absent). The reduced-M2 drivers now
  pin `mix_scheme_nmb=2` (the new dispatch defaults to KPP=1). `FESOM3_MIX_KPP` env in `fesom_lifecycle` (mix_scheme_nmb=1
  + `use_sw_pene=.false.` + work_core DOUBLES + alloc + `oce_mixing_kpp_init`); oracle `run_lifecycle_core2.sh` `MIX_KPP=1`
  KEEPS `mix_scheme='KPP'` (NO oracle source change). **ROOT-CAUSE (L45): the ONE step beyond M5a-4 was
  `smooth_blmc=.true.`** (oracle `:439-449`: 3 area-weighted `smooth_nod` sweeps × the 3 `blmc` channels, BEFORE the
  combine) — the M5a-4 isolated gate called `blmix`/`enhance`/combine DIRECTLY so never exercised it; omitting it diverged
  ~1e-3 in Kv at the kbl-1 BL level. Fixed by reusing the M2.1 byte-proven `smooth_nod` (made public in `oce_pressure_bv`;
  == FESOM2 `gen_support.F90` `smooth_nod3D`). **No regression:** PP step-65 1-rank + dist_2 + PP lifecycle 195 +
  CORE2/pi pressure (95) + ctest 13/13 all `max|Δ|=0`. *(Historical M5b plan follows.)*
- **M5b PLAN (done above):** the isolated KPP module is fully byte-proven — now ASSEMBLE the driver + wire it into the step.
  Build the real `oce_mixing_KPP` driver in `oce_mixing_kpp.F90` (prestep dVsq/ustar/Bo + `ri_iwmix` → `bldepth` →
  `blmix_kpp` → `enhance` → combine → node→elem viscAE average — all four sub-kernels already byte-proven, so this is
  assembly, not new physics) taking `viscAE`(elem)+`Kv_double`(node,ntr); add KPP arrays (`Kv_double`/`ghats`/`blmc`/
  `hbl`/`kbl`/`dbsfc`/`dVsq`/`ustar`/`Bo`/`bfsfc`/`stable`/`caseA`/`dkm1`) to `dyn%work`, allocated only when KPP.
  Dispatch in `mod_step_oce` (`oce_ale.F90:3711`): after the producers call `oce_mixing_KPP`, then `Kv(:,n)=
  Kv_double(:,n,1)` (single T-channel Kv → the tracer TDMA is UNCHANGED), then `mo_convect` (already ported). `Av`=
  viscAE stays element-based → `impl_vert_visc` UNCHANGED. `dbsfc` is already the optional pressure_bv output (M5a-3);
  wire it through `dyn%work`. Gate the UNFORCED CORE2 lifecycle `max|Δ|=0` (unforced ⇒ zero surface flux + sw ⇒
  ghats=0, sw_3d=0 ⇒ tests the KPP Kv/Av integration + mo_convect-after-KPP). M5c adds the ghats + sw_3d tracer terms
  (forced/native); M5d multi-rank (the kernels are optional-`partit` from the start → pure wiring, the M4f lesson).

**M4 (GM/Redi) ✅ COMPLETE + COMMITTED (tag `m4`).** Byte-exact `max|Δ|=0` 1-rank AND multi-rank (CORE2 dist_2/dist_8),
unforced AND forced/fully-native, GM-only AND GM+Redi, BOTH whichEVP. Gates: `run_lifecycle_gm_gate_core2.sh` (M4c),
`run_lifecycle_redi_gate_core2.sh` (M4d), `run_lifecycle_gmredi_native_gate_core2.sh` (M4e),
`run_lifecycle_gmredi_native_gate_multirank.sh` (M4f). M4 enabled by env: FESOM3 drivers read `FESOM3_FER_GM`/
`FESOM3_REDI`; oracle runners read `FER_GM=1`/`REDI=1` to KEEP work_core `Fer_GM`/`Redi=.true.`.

**M4f KEY FINDING (the M4 capstone):** the M4 routines were ALREADY MULTI-RANK READY — every kernel was transcribed
(M4a–M4d) with the M2.12/M3 optional-`partit` pattern (owned-loop bounds + the FESOM2 halo exchanges), and
`mod_step_oce::step_oce` already threads the optional `partit` to all of them. So M4f was a pure WIRING task (the
GM/Redi config + LOCAL-sized array allocation into `fesom_lifecycle_native_mr` + the gate env) — the partit-present
path that M4a–M4d wrote but never gated just worked, FIRST try, NO halo extension, NO `dist_N` invariant check needed
(the GM/Redi producers reuse the SAME mesh ops — `nod_in_elem2D`+`gradient_sca` over owned elements, owned-edge loops —
the M2.12 dynamics already proved invariant). `sw_alpha_beta`/`tracer_gradient_z` compute over owned+halo (no exchange);
`compute_sigma_xy`/`compute_neutral_slope`/`init_Redi_GM`/`fer_solve_Gamma`/`fer_gamma2vel`/`vert_vel_ale` exchange
their outputs; the Redi diff terms loop owned edges/nodes reading exchanged/halo-valid inputs.

**M5 = KPP vertical mixing (the real production CORE2 config: `mix_scheme='KPP'` instead of the reduced-M2 PP).**
The production `work_core` runs KPP + Fer_GM + Redi + sw_pene + `which_ALE='zlevel'`. M4 added GM+Redi; M5 swaps PP→KPP
and turns on `use_sw_pene` (the shortwave penetration term in the tracer TDMA, currently `.false.` to dodge the
unallocated `sw_3d`). Mirror the M4 method: port `oce_ale_mixing_kpp.F90` (the KPP boundary-layer + interior scheme),
gate `max|Δ|=0` vs FESOM2 with `mix_scheme='KPP'`, then the forced/fully-native lifecycle, then multi-rank. The oracle
has work-dirs staged for the future targets (`work_zstar_kpp`/`work_tke_dump`). See "Milestone ordering" below + the
plan doc `2026-06-18-fesom3-architecture.md`.

⚠️ multi-rank levante needs the `env.sh` KNEM flag (L35). ⚠️ Gate every consumed intermediate (L29). ⚠️ **Read the
ACTUAL FESOM2 `.F90` + `work_core` namelists, NOT the plan summaries** (they were wrong on `MLD1_ind` + `K_hor` in M4).
M3 ✅ COMMITTED (`b591153`, tag `m3`); M4 ✅ COMMITTED (`7a088de`, tag `m4`); the M2 baseline stays at `m2-mvp`.
**M5a+M5b ✅ COMMITTED (`30760c3`, UNTAGGED — M5c forced/native + sw_pene and M5d multi-rank remain before tag `m5`).**

**Milestone ordering (plan doc `2026-06-18-fesom3-architecture.md`):** M4 = GM/Redi → M5 = KPP + production
multi-year (paper-parity) → M6 = beyond-paper (zstar/zlevel ALE, TKE, aEVP). mEVP (whichEVP=1) was already
delivered in M3. **The real production CORE2 config (`work_core`) runs KPP + Fer_GM + Redi + sw_pene + which_ALE='zlevel'**
— the reduced-M2 gate sed-downgrades ALL of these (linfs/PP/no-GM/no-Redi/no-sw_pene, `run_lifecycle_forced_core2.sh:44-53`);
the oracle has work-dirs staged for the future targets (`work_zstar_dump`/`work_zstar_kpp`/`work_zstar_tke`/`work_tke_dump`).

**M4 = GM/Redi (turn ON `Fer_GM`+`Redi` over the existing PP+linfs core; gate `max|Δ|=0` vs FESOM2 with the matching
config).** Decomposed M4a–M4f (full detail + oracle map + SIMD-divide watch-list in the plan file):
- **M4a ✅ DONE** — producers (`mesh_resolution` compute, `sw_alpha_beta`, `compute_sigma_xy`,
  `compute_neutral_slope`→`slope_tapered`/`fer_tapfac`); feed-forward, the safest first gate. (`MLD1_ind` turned out
  DEAD under work_core — see M4b.)
- **M4b ✅ DONE** — GM diffusivity + streamfunction + bolus velocity (new `oce_fer_gm.F90`: `init_Redi_GM`,
  `fer_solve_Gamma` TDMA, `fer_gamma2vel`; `fer_w` in `vert_vel_ale`). `Fer_GM=T`, `Redi=F`. 70-field pressure gate
  `max|Δ|=0`, first try.
- **M4c ✅ DONE** — GM bolus into tracer advection (add/subtract `fer_uv`/`fer_w` in `solve_tracers_ale` + GM chain
  wired into `step_oce`); unforced lifecycle gate `max|Δ|=0` (195+325 records). `Fer_GM=T`, `Redi=F`.
- **M4d ✅ DONE** — Redi isopycnal diffusion (`Ki` real field; K13/K23 in `diff_part_hor_redi`, NEW
  `diff_ver_part_redi_expl` K31/K32, K33 in `diff_ver_part_impl_ale`). Gated `Fer_GM=T, Redi=T` (the plan's `Fer_GM=F`
  is vacuous — work_core `K_hor=0`); unforced lifecycle `max|Δ|=0` (195+325). Merges with M4e-unforced.
- **M4e ✅ DONE** — both on in the FORCED / fully-native lifecycle (production tracer physics, + ice + native forcing).
  `Fer_GM=T`, `Redi=T`. Byte-exact 195+325 records, BOTH whichEVP (`tools/run_lifecycle_gmredi_native_gate_core2.sh`).
- **M4f ✅ DONE** — multi-rank (the M4 routines were ALREADY MR-ready; pure wiring + gating). CORE2 dist_2/dist_8,
  195+325 records, BOTH whichEVP, `max|Δ|=0` (`tools/run_lifecycle_gmredi_native_gate_multirank.sh`). Tag `m4`.

**M4 ✅ COMPLETE + COMMITTED (tag `m4`). M4a producers + M4b GM diffusivity/streamfunction/bolus-velocity + M4c
GM-bolus-into-advection + M4d Redi isopycnal diffusion + M4e GM+Redi in the FORCED/fully-native lifecycle + M4f
multi-rank ✅ DONE + byte-gated 2026-06-23/24 (`max|Δ|=0`, CORE2 1-rank AND dist_2/dist_8). ACTIVE NEXT = M5 (KPP).**
- **M4a** — `mesh_resolution` (`mod_mesh_areas.F90::compute_mesh_resolution`, rides the geom gate → 20 fields) +
  `sw_alpha_beta`/`compute_sigma_xy`/`compute_neutral_slope` (`oce_pressure_bv.F90`): the 6 fields
  `sw_alpha`/`sw_beta`/`sigma_xy`/`neutral_slope`/`slope_tapered`/`fer_tapfac` appended to the CORE2 pressure gate.
  Both harnesses force the `work_core` taper config around the producer calls (`Fer_GM=Redi=Redi_Ktaper=.true.`,
  `ODM95_Scr=0.2e-2` — the work_core override, not the 1e-2 default), then restore; the oracle snapshots
  `slope_tapered` before the shim zeros it at `:483`.
- **M4b** — new `src/oce/oce_fer_gm.F90` (`init_Redi_GM` + `fer_solve_Gamma` TDMA + `fer_gamma2vel`, explicit-dataflow
  + optional-`partit`) + `fer_w` added to `vert_vel_ale` (`oce_ale.F90`, `if(Fer_GM)`-guarded → byte-neutral at
  `Fer_GM=.false.`). The 6 fields `fer_K`/`fer_c`/`fer_scal`/`fer_gamma`/`fer_uv`/`fer_w` appended to the CORE2
  pressure gate → **70 PASS / 0 FAIL** (FIRST try, NO NOVECTOR needed). FESOM3 `fesom_pressuredump` + oracle
  `fesom_pressure_dump` each force `Fer_GM=T, Redi=F` + work_core GM config (`K_GM_max=1000`/`min=2`/`cm=3`/`cmin=0.1`/
  `resscalorder=2`/`rampmax=rampmin=-1`/`Ktaper=F`/`scaling_resolution=T`/`scaling_GMzexp=T`/`zref=500`/`smin=0.6`)
  around the GM calls + the existing `vert_vel_ale` (so it fills `fer_w`), then restore; fer arrays init'd to the
  FESOM2 `oce_setup_step.F90:962-967` values (`fer_K=500`/`fer_c=1`/`fer_scal=0`/`fer_gamma=0`, `fer_uv/fer_w=0`) so
  below-bottom byte-matches. **KEY: `MLD1_ind` is DEAD under work_core** (read ONLY in `init_Redi_GM`'s
  `scaling_Ferreira` branch; work_core uses `scaling_GMzexp` → GMzexp `exp(-|zbar_3d_n|/zref)`, needs no MLD1_ind/
  bvref) → NOT implemented (L33, the plan's "fold into M4b" assumed Ferreira); `fer_tapfac`/`neutral_slope` likewise
  unconsumed by GM (`K_GM_Ktaper=F`, `scaling_FESOM14=F`). Only the work_core GM path is ported — `init_Redi_GM`
  `error stop`s on `scaling_Ferreira/Rossby/GINsea/FESOM14/K_GM_Ktaper/Redi` (the **Redi Ki path is M4d**, which
  removes Redi from the guard + passes a `Ki` arg). `zbar_n_bot(n)` (absent in FESOM3's mesh) ==
  `mesh%zbar_3d_n(nlevels_nod2D(n),n)` (FESOM2 `oce_ale.F90:550`); `fer_solve_Gamma` mirrors the oracle's
  `tr => fer_gamma(:,:,n)` POINTER (dodges the L29 SIMD-divide trap). No regression (pressure 70 + step-65 1-rank +
  step-65 MR dist_2 + ctest 13/13 `max|Δ|=0`).
- **M4c** — GM bolus wired INTO the step + lifecycle (`mod_step_oce`: producers `sw_alpha_beta`→`compute_sigma_xy`
  after PGF, then `init_Redi_GM`→`fer_solve_Gamma`→`fer_gamma2vel` after `update_eta_n` before `vert_vel_ale`; bolus
  add/subtract in `solve_tracers_ale` — once around the tracer loop over owned+halo, `uv += fer_uv`/`w_e += fer_w`/
  `w += fer_w` before advection, subtracted before the salinity clamp; all `if(Fer_GM)`). GM aux fields added to
  `t_dyn_work` (non-serialized; allocated only when `Fer_GM`). The UNFORCED CORE2 lifecycle with GM ON is **`max|Δ|=0`
  — 195 records (3-step) AND 325 records (5-step)** (`tools/run_lifecycle_gm_gate_core2.sh`, FIRST try). Gate: the
  oracle base namelist IS work_core (GM params already correct) so `run_lifecycle_core2.sh` gained `FER_GM=1` to KEEP
  `Fer_GM=.true.` (Redi sed off); FESOM3 `fesom_lifecycle` gained `FESOM3_FER_GM` (sets `Fer_GM` + work_core GM config
  + allocates/inits the GM arrays). NO oracle source change. `compute_neutral_slope` SKIPPED (Redi/Ktaper-only,
  unconsumed at M4c → M4d). `fer_w` into BOTH `w` and `w_e` (faithful; `w_e` unused at `use_wsplit=.false.`). No
  regression: GM-OFF lifecycle 195 `|Δ|=0` (the `Fer_GM`-guarded edits are byte-neutral), step-65 1-rank `|Δ|=0`,
  ctest 13/13.
- **M4d** — Redi isopycnal diffusion, gated `Fer_GM=T + Redi=T` (production GM+Redi) UNFORCED lifecycle: **`max|Δ|=0`,
  195 (3-step) + 325 (5-step) records** (`tools/run_lifecycle_redi_gate_core2.sh`, FIRST try, NO NOVECTOR). **KEY (from
  the actual Fortran): the plan's "Fer_GM=F, Redi=T" isolating config is VACUOUS** — work_core `namelist.tra` sets
  `K_hor=0` so `init_Redi_GM` F1 `Ki(nzmin)=K_hor**(..)=0`; the only non-zero `Ki` is the `if(Redi.and.Fer_GM)` coupling
  `Ki=max(fer_scal*Redi_Kmax,K_GM_min)` (needs Fer_GM=T). So M4d gates BOTH on (merges plan M4d+M4e-unforced); M4c
  already pinned GM-only so a pass isolates the Redi terms. Implemented: `init_Redi_GM` Redi path (Redi removed from the
  guard; `Redi_Kmax<=0` sync; F1/coupling/F2 Ki + `Redi_Ktaper` sqrt-split; `Ki`/`fer_tapfac` optional args);
  `compute_neutral_slope` into `step_oce` (`if(Fer_GM.or.Redi)`); `diff_part_hor_redi` K13/K23 (`Tz`/`SxTz`/`SyTz`,
  `replace_all` across the 5 level-ranges); NEW `diff_ver_part_redi_expl` K31/K32; `diff_ver_part_impl_ale` K33
  (`slope_tapered(3)^2*Ki` augments `Kv` in the tridiag); `tracer_gradient_z`→`tr_z` (new `t_tracer_work` field).
  Redi aux (`Ki`/`neutral_slope`/`slope_tapered`/`fer_tapfac`) in `t_dyn_work`, allocated only when Redi; ALL Redi reads
  `if(Redi)`-guarded so the off configs (arrays unallocated) are byte-neutral (`-(Kv+0)*..==-Kv*..` by IEEE +0). Gate:
  `run_lifecycle_core2.sh` `REDI=1` keeps work_core `Redi=.true.`; `fesom_lifecycle` `FESOM3_REDI` enables it. NO oracle
  source change; Redi non-vacuous (max|uv| differs from GM-only). No regression: M4c GM-only 195 `|Δ|=0`, GM-OFF 195
  `|Δ|=0`, step-65 `|Δ|=0`, ctest 13/13.
- **M4e** — GM+Redi in the FORCED / fully-native lifecycle (production tracer physics + ice + native forcing): **`max|Δ|=0`,
  195 (3-step) + 325 (5-step) records, BOTH whichEVP=0 (std EVP) AND whichEVP=1 (mEVP)**
  (`tools/run_lifecycle_gmredi_native_gate_core2.sh`, FIRST try, 2026-06-24). **Pure WIRING — NO new physics code.** The
  M4a–M4d kernels (`step_oce`'s `Fer_GM`/`Redi`-guarded chain) are unchanged; M4e just enables them in the native run:
  (1) FESOM3 `fesom_lifecycle_native.F90` got the `FESOM3_FER_GM`/`FESOM3_REDI` env reads + work_core GM+Redi config +
  `dyn%work` GM/Redi array allocation block, COPIED VERBATIM from `fesom_lifecycle.F90`; `step_oce` (1-rank, no `partit`)
  reads the module flags + consumes the `dyn%work`/`fer_uv`/`tr_z` arrays → the GM+Redi chain runs on the FORCED state
  with no further change (first GM/Redi exercise on a non-zero-flux ocean). (2) Oracle `run_lifecycle_forced_core2.sh`
  gained `FER_GM`/`REDI` env (mirror `run_lifecycle_core2.sh`) — `=1` KEEPS work_core `Fer_GM`/`Redi=.true.`; NO oracle
  source change. (3) `run_lifecycle_fullynative_gate_core2.sh` gained `FER_GM`/`REDI` env (passes to oracle + exports
  `FESOM3_FER_GM`/`FESOM3_REDI`); thin wrapper `run_lifecycle_gmredi_native_gate_core2.sh` (FER_GM=1 REDI=1). The native
  forcing + flux self-checks print `max|Δ|=0` at every step under GM+Redi. No regression: M3f GM/Redi-OFF fully-native
  195 `|Δ|=0` (the shared-oracle-runner + gate-script edits are byte-neutral OFF), M4d unforced GM+Redi 195 `|Δ|=0`,
  ctest 13/13. Build incremental (only `fesom_lifecycle_native` relinked — no core-physics/new-file change).
- **M4f** — multi-rank GM+Redi (production tracer physics at MULTI-RANK): **`max|Δ|=0`, CORE2 dist_2 AND dist_8,
  195 (3-step) + 325 (5-step) records, BOTH whichEVP** (`tools/run_lifecycle_gmredi_native_gate_multirank.sh`,
  FIRST try, 2026-06-24). **Pure WIRING — NO new physics, NO halo extension, NO `dist_N` invariant check needed.** The
  M4 routines were ALREADY MR-ready (transcribed M4a–M4d with the M2.12/M3 optional-`partit` pattern; `step_oce` threads
  `partit` to all of them — the partit-present path was written but never gated). The exchanges already in place:
  `compute_sigma_xy`→`exchange_nod(sigma_xy)` (the rank-3 superset of the oracle's 3×`MPI_BARRIER` per-component dance =
  same bytes), `compute_neutral_slope`→`exchange_nod(neutral_slope/slope_tapered)`, `init_Redi_GM`→`exchange_nod(fer_c/
  fer_K/Ki)`, `fer_solve_Gamma`→`exchange_nod(fer_gamma)`, `fer_gamma2vel`→`exchange_elem(fer_uv)`, `vert_vel_ale`→
  `exchange_nod(fer_w)`. `sw_alpha_beta`+`tracer_gradient_z` compute owned+halo (no exchange — the compute-at-halo idiom).
  The Redi diff terms loop OWNED edges/nodes reading exchanged `slope_tapered`/`Ki` + halo-valid `tr_z` (reusing the
  M2.12c-3 owned-edge invariant). Only edit: the GM/Redi config + LOCAL-sized (nNodL/nElemF) array allocation block into
  `fesom_lifecycle_native_mr.F90` (`fer_uv` sized nElemF like `dyn%uv`) + the gate env. No regression: GM/Redi-OFF MR
  fully-native dist_2 195 `|Δ|=0`, ctest 13/13.
**M4 ✅ COMMITTED (tag `m4`, M4a–M4f together).** FESOM3 changed: `oce_fer_gm.F90` (NEW)/`oce_ale.F90`/
`oce_ale_tracer.F90`/`oce_pressure_bv.F90`/`oce_tracer_grad.F90`/`mod_dyn.F90`/`mod_tracer.F90`/`mod_step_oce.F90`/
`mod_mesh_areas.F90`/`mod_geom_dump.F90`/`fesom_pressuredump.F90`/`fesom_lifecycle.F90`/`fesom_lifecycle_native.F90`/
`fesom_lifecycle_native_mr.F90` + `tools/run_lifecycle_core2.sh`/`run_lifecycle_gm_gate_core2.sh`/
`run_lifecycle_redi_gate_core2.sh`/`run_lifecycle_forced_core2.sh`/`run_lifecycle_fullynative_gate_core2.sh`/
`run_lifecycle_fullynative_gate_multirank.sh`/`run_lifecycle_gmredi_native_gate_core2.sh` (NEW)/
`run_lifecycle_gmredi_native_gate_multirank.sh` (NEW). Oracle (UNCOMMITTED working-tree, by design):
`fesom_pressure_dump.F90`/`fesom_geom_dump.F90` + `libfesom.so` rebuilt.
**NEXT: M5 — KPP vertical mixing + production multi-year (paper-parity).** Plan:
`docs/plans/2026-06-23-m4-gm-redi.md` (M4 log) + `2026-06-18-fesom3-architecture.md` (milestone ordering). The PRIOR
resume target follows.

**↳ (prior) M3f — whole ice step + forced lifecycle (NATIVE fluxes) + multi-rank, tag `m3`.**
M3a + M3b + M3c + M3d + M3e are ALL DONE — the ice foundation, the EVP dynamics (ocean2ice + standard EVP + mEVP), the
ice FCT advection (`ice_TG_rhs` + `ice_fct_solve`), the ice thermodynamics (`cut_off` + `thermodynamics` →
`flx_h`/`flx_fw`/`t_skin`), AND the air-sea coupling-out (`oce_fluxes_mom` + `oce_fluxes` →
`heat_flux`/`water_flux`/`virtual_salt`/`relax_salt`/`stress_surf`) all byte-match FESOM2 `max|Δ|=0`. The M3 oracle
machinery is PROVEN + EXTENDED (`fesom_ice_dump` shim now has `ice_dump_write` [M3a], `evp_dump_write` [M3b],
`ice_fct_dump_write` [M3c], `ice_thermo_dump_write` [M3d], AND `ice_flux_dump_write` [M3e]; runners
`run_ice*`/`run_evp*`/`run_icefct*`/`run_icethermo*`/`run_iceflux*`). M2 stays closed (tag `m2-mvp`); no-regression
baseline holds (13/13 ctest + step-65 1-rank/multirank + M3d-thermo-6×2 + M3e-flux-5×2 all `max|Δ|=0`). **M3f assembles
the native ice step into the lifecycle and drops the prescribed fluxes:** build a real `ice_timestep` (ocean2ice →
EVPdynamics_solve → ice_TG_rhs → ice_fct_solve → cut_off → thermodynamics) + `oce_fluxes_mom` + `oce_fluxes`, wire them
into a new forced lifecycle driver (the analog of `fesom_lifecycle` but with NATIVE fluxes replacing the
`FESOM3_FLUX_FILE` prescription), and byte-gate the multi-step forced CORE2 lifecycle vs `run_lifecycle_forced_*`
(195+ records). Then multi-rank: thread the optional `partit` + add the ice halo exchanges (the momentum-stress /
surface-flux halo reads in `oce_fluxes_mom`/`oce_fluxes` — `a_ice`/`uice`/`vice`/`srfoce`/`stress_atmoce` at the halo;
`oce_fluxes_mom` writes owned elements only so its node loop needs valid halo `stress_node_surf`; the `integrate_nod`/
`ocean_area` `allreduce_sum` is already wired guarded by `is_multirank` — and `mesh%ocean_area` needs the cross-rank
`allreduce_sum` added in `mod_mesh_areas.F90`, currently local-only, see its TODO comment). ⚠️ The real atmospheric
forcing read (NCAR bulk → `stress_atmoce_x/y` + `Ch/Ce_atm_oce_arr` + prec/runoff) must feed the native thermo/fluxes
— the M2.10 forcing read is proven, but the bulk-formula `stress_atmoce` assembly (`gen_bulk`/`fesom_forcing_dump.F90`
path) is NOT yet ported (M3f scope). ⚠️ **Any multi-rank run on levante needs the env.sh MPI flag** (sourcing `env.sh`).

Regression sanity (re-run any time — e.g. after a rebuild; all `max|Δ|=0`/green as of 2026-06-24):
```bash
./configure.sh --compiler intel --precision dp --build           # incremental (edits-only M4b–M4e);
                                                                 # --clean only when NEW files were added (L19)
bash tools/run_step_gate.sh                      # 65 records, worst |Δ|=0 (1-rank)
bash tools/run_step_gate_multirank.sh 2          # 65 records, worst |Δ|=0
bash tools/run_pressure_gate.sh                  # pi 1-rank, 95 fields, max|Δ|=0 (incl. M5a-2/3/4 kpp_* + fer_tapfac fix)
bash tools/run_pressure_gate_core2.sh            # CORE2 95 fields incl. M4a/b fer_* + M5a-2 ri_* + M5a-3/4 kpp_* (bldepth+blmix+enhance+viscAE), max|Δ|=0
bash tools/run_lifecycle_gate_core2.sh           # GM/Redi-OFF lifecycle 195 records (no-regression baseline)
bash tools/run_lifecycle_kpp_gate_core2.sh 3     # M5b: KPP-alone unforced lifecycle 195, max|Δ|=0 (CORE2 1-rank)
FER_GM=1 REDI=1 bash tools/run_lifecycle_kpp_gate_core2.sh 3  # M5b: KPP+GM+Redi 195, max|Δ|=0 (production combo)
bash tools/run_lifecycle_gm_gate_core2.sh        # M4c: GM bolus lifecycle 195 records, max|Δ|=0 (CORE2 1-rank)
bash tools/run_lifecycle_redi_gate_core2.sh      # M4d: GM+Redi lifecycle 195 records, max|Δ|=0 (CORE2 1-rank)
bash tools/run_lifecycle_fullynative_gate_core2.sh        # M3f: GM/Redi-OFF fully-native lifecycle 195 (no-regr)
bash tools/run_lifecycle_gmredi_native_gate_core2.sh 3 0  # M4e: GM+Redi fully-native lifecycle 195, max|Δ|=0 (whichEVP 0)
bash tools/run_lifecycle_gmredi_native_gate_core2.sh 3 1  # M4e: ... whichEVP 1 (mEVP)
bash tools/run_lifecycle_gmredi_native_gate_multirank.sh 2 3 0  # M4f: GM+Redi MR fully-native dist_2 195, max|Δ|=0
bash tools/run_lifecycle_gmredi_native_gate_multirank.sh 8 3 1  # M4f: ... dist_8 whichEVP 1 (needs env.sh KNEM flag, L35)
bash tools/run_icefct_gate_core2.sh              # M3c: 6 fields × whichEVP 0/1, max|Δ|=0 (CORE2 1-rank)
bash tools/run_iceflux_gate_core2.sh             # M3e: 5 fields × whichEVP 0/1, max|Δ|=0 (CORE2 1-rank)
cd build_intel_dp && ctest                       # 13/13
```
**Before assuming any kernel needs a halo extension, run the cheap `dist_N`-file invariant check first (L33).** The
optional-`partit`-absent path is byte-for-byte the proven 1-rank code, so 1-rank callers never regress.

✅ **RESOLVED — the free-surface CG "reproducibility floor" was a bug, now fixed (2026-06-21; LESSONS L29).** The L28
"iterative-solver floor" diagnosis was WRONG (a CG with byte-identical `A`/`b`/`x0`/kernels cannot manufacture
divergence from iteration count). Root cause: the preconditioner off-diagonal divide in
`oce_ssh_solve.F90::ssh_solve_preconditioner` was **auto-vectorised** (`divpd`) by the compiler while the FESOM2 oracle
compiles it SCALAR (`divsd`) — packed vs scalar division differ ~1 ULP under `-no-prec-div -fimf-use-svml`, so the
**un-gated `pr_values` array** drifted in 1299/870146 entries, seeding `z=M⁻¹r` and surfacing in the CG residual at
~iter 6 (CORE2 136 iters; pi's 37 iters stayed below the bit). The asymmetry: the oracle writes via a LOCAL POINTER
(compiler can't disprove aliasing → scalar), F3 wrote the DERIVED-TYPE COMPONENT (provably distinct → vectorised).
**Fix: one `!DIR$ NOVECTOR`** on the precond loop → scalar `divsd`, `pr_values` `max|Δ|=0`. **Verified (production,
un-instrumented binaries):** `tools/run_pressure_gate_core2.sh` PASS `max|Δ|=0` on all 28 fields incl.
`d_eta`/`eta_n`/`uv_upd`/`hbar`; and the previously-blocked `tools/run_lifecycle_gate_core2.sh` now MATCHes — **195
records (13 substeps × 5 probes × 3 steps), worst |Δ| = 0**. The whole dynamical core is byte-exact across multiple
steps on CORE2. Method + generalisable lessons (gate every consumed intermediate; a divide byte-matches only when
operands AND SIMD width match) in LESSONS L29. *Historical scoping notes + the DONE sub-gates follow.*

M2.11 was decomposed into byte-gated sub-tasks (mirroring M2.10a/b/c); scoping done **2026-06-20** (3 oracle-investigation
passes — all input paths verified: CORE2 mesh `/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2/`, IC `phc3.0_winter.nc`):

- **M2.11a — CORE2 geometry byte-gate — ✅ DONE (2026-06-20; `max|Δ|=0`, 19 fields, `tools/run_geom_gate_core2.sh`).**
  Scaled mesh-read ~40× (CORE2 `nod2D=126858`/`elem2D=244659`/`edge2D=371644`/`edge2D_in=362333`/`nl=48`) and **CLOSED
  the L8 CW-swap deferral**: FESOM3 made **244654/244659** `enforce_cw_orientation` swaps (≈100%, mesh stored CCW; pi
  had 0) and the post-swap `elem2D_nodes`/centroid/`elem_area`/`gradient_sca`/`area` are all `max|Δ|=0` vs the oracle's
  runtime `test_tri`. **`tools/make_dist1.py`** hand-crafts dist_1 (validated vs the pi template); `fesom_geomdump`
  needed NO code change (just `FESOM3_MESH_DIR`); the geom shim stops at mesh_setup (~2 s, mesh-only). See LESSONS L26.
  ↳ **ACTIVE NEXT: M2.11b.** *(historical scoping notes follow.)* **No `dist_1` existed** (only dist_2/4/…/512) →
  hand-crafted from the pi template
  (`tests/data/MESHES/pi/dist_1/`: `rpart.out`=npes/nod2D/identity-list, `my_list00000.out`=identity,
  `com_info00000.out`=empty halo w/ blank lines). Run the existing FESOM2 geom-dump shim (`fesom_geom_dump.F90`, fires
  at `mesh_setup`, STOPS before forcing → NO lifecycle/output-crash exposure) + `fesom_geomdump` on CORE2 1-rank;
  extend `tools/run_geom_gate.sh` for CORE2 → `max|Δ|=0`. **Safe deferrals confirmed:** min `nlevels=5` (nodes+elems)
  → NO single-layer columns (L18 TDMA OOB guard NOT needed); **cavity + partial-cell OFF**.
- **M2.11b — initial conditions (`do_ic3d`) — ✅ DONE (2026-06-21; `max|Δ|=0`, 3 fields [Z_3d_n/ic_temp/ic_salt], CORE2
  1-rank, `tools/run_ic_gate_core2.sh`).** Ported `oce_initial_state`→`do_ic3d` (`gen_ic3d.F90:493`) into
  `src/oce/oce_initial_state.F90` = netCDF 3D-double read (`nc_get_var3d_dp` in `mod_io_netcdf` + reused
  `forcing_binarysearch`/bilinear from `mod_forcing_read`) + vertical LINEAR interp onto `Z_3d_n` + `extrap_nod3D`
  (`gen_support.F90:400`; fill `dummy=1e10` by Gauss-Seidel neighbour-average — the partition-order step, deterministic at
  1-rank) + Kelvin guard + `insitu2pot` (`ptheta` Bryden-1973 RK4 + `atg`, in `oce_pressure_bv.F90`). **IC file CONFIRMED:
  `phc3.0_winter.nc`** (360×180×33, the work_core production file — NOT woa18; the test harness path also has it). phc3.0 land
  = NaN (no `_FillValue`) → missing-value mask = `ieee_is_nan`. `idlist=2,1` ⇒ salt→data(2) FIRST, temp→data(1) SECOND,
  `t_insitu=.true.`. Oracle shim `fesom_ic_dump.F90` dumps the live `Tclim`/`Sclim` (no prescribe). See LESSONS L27.
- **M2.11c-1 — unforced multi-step lifecycle + CORE2 dynamics-kernel gate — ✅ DONE (2026-06-21; whole dynamical core
  incl. CG `d_eta` `max|Δ|=0` over multiple steps post-L29).** Built: FESOM3 `src/drivers/fesom_lifecycle.F90` (cold-start CORE2 + `do_ic3d`
  IC + N-step `step_oce` loop, multi-step AB2 in place), oracle `tools/run_lifecycle_core2.sh` (REAL multi-step
  `oce_timestep_ale` + built-in dump_shim, `use_ice=.false.` unforced — `forcing_setup` is then a no-op, all surface
  fluxes 0), `tools/run_lifecycle_gate_core2.sh` (compares 13 NODE substeps × 5 probes × N steps), AND the NEW
  `tools/run_pressuredump_core2.sh`+`tools/run_pressure_gate_core2.sh` (the M2.1-M2.9 per-kernel gate on CORE2 —
  `fesom_pressuredump` now env-configurable `FESOM3_STEP_PER_DAY=48`). Oracle fixes (uncommitted): `io_meandata.F90::output`
  + `io_restart.F90::write_initial_conditions` both `if(partit%npes==1) return` (the 1-rank `io_gather init_nod2D_lists`
  bug hits BOTH output and restart-write); `use_sw_pene=.false.` in the run dir (else `sw_3d` unallocated → segfault, L28).
  **After the L29 fix the CORE2 pressure gate is `max|Δ|=0` on all 28 dynamical-core fields INCLUDING `d_eta`, and the
  multi-step `run_lifecycle_gate_core2.sh` MATCHes (195 records, worst |Δ|=0).** Debug `-check all` clean. See LESSONS L29.
- **M2.11c-2 — FORCED lifecycle — ✅ DONE (2026-06-21; forced dynamical core byte-exact on CORE2, CG floor resolved by L29, re-verify CONFIRMED).**
  The REAL forced FESOM2 lifecycle (`use_ice=.true.`, real CORE2 NCAR forcing at the 1948 stubs + pool runoff/SSS, the ice
  EVP + `oce_fluxes` producing the air-sea fluxes) runs cleanly at 1-rank on CORE2 — **the ice-at-1-rank unknown is
  DE-RISKED** (`tools/run_lifecycle_forced_core2.sh`). Built: oracle `port2/fesom2/src/fesom_flux_dump.F90` (per-step
  full-field dump of `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt` [nod2D] + `stress_surf` [2,elem2D] BEFORE
  `oce_timestep_ale`, env `FESOM_FLUX_DUMP`, wired in `fesom_module.F90` runloop), FESOM3 `fesom_lifecycle` reads them per
  step (`FESOM3_FLUX_FILE`) and prescribes into `step_oce` (the M2.5/M2.8 prescribe-the-unsourced-input pattern), gate
  `tools/run_lifecycle_forced_gate_core2.sh`. **Result identical to unforced — now byte-exact on ALL substeps:** the
  re-verify ran 2026-06-21 POST-L29 and `tools/run_lifecycle_forced_gate_core2.sh` MATCHes — **195 records (13 substeps ×
  5 probes × 3 steps), worst |Δ| = 0** (the pre-L29 `d_eta` drift 5.5e-17 → 3.5e-6 is GONE; the NOVECTOR fix makes the
  forced CG byte-exact exactly like the unforced). FESOM3 step diagnostics match the oracle to all printed digits
  (step 1 uv=0.23283/eta=0.34914 … step 3 uv=0.45194/eta=0.67866). So **the M3-gap flux prescription is byte-exact.**
  Calendar: CORE (noleap) forcing REQUIRES `include_fleapyear=.false.`; `use_sw_pene=.false.` (matches the ported step_oce
  — no `sw_3d` term). Debug `-check all` clean (incl. the flux-read path).
- **SCOPE (LESSONS L25):** `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt` air-sea budget = **M3** (the `obudget`
  in `ice_thermo_oce.F90` + `oce_fluxes`) → prescribed from the oracle dump at M2.11c. The SSS/runoff/chl climatology
  reads + the JRA55 (gregorian) forcing variant are deferred (CORE2 NCAR stubs suffice). M1's multi-rank advection gate
  + the local-mesh remap ride **M2.12**.

- **The reduced-M2 oracle namelist was NEVER needed** (through M2.9b). Every M2.x gate runs the prescribe-and-stop
  shim at end-of-`ocean_setup`, which drives the REAL kernels (up to the FULL `oce_timestep_ale` at M2.9b) on
  PRESCRIBED state — `oce_timestep_ale` reads forcing ARRAYS (set by the shim), not files, so the 1-rank forcing-file
  hang (L8) is bypassed. The shim FORCES the reduced-M2 dispatch (`mix_scheme_nmb=2`/`Fer_GM=.false.`/`Redi=.false.`/
  `opt_visc=7`/`use_wsplit=.false.`) over the shipped pi namelist (which keeps all GM/Redi/KPP arrays allocated).
  M2.10 forcing is where reading real files + running PAST the shim's stop first becomes necessary.
- **`use_wsplit=.true.` / `adv_tra_vert_impl` is an unported kernel** (do_oce_adv_tra's FCT implicit-vertical-advection
  correction; guarded with `error stop`). M2.9b forced `use_wsplit=.false.` (M1.4 precedent). Port it + integrate
  the M2.7 explicit/implicit w-split into a real multi-step run in a later dynamics pass.

**M1 multi-rank gate (folded into M2.12b):** the local-mesh remap (global→local
numbering/connectivity/geometry, com-structs) is ✅ DONE (M2.12a — geometry byte-gated `max|Δ|=0` on dist_2+dist_8;
the multi-rank `enforce_cw_orientation` swap path is CLOSED). **M2.12b** then lifts M1.1–M1.4 advection to multi-rank
+ adds the `find_neighbors` halo `nod_in_elem2D` dance, and **M2.12c** gates the WHOLE model on 1/8/32-rank. For
M2.12b re-confirm the M1.1–M1.4 kernels' dropped halo exchanges (tr_xy/edge_up_dn_grad/fct_LO/fct_plus_minus/
del_ttf) + loop bounds (myDim vs myDim+eDim) are correctly lifted; the gate target is the post-exchange OWNED values
on the SAME partition (L8 accumulation-order caveat). Still deferred to a richer mesh: the FCT
`AUX`/`edge_up_dn_grad`-scratch cavity caveat (L11).

## M2.12 entry notes (scoped 2026-06-21; M2.12a + M2.12b DONE — the rest carries to M2.12c)

The full decomposition (M2.12a/b/c) lives in `docs/plans/2026-06-18-fesom3-architecture.md` Task M2.12.
**M2.12a (local-mesh remap + geometry gate) AND M2.12b (multi-rank tracer advection gate) are ✅ DONE** (both
`max|Δ|=0` pi dist_2+dist_8, see "Next task"). **Active next = M2.12c** (whole-model: lift the dynamics + tracer
SOLVE, threading the optional `partit` exactly like M2.12b). The notes below are now mostly reference; the
M2.12b-flagged items (the `find_neighbors` halo dance, `exchange_elem_full`, the area/elem_area exchanges) are DONE.

- **Foundation EXISTS — do NOT rebuild it.** `src/infra/mod_halo.F90` is a real multi-rank MPI exchange
  (manual pack → Isend/Irecv → `MPI_Waitall` → broadcast-only owner→halo unpack; `exchange_nod` 2D/3D-real +
  2D-int, `exchange_elem` 2D/3D-real), and `src/infra/mod_partitioning.F90::read_dist_partition` already loads
  `myList_nod2D/elem2D/edge2D` (local→global) + `com_nod2D`/`com_elem2D`/`com_elem2D_full`. Both tested
  1/2/8-rank (`test_partit`, `test_halo`). FESOM2 uses precompiled MPI_TYPE_INDEXED datatypes vs FESOM3's
  manual pack — same bytes moved, no arithmetic → byte-identical. **The geometry pipeline
  (`compute_geometry`/`setup_vertical`/`enforce_cw_orientation`) is already partition-agnostic** (loops over
  `mesh%nod2D`/`elem2D`/`edge2D`, reads only local arrays) → it runs UNCHANGED once the local arrays are built.
- **Gate rule (L8) — the load-bearing constraint.** Compare FESOM3 `dist_N` vs FESOM2 `dist_N` PER-RANK on
  OWNED entries (1..myDim). NOT vs 1-rank global: both `dist_N` runs share the same `myList` order → local idx
  i ↔ same global id AND the per-node area sum order matches → byte-identical; the 1-rank global uses a
  different element permutation (non-associative FP) → would mismatch at ~ULP.
- **Algorithm (FESOM2 `port2/fesom2/src/oce_mesh.F90`):** ⟶DONE in M2.12a (`read_mesh_local`): `read_mesh` :212
  (scatter — implemented with a full-size inverse map, not FESOM2's chunked `mapping`; identical result); edges read
  :1787 (scatter + negative-localize); owned `nod_in_elem2D` build (`find_neighbors` :2021-2049, OWNED part);
  `setup_vertical` `exchange_nod(nlevels_nod2D_min)` :1669. `elem2D_nodes` is `(MAX_NV,myDim_elem2D)` — OWNED only.
  ⟶DONE in M2.12b: the `find_neighbors` HALO dance — `exchange_nod(num)` :2051, pack global elem-ids :2057
  → `exchange_nod` → re-localize through **eXDim** :2066-2073 (`complete_nod_in_elem_halo` in `read_mesh_local`,
  completes halo `nod_in_elem2D` for MUSCL/advection); and `mesh_areas` :2162 — `compute_geometry` now does
  `exchange_nod(area/areasvol/*_inv)` :2322 + `exchange_elem_full(elem_area)` :2220 at npes>1 (the owned-node area
  accumulation stays owned-only, the exchange fills the halo). `elem_neighbors`/`elem_edges` still SKIPPED (not
  needed for advection; may be needed by the M2.12c dynamics — re-check when lifting `momentum_adv`/visc).
- **Infra note (M2.12a + M2.12b):** the owned-entry GEOMETRY gate (M2.12a) did NOT need `exchange_elem_full` — the
  fix there was precompute owned element centers + `exchange_elem` them. **M2.12b ADDED `exchange_elem_full`**
  (`com_elem2D_full`, eDim+eXDim; 2D-r/2D-i/3D-r variants + `core_blk_r`) — needed for `elem_area` AND `tr_xy` AND
  the `coord_elem`/`e_nodes` of `find_up_downwind_triangles` (all read at eXDim through a halo node's element list).
- **Oracle (⟶DONE):** `port2/fesom2/src/fesom_geom_dump.F90` extended to npes>1 (per-rank LOCAL OWNED dump,
  `<path>.<mype5>`); `libfesom.so` rebuilt. **Gate (⟶DONE):** `tools/run_geom_gate_multirank.sh [np]`
  (FESOM2 dist_N dump vs FESOM3 dist_N dump, per-rank owned, `geom_diff.py`); `max|Δ|=0` on dist_2 + dist_8.
  CLOSED the deferred multi-rank `enforce_cw_orientation` swap path. Lessons in **L30**.
- **Build discipline:** clean-rebuild after adding NEW files (L19); `ulimit -s unlimited` for the dump writer (L20).

## M3 entry notes (scoped 2026-06-22)

The full M3 decomposition (M3a–M3f) is in "Next task". M3 is sea ice (standard EVP) → the air-sea budget. The
FESOM2 source is the same v2.7.3 oracle: `ice_EVP.F90` / `ice_fct.F90` / `ice_thermo_oce.F90` / `ice_oce_coupling.F90`
/ `ice_setup_step.F90` / `MOD_ICE.F90`. **M3a–M3e ✅ DONE (all `max|Δ|=0`); active next = M3f** (assemble the native
ice step + `oce_fluxes` into a forced lifecycle, drop the prescribed fluxes, then multi-rank; tag `m3`).

- **Call sequence (FESOM2 `fesom_module.F90` runloop):** `ocean2ice` (:677, set ice's view of the ocean surface) →
  `ice_timestep` (:707 = `ice_setup_step.F90:96`: `EVPdynamics`(:215) → `ice_TG_rhs`(:258) → `ice_fct_solve`(:261) →
  `cut_off`(:295) → `thermodynamics`(:314)) → `oce_fluxes_mom`(:715) → `oce_fluxes`(:716). Ice runs at `ice_dt =
  ice_ave_steps·dt` (=`dt` here; `Tevp_inv=3/ice_dt`); `evp_rheol_steps=120` subcycles; the σ stress tensor is
  PROGNOSTIC (elastic memory persists across steps — already a `t_ice` field in `mod_ice.F90`).
- **Oracle = the proven forced lifecycle.** `tools/run_lifecycle_forced_core2.sh` ALREADY runs the REAL FESOM2 ice
  EVP+advection+thermo+oce_fluxes at 1-rank on CORE2 (`use_ice=.true.`, real NCAR forcing) and dumps the air-sea fluxes
  via `fesom_flux_dump.F90`. M3 per-kernel gates extend this: add a prescribe-and-stop `fesom_ice_dump` shim (mirror
  `fesom_step_dump`) that runs `ice_setup` + one `ice_timestep` + `oce_fluxes` on the do_ic3d-IC'd ocean and dumps each
  intermediate (`mod_dump` byte format), per-kernel gated like M2.1–M2.9. The flux dump is the M3e/M3f target.
- **Reuses (don't rebuild):** `ice_mass_matrix_fill` rides the **M2.6 `ssh_stiff` CSR** (rowptr/nn_num/nn_pos); the
  ice FCT advection reuses the M1/M2.12b tracer-FCT pattern; the cold-start IC reads `SST` from the M2.11b `do_ic3d`
  state + `geo_coord_nod2D`; dumps reuse `mod_dump`/`mod_binary_arrays`. `t_ice` skeleton already in
  `src/types/mod_ice.F90` (sigma11/12/22 prognostic, a_ice/m_ice/m_snow, uice/vice, work) — flesh it out from `MOD_ICE.F90`.
- **Config (reduced, like reduced-M2):** standard EVP only — `whichEVP=0`; NO `__icepack`/`use_meltponds`/`use_cavity`/
  `__oasis`/`__yac`/`__oifs`; `num_itracers=3` (a_ice/m_ice/m_snow). Namelists `&ice_dyn` + `&ice_therm` (defaults in
  `MOD_ICE.F90`: `pstar=30000`, `ellipse=2`, `c_pressure=20`, `delta_min=1e-11`, `evp_rheol_steps=120`,
  `ice_gamma_fct=0.25`, `ice_diff=10`, `cd_oce_ice=5.5e-3`, `Sice=4`, `h0=0.5`, densities rhoice=910/rhosno=290/...).
- **Gate rule (same L8):** per-rank owned dumps vs same-partition FESOM2; ice fields are nod2D (a/m/uv/flx) or elem2D
  (sigma). Watch (plan §M3): `bc_index_nod2D` multi-rank safety, NaN-vs-0 ice masking in diagnostics, ice_dt sync.

## Open notes / risks

- The FESOM2 oracle is built + PROVEN (M0.7 geometry + M1.1–M1.4 advection + M2.1–M2.9a
  pressure/EOS/N²/PGF/vel_rhs/momadv/viscosity/ivertvisc/SSH/ALE-update/PP-mixing/convective-adjustment/tracer-solve gates
  all `max|Δ|=0`, **57 fields**, `tools/run_pressure_gate.sh`) + M2.9b the WHOLE assembled `oce_timestep_ale`
  (`tools/run_step_gate.sh`, **65 records** vs the REAL step's built-in dumps, `--ignore-substep=2`).
  Self-tests run standalone (13/13 ctest). M2+ kernel gates extend the proven 1-rank end-of-`ocean_setup` shim
  pattern (`fesom_pressure_dump.F90` drives the isolated kernels; `fesom_step_dump.F90` drives the REAL
  `compute_vel_nodes` + `oce_timestep_ale` on prescribed state). The prescribe-and-stop shim STILL sufficed for
  M2.9b — the assembled step reads forcing ARRAYS (set by the shim), not files, so the reduced-M2 namelist + the
  forcing-file read were never needed. **M2.10 forcing is next** (the FIRST netCDF I/O), see "Next task".
- **⚠️ CLEAN-rebuild rule (L19):** when a milestone adds NEW `src/**/*.F90` files (not just edits), run
  `./configure.sh --compiler intel --precision dp --clean --build` before trusting the byte-gate — an
  incremental build re-runs CMake configure (GLOB `CONFIGURE_DEPENDS`) and can link objects against a MIX
  of stale/fresh `.mod` interfaces, ULP-drifting EVERY field (the M2.6 first-gate red herring).
- **M2.12 (multi-rank, MVP exit) is ✅ DONE** (a/b/c-1/c-2/c-3). The local-mesh remap + the whole-model multi-rank
  byte-match (all 65 substeps, dist_2+dist_8) landed cleanly. The deferred cavity caveat (L11) + dist_32 (pi has no
  such partition) ride a richer mesh; M3 is next.
- **M2.12c-1/c-2/c-3 uncommitted artifacts:** FESOM3 — drivers `src/drivers/fesom_stepdump_mr.F90` (c-1/c-2: dynamics
  chain → d_eta, direct kernel calls) + **`src/drivers/fesom_stepfull_mr.F90` (c-3 NEW: the WHOLE step via `step_oce`
  through `partit`)**; gates `tools/run_stepdyn_gate_multirank.sh` (c-1/c-2) + **`tools/run_step_gate_multirank.sh`
  (c-3 NEW: all 65 substeps, ignore only SW_AB id=2)**; `mod_halo.F90` (c-1 rank-3 `exchange_nod`; c-2 `allreduce_sum`),
  `mod_part_bounds.F90` (`local_dims`); the c-1 lifted dynamics kernels (`oce_pressure_bv`/`oce_pgf`/`oce_ale_mixing_pp`/
  `oce_mo_conv`/`oce_dyn_velrhs`/`oce_dyn_visc`/`oce_dyn_ivertvisc`), c-2 `oce_ssh_rhs.F90`+`oce_ssh_solve.F90`, and
  **c-3 `oce_ale.F90` (update_vel/compute_hbar_ale/update_eta_n/vert_vel_ale/compute_CFLz/compute_Wvel_split/
  update_thickness_ale optional partit) + `oce_ale_tracer.F90` (solve_tracers_ale/diff_tracers_ale/diff_part_hor_redi/
  diff_ver_part_impl_ale optional partit) + `mod_step_oce.F90` (optional partit threaded to every kernel)**. FESOM2
  oracle — **UNCHANGED for c-3** (the c-1 `fesom_step_dump.F90` npes>1 extension — dropped the `npes/=1` return; nodes
  owned+halo, element velocity OWNED + `exchange_elem(UV)`; ⚠️ do NOT prescribe UV at eXDim elements,
  `elem2D_nodes(1)` can be beyond the eDim node halo → `coord_nod2D` OOB — already drives the whole `oce_timestep_ale`
  and dumps ALL substeps). Rebuild the oracle with `make -C port2/fesom2/build fesom.x` only if its tree changed.
- The FESOM2 oracle shim edits (`fesom_advhor_dump.F90`, `fesom_pressure_dump.F90`, `fesom_step_dump.F90` [now npes>1],
  `fesom_forcing_dump.F90`, `fesom_ic_dump.F90` [M2.11b], `oce_setup_step.F90`, `fesom_module.F90`, geom/ale shims) **+
  the M2.10a `io_netcdf_workaround_module.F90` np=1 fix** live as UNCOMMITTED working-tree instrumentation in
  `port2/fesom2` (not committed there, by design); `libfesom.so` must be rebuilt (`make -C port2/fesom2/build
  fesom.x`, then re-run `cmake .` first if a NEW shim file was added so the GLOB picks it up) after editing a
  shim. The `next_io_rank` np=1 fix is behavior-preserving (sequential I/O, value-identical) and arguably a real
  bug fix (the code's own TODO admits the recursion never ends at 1 rank).
