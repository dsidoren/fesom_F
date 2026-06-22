# FESOM3 — Handoff (durable state across sessions)

Single source of truth for "where are we / what's next". Update at the end of every task.
Full pre-M2.12 milestone detail + the per-gate recipes live in [`HANDOFF-archive.md`](HANDOFF-archive.md).

## Where we are

- **Milestone:** M2 (minimal ocean dynamical core) — **multi-rank MVP byte-match ACHIEVED.**
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
- **ACTIVE NEXT: M3 — sea ice (EVP).** The prescribed `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt` become the
  real air-sea budget (sea ice EVP + thermo → `oce_fluxes` → the obudget in `ice_thermo_oce.F90`). Scope sketch in the
  Roadmap (plan §M3).

**↳ RESUME HERE (next session): M3 (sea ice EVP).** The multi-rank free-running lifecycle is byte-exact; the m2-mvp
single-step gates stay `max|Δ|=0`/green (no-regression baseline). The production-validation work (UV eDim fix +
`do_ic3d` MR lift + `fesom_lifecycle_mr` + the `env.sh` `single_copy_mechanism=none` MPI fix) is **committed to
`main`** (the M2.12 production-validation commit, immediately after `m2-mvp`), F3-side §6 diagnostics stripped. (The FESOM2 oracle's inert diagnostic shims remain in its own repo —
out of scope here, harmless.) ⚠️ **Any multi-rank run on levante needs the env.sh MPI flag** (sourcing `env.sh` sets
it) — without it, large vader messages silently corrupt.

Regression sanity (re-run any time — e.g. after a rebuild; all `max|Δ|=0`/green as of 2026-06-22):
```bash
./configure.sh --compiler intel --precision dp --clean --build   # clean rebuild (L19) if NEW files were added
bash tools/run_step_gate_multirank.sh 2          # 65 records, worst |Δ|=0
bash tools/run_step_gate_multirank.sh 8          # 65 records, worst |Δ|=0
bash tools/run_step_gate.sh                      # 65 records, worst |Δ|=0 (1-rank)
bash tools/run_stepdyn_gate_multirank.sh 2       # 30 records (c-1/c-2), worst |Δ|=0
bash tools/run_stepdyn_gate_multirank.sh 8       # 30 records, worst |Δ|=0
bash tools/run_advhor_gate_multirank.sh 8        # all ranks max|Δ|=0
bash tools/run_pressure_gate.sh                  # 57 fields, max|Δ|=0 (1-rank)
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
