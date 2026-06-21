# FESOM3 — Handoff (durable state across sessions)

Single source of truth for "where are we / what's next". Update at the end of every task.
Full pre-M2.12 milestone detail + the per-gate recipes live in [`HANDOFF-archive.md`](HANDOFF-archive.md).

## Where we are

- **Milestone:** M2 (minimal ocean dynamical core) — **IN PROGRESS.**
  **CURRENT STATUS (2026-06-21):** the whole M2 dynamical core is BYTE-COMPLETE — `max|Δ|=0` vs FESOM2 on pi 1-rank
  (per-kernel + whole-step) AND across multi-step CORE2 (unforced + forced lifecycle, incl. the free-surface CG).
  Now in **M2.12 (multi-rank, MVP exit): M2.12a — the local-mesh remap + per-rank geometry byte-gate — is ✅ DONE
  (`max|Δ|=0` on pi dist_2 + dist_8); ACTIVE NEXT = M2.12b** (multi-rank advection). See "Next task".

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
  fields, every rank, on pi **dist_2 (2 ranks) AND dist_8 (8 ranks)** vs same-partition FESOM2
  (`tools/run_geom_gate_multirank.sh`). Built `mod_mesh_read.read_mesh_local` (global→local scatter via full-size
  inverse maps + owned `nod_in_elem2D` + `enforce_cw`(owned) — CLOSES the multi-rank CW-swap caveat) and made
  `mod_mesh_areas.compute_geometry` partition-aware (`local_bounds`; owned element CENTERS precomputed +
  `exchange_elem`'d so owned-edge `edge_cross_dxdy` resolves halo-element centers; owned-node areas computed
  locally). Per-rank owned dump in `mod_geom_dump` + the oracle shim `fesom_geom_dump.F90` (npes>1). **No
  regression:** 1-rank geom gate (19), 1-rank pressure gate (57), 13/13 ctest all still `max|Δ|=0`/green.
  **DEFERRED to M2.12b** (not needed for the owned-entry geometry gate; needed when dynamics consume halos): the
  `find_neighbors` `nod_in_elem2D` halo dance (eXDim re-localize) + the area/elem_area halo exchanges.
- **ACTIVE NEXT: M2.12b** (lift M1.1–M1.4 advection to multi-rank: exchange_nod/elem of the intermediates + the
  `find_neighbors` dance; gate `del_ttf` on pi dist_2/dist_8). M3 (ice + `oce_fluxes` air-sea budget) remains the
  other unblocked milestone if priorities change.

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

## M2.12 entry notes (scoped 2026-06-21; M2.12a DONE — the rest carries to M2.12b/c)

The full decomposition (M2.12a/b/c) lives in `docs/plans/2026-06-18-fesom3-architecture.md` Task M2.12.
**M2.12a (local-mesh remap + per-rank geometry gate) is ✅ DONE** (`max|Δ|=0` pi dist_2+dist_8, see "Next task").
**Active next = M2.12b** (multi-rank advection gate). The notes below still hold; what M2.12a consumed vs what
carries to M2.12b is flagged inline (⟶done / ⟶M2.12b).

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
  ⟶M2.12b (still to port): the `find_neighbors` HALO dance — `exchange_nod(num)` :2051, pack global elem-ids :2057
  → `exchange_nod` → re-localize through **eXDim** :2066-2073 (completes halo `nod_in_elem2D` for MUSCL/advection);
  and `mesh_areas` :2162 area over owned+halo nodes :2258 with `exchange_nod(area)` :2322 + `exchange_elem(elem_area)`
  :2220 (fills the HALO area that the dynamics read).
- **Infra note (corrected by M2.12a):** the owned-entry GEOMETRY gate did NOT need `exchange_elem_full` — owned-edge
  `edge_cross_dxdy` only reads an eDim neighbour's CENTER, covered by the existing `exchange_elem` (com_elem2D); the
  fix was precompute owned element centers + `exchange_elem` them (centers are WP; route the MP `elem_cos` via a WP
  scratch). `exchange_elem_full` (com_elem2D_full, eDim+eXDim) may still be needed in M2.12b for `elem_area` over
  eXDim. SKIP `elem_neighbors`/`elem_edges` (not needed for geometry; M2.12b may need `elem_neighbors`).
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
- **M2.12 is decomposed (M2.12a/b/c) and underway.** M2.12a (local-mesh remap + geometry gate) is DONE; the
  remaining weight is M2.12b (folded M1 advection multi-rank gate + the `find_neighbors` halo dance) + M2.12c
  (whole-model 1/8/32-rank byte-match + the deferred cavity caveat). The remap landed cleanly without splitting.
- The FESOM2 oracle shim edits (`fesom_advhor_dump.F90`, `fesom_pressure_dump.F90`, `fesom_step_dump.F90`,
  `fesom_forcing_dump.F90`, `fesom_ic_dump.F90` [M2.11b], `oce_setup_step.F90`, `fesom_module.F90`, geom/ale shims) **+
  the M2.10a `io_netcdf_workaround_module.F90` np=1 fix** live as UNCOMMITTED working-tree instrumentation in
  `port2/fesom2` (not committed there, by design); `libfesom.so` must be rebuilt (`make -C port2/fesom2/build
  fesom.x`, then re-run `cmake .` first if a NEW shim file was added so the GLOB picks it up) after editing a
  shim. The `next_io_rank` np=1 fix is behavior-preserving (sequential I/O, value-identical) and arguably a real
  bug fix (the code's own TODO admits the recursion never ends at 1 rank).
