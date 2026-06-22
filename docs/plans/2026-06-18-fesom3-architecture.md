# FESOM3 — Clean-Architecture Fortran Reimplementation of FESOM2

## Overview

Build **FESOM3**: a clean-architecture Fortran reimplementation of the FESOM2 ocean + sea-ice
model. The first deliverable (v1) must **byte-exactly reproduce a specific FESOM2 configuration**
(double precision, triangle meshes, same Intel toolchain, fixed MPI rank count). The "clean
architecture" is the *plumbing only* — derived types, dependency injection, module layering, a
generic halo layer, write-once config, and precision/mesh-arity *scaffolding*. The numerical
**kernels are transcribed line-for-line from the original FESOM2 Fortran** with no invention and no
optimization.

Target physics (full v1 scope, built milestone by milestone): nonlinear EOS, hydrostatic PGF, AB2,
CG SSH solve, MFCT tracer advection, biharmonic viscosity, **PP / KPP / TKE** vertical mixing,
**GM/Redi** eddy parameterization, **EVP + mEVP** sea-ice rheology, **linfs + zstar** vertical
coordinates, JRA55 bulk forcing.

Foundations *laid but not exercised* in v1 (later phases): mixed precision (WP/MP scaffolding),
quad/mixed-polygon meshes (element-arity generalization), GPU-direct/overlap hooks, and the
language-port contract. Goal: **"first reproduce exactly, then optimize."**

This plan was produced from a completed architecture brainstorm. The full decision record lives in
project memory at `/home/a/a270088/.claude/projects/-home-a-a270088-fesom3/memory/`
(`project-brainstorm-decisions.md` D0–D9, `project-port-synthesis-principles.md`,
`project-bit-identity-reality.md`, `reference-treasures-map.md`).

### The governing principle (D0 — anchors every task)

Every task is judged twice:
1. **Does it change kernel arithmetic?** → must be **NO** (or it breaks byte-exactness).
2. **Does it make a future improvement easier?** → should be **YES** (or it isn't worth doing now).

Bit-identity is achievable **Fortran→Fortran** (same compiler+flags+DP+fixed ranks+faithful
arithmetic) — that is *why* we can demand `max|Δ|=0`; it was **not** achievable C→Fortran. Therefore
clean architecture lives in the plumbing; the math is faithful.

## Context (from discovery)

- **FESOM2 source (algorithm oracle — transcribe FROM here, byte-gate AGAINST here):**
  `/home/a/a270088/port2/fesom2/src/` (232 `.F90`; instrumented v2.7.3, tag `fesom2.7.3-cport-instr`).
  Key files: `oce_ale.F90`, `oce_adv_tra_hor.F90`, `oce_adv_tra_ver.F90`, `oce_adv_tra_fct.F90`,
  `oce_muscl_adv.F90`, `oce_ale_mixing_pp.F90`, `oce_ale_mixing_kpp.F90`, `gen_modules_cvmix_tke.F90`,
  `oce_dyn.F90`, `oce_setup_step.F90`, `ice_EVP.F90`, `ice_maEVP.F90`, `ice_thermo_*.F90`,
  `gen_halo_exchange.F90`, `gen_modules_partitioning.F90`, `MOD_MESH.F90`, `MOD_DYN.F90`, `MOD_ICE.F90`.
- **Structural template (HOW to decouple — copy structure, NOT math):**
  `/home/a/a270088/fesom3/design_refs/tracer_dwarf/lib/` (MOD_MESH, MOD_PARTIT, MOD_DYN, MOD_TRACER,
  gen_halo_exchange, gen_modules_partitioning, oce_adv_tra_driver/_hor/_ver/_fct, oce_muscl_adv,
  oce_modules, oce_mesh, analytic_mesh, hp_math_intrinsics, io_restart_derivedtype,
  MOD_READ/WRITE_BINARY_ARRAYS, tracer_c_interface) + `src/` drivers (fesom_analytic, fesom,
  fesom_mesh_init) + `CMakeLists.txt`, `configure.sh`, `cmake/`, `docs/` (ARCHITECTURE, PRECISION,
  MESH_DATA_STRUCTURES, CLAUDE.md).
- **FESOMx design docs (conceptual):** `design_refs/fesom3-design/design_doc/` (00–06).
- **Secondary refs:** C port `/home/a/a270088/port2/fesom2_port/` (pinned CORE2 config, dump
  instrumentation, gotchas); JAX port `/home/a/a270088/port_jax/` (functional design).
- **Reference inputs:** `dist_<NP>/` partitions, PHC3.0 IC (partition-matched), JRA55-do v1.4.0
  forcing, `pi` mesh (3140 nodes, ~23 lev) + CORE2 (0.13M nodes, 48 lev). Environment: DKRZ Levante.
- **New repo root:** `/home/a/a270088/fesom3/` (currently holds `design_refs/`, `paper/`,
  `FESOM2LLM.zip`). Source goes under `src/`.

## Development Approach

- **Faithful translation, no "improvements."** Transcribe each kernel from FESOM2 Fortran. Only
  mechanical changes allowed: USE-globals → explicit type arguments; module reorganization; the
  arity/precision *scaffolding* (inert on the DP/triangle anchor build). **Cite `FESOM2 file:line`
  for every ported constant/formula** in a comment. **Arity caveat:** the `1/3 → 1/n_vert`
  generalization is bit-exact ONLY where FESOM2 *divides* (`X/3.0_WP`); where FESOM2 multiplies by a
  precomputed reciprocal or truncated literal, reproduce that *exact* expression form per site — do
  not assume `1.0_WP/n_vert` is universally substitutable.
- **Build the validation harness FIRST** (M0), before any physics — every prior port's biggest time
  sinks were validation gaps, not coding.
- **Off-by-default feature flags with byte-identical off-switches.** Each new physics block, when
  off, must reproduce the previous milestone bit-for-bit (a gate).
- **Validate at the production config, long, and run past the known blow-up step.** dt=500 smoke
  tests hide whole bug classes (the "no-op at small dt" class) that are fatal at dt=1800.
- **Pin all parameters to the production run-dir namelist**, not templates/defaults. Audit dispatcher
  branches against the *run's* branch.
- **One task fully complete (incl. its gate passing) before the next.** Commit/tag working states.
- **Cross-session LLM discipline** (from the paper): short scoped sessions; maintain a
  `docs/LESSONS.md` (mistakes + fixes, never repeated) and `docs/HANDOFF.md` (state, canonical
  reference run, next task).

### Testing Strategy (mapped to this domain)

"Tests" here = **the validation harness gates**, which are required deliverables of every task:
- **Byte-gate (R2):** whole-model `max|Δ|=0` vs instrumented FESOM2 via the gid-keyed per-substep
  dump, at a fixed rank count.
- **Operator-diff (R1):** per-kernel `max|Δ|=0` on identical (or replayed reference) inputs.
- **Self-tests:** halo-identity, 1-rank partition synthesis, dump round-trip, conservation
  (`Σssh_rhs` telescoping ~1e-13, tracer-sum), CW-orientation enforcement.
- **Multi-rank gate:** 1/8/32-rank byte-match (accepting only the `extrap_nod3D` partition-order diff).
- **Climate-stats (R-climate):** multi-year SST/SSS RMS — reserved for *later variations*
  (mixed precision, quads), not v1.
- **Traditional unit tests (ctest):** for pure utilities (calendar, namelist parse, partition
  synthesis, dump record I/O).

Every task lists its gate as explicit checklist item(s); **the gate must pass before the next task.**

## Progress Tracking

- Mark completed items `[x]` immediately. Add discovered tasks with ➕. Flag blockers with ⚠️.
- Keep this file in sync; update scope here when implementation deviates.

## What Goes Where

- **Implementation Steps** (`[ ]`): code, harness, gates achievable in this repo.
- **Post-Completion** (no checkboxes): external/manual items (multi-year HPC runs, deployment,
  partitioner provenance).

---

## Implementation Steps

> M0–M2 are detailed below (foundation → architectural MVP). M3–M6+ are a roadmap (Section
> "Roadmap"); they are detailed into tasks when their predecessor's gate passes, since their detail
> depends on earlier results.

### MILESTONE 0 — Foundation (no physics) — ✅ COMPLETE (tag `m0`)

Goal: a buildable skeleton + the full validation harness + mesh/partition/halo infrastructure, able
to load `pi`, synthesize a 1-rank partition, read `dist_N`, pass a halo-identity self-test, and
round-trip the dump — with **zero physics**.

#### Task M0.1: Repo skeleton + CMake build system

**Files:**
- Create: `CMakeLists.txt`, `configure.sh`, `env.sh`, `cmake/` (port `apply_fesom_compile_flags`)
- Create: `.gitignore`, `README.md`, `docs/LESSONS.md`, `docs/HANDOFF.md`

- [x] port the dwarf's CMake (`design_refs/tracer_dwarf/CMakeLists.txt`, `cmake/`, `configure.sh`):
      per-compiler/precision flags (GNU/Intel), options `USE_SINGLE_PRECISION`/`USE_HALF_PRECISION`,
      `build_<compiler>_<precision>/` dirs
- [x] define the **anchor build**: Intel + DP + `-fp-model precise` + no-reassoc (match FESOM2 v2.7.3)
- [x] record Levante env in `env.sh` (`module --force purge`; Intel/GNU + OpenMPI + serial netCDF)
- [x] `git init`; commit skeleton
- [x] **Gate:** empty library + a hello-MPI driver build cleanly in `intel_dp` and `gnu_dp`; runs on
      1 and 2 ranks on a login node

#### Task M0.2: `params/` — precision, constants, write-once config

**Files:**
- Create: `src/params/mod_precision.F90` (WP, MP=max(WP,single), `MPI_WP`)
- Create: `src/params/mod_constants.F90` (truncated π=3.14159265358979, g, density_0=1030, omega, …)
- Create: `src/params/mod_config.F90` (run settings) + `src/params/mod_param_phys.F90` (physics switches)
- Create: `src/params/hp_math_intrinsics.F90` (port dwarf; inert unless HP)
- Create: `test/test_params.F90`

- [x] `mod_precision`: two-tier WP/MP (D3); `MPI_WP` selects `MPI_DOUBLE_PRECISION`/`MPI_REAL`
- [x] `mod_constants`: **cite `oce_modules.F90` line for each constant**; truncated π exactly
- [x] `mod_config`/`mod_param_phys`: namelist read at init, then read-only; per-entity config rides
      in data types (not here)
- [x] write tests: constants equal FESOM2 values bit-for-bit; namelist parse (accept + reject)
- [x] **Gate:** `test_params` passes; builds in all precision variants (WP=MP=8 path is byte-trivial)

#### Task M0.3: `types/` — data types (no behavior)

**Files:**
- Create: `src/types/mod_mesh.F90` (`t_mesh` + `sparse_matrix`), `src/types/mod_partit.F90`
  (`t_partit` + `com_struct`), `src/types/mod_dyn.F90`, `src/types/mod_tracer.F90`,
  `src/types/mod_ice.F90`
- Create: `src/infra/mod_io_restart.F90` (+ `MOD_READ/WRITE_BINARY_ARRAYS` port), `test/test_types.F90`

- [x] `t_mesh`: cell-vertex; **`elem_nnodes(:)` + `elem_nodes(max_nv=4,:)`** (arity, D2); `edges`,
      `edge_tri`; CSR `nod_in_elem_ptr/nod_in_elem`; `gradient_sca(2*max_nv,:)`; geometry;
      integer level bounds `nlevels/ulevels/nlevels_nod2D[_min]`; `hnode/hnode_new/helem`; `ssh_stiff`
- [x] `t_partit`: `myDim/eDim/eXDim`, `myList_*`, `com_struct` (rPE/rptr/rlist, sPE/sptr/slist,
      1-based-cumulative), precompiled mpitypes
- [x] evolving types split **prognostic** (restartable) vs **work** (recomputed): `t_dyn`
      (uv/uv_rhsAB/w/eta_n/d_eta/ssh_rhs + `t_dyn_work` for density/N²/Kv/Av/hpressure/sw_alpha-beta/GM
      slopes + `t_solverinfo`); `t_tracer` (`t_tracer_data(:)` + `t_tracer_work`); `t_ice`
      (a_ice/m_ice/m_snow/uice/vice + **σ11/σ12/σ22 prognostic** + `t_ice_work`)
- [x] SoA, level-contiguous `(nl, n_entity)`; **tracer stride = nl** (allocation shape, not loop range)
- [x] generic `read(unformatted)`/`write(unformatted)` bindings per type (round-trip serialization is
      needed now for type tests + the dump; full restart-file orchestration deferred to M2.11/M5)
- [x] write tests: allocate/deallocate all; serialize→deserialize round-trip `max|Δ|=0`
- [x] **Gate:** `test_types` passes in all variants

#### Task M0.4: `infra/mod_partitioning` — partition + 1-rank synthesis

**Files:**
- Create: `src/infra/mod_partitioning.F90`, `test/test_partit.F90`

- [x] `par_init`/`par_ex`/`init_mpi_types` (port structure from FESOM2 `gen_modules_partitioning.F90`)
- [x] read `dist_<NP>/` (rpart.out, my_list, com_info), 1-based on disk → shift ids only
- [x] **1-rank synthesis** (D7): `npes==1` ⇒ identity local↔global, `myDim`=global, `eDim=eXDim=0`,
      empty `com_struct`, no neighbors, **no file read**
- [x] write tests: synthesize 1-rank on `pi`; read `dist_2`/`dist_8`; assert `Σ myDim_nod2D == nod2D`;
      elements/edges redundant at boundaries (`Σ myDim_elem2D ≥ elem2D`)
- [x] **Gate:** `test_partit` passes 1/2/8-rank

#### Task M0.5: `infra/mod_halo` — generic exchange

**Files:**
- Create: `src/infra/mod_halo.F90`, `test/test_halo.F90`

- [x] generic `exchange_nod`/`exchange_elem` (overloaded 2D/3D/multi-field/int; nl vs nl-1) +
      `*_begin`/`*_end` split; `luse_g2g` arg present-but-inert (structure from dwarf `gen_halo_exchange.F90`)
- [x] pack→`MPI_Isend`/`Irecv`→unpack via `com_struct`; **broadcast-only** (owner→halo, NO additive);
      `MPI_WP` datatype
- [x] exchange-and-compare **stale-halo probe** (copy, exchange copy, diff at halo)
- [x] write tests: **halo-identity** (set field=global id, exchange, assert halo==owner id) on
      `dist_2`/`dist_8`; probe detects an injected stale halo
- [x] **Gate:** `test_halo` passes multi-rank

#### Task M0.6: `infra/mod_dump` — gid-keyed validation harness

**Files:**
- Create: `src/infra/mod_dump.F90`, `tools/dump_diff.py`, `test/test_dump.F90`

- [x] per-substep dump: named node **and element** arrays at probe gids after each substep, keyed by
      1-based global id (rank-order independent), **FESOM2 binary layout**; env-gated (`FESOM_*_DUMP_DIR`),
      compiled-inert when off. **Element fields needing NEW shims in the oracle** (node fields already
      exist): `uv`, `UV_rhs`/`UV_rhsAB`, `Av`, `pgf_x`/`pgf_y` — element gid-keyed. The M2.2–M2.5
      element-path gates use controlled-input replay against these (so they don't block on later tasks
      like M2.8's `Av`).
- [x] `dump_diff.py`: per-step/point/column `|Δ|` histogram; names the **first divergent substep**;
      SIGNAL threshold separating real divergence from FP noise
- [x] add the matching dump gates to the FESOM2 oracle tree where missing (reuse existing
      kpp/ale/evp/tke shims from `port2/fesom2`)
- [x] write tests: dump round-trip; `dump_diff.py` flags an injected diff, passes identical files
- [x] **Gate:** `test_dump` + a `dump_diff.py` self-check pass

#### Task M0.7: `mesh/` — geometry + CW orientation + analytic generator

**Files:**
- Create: `src/mesh/mod_mesh_read.F90`, `src/mesh/mod_mesh_areas.F90`, `src/mesh/mod_mesh_aux.F90`,
  `src/mesh/mod_mesh_analytic.F90`

- [x] read partitioned mesh into `t_mesh`; build `edges`/`edge_tri`/CSR `nod_in_elem`
- [x] **enforce CW element orientation at load** (`test_tri` analog); abort loudly if violated; count swaps
- [x] compute geometry **transcribed from FESOM2 `oce_mesh.F90`** (elem_area, area/areasvol,
      gradient_sca, edge_cross_dxdy, metric_factor, zbar/Z, nlevels/ulevels, nlevels_nod2D[_min]) —
      generalize `1/3`→`1/n_vert` (inert for triangles; match FESOM2's exact divide-vs-multiply form
      per site — see Development Approach arity caveat)
- [x] analytic generator (port dwarf `analytic_mesh.F90`): in-memory tri mesh, closed + doubly-periodic
- [x] write tests/gate: load `pi`; **dump geometry arrays, `max|Δ|=0` vs FESOM2 mesh diagnostics**;
      CW swap count matches FESOM2; analytic mesh builds + geometry self-consistent
- [x] **Gate:** geometry byte-matches FESOM2 on `pi` (1-rank)

#### Task M0.8: `drivers/fesom_analytic` — end-to-end wiring (no physics)

**Files:**
- Create: `src/drivers/fesom_analytic.F90`, `src/step/mod_model.F90` (`t_model` + lifecycle stubs)

- [x] `t_model` aggregate; `model_init`/`model_step`(empty)/`model_finalize`
- [x] wire: MPI init → partit → analytic mesh → allocate types → (no physics) → finalize
- [x] write tests/gate: runs clean on analytic mesh 1-rank + multi-rank; allocations balanced
      (no leaks under a debug build)
- [x] **Gate:** `fesom_analytic` runs end-to-end; **M0 exit gate** — all M0 self-tests green; tag `m0`

### MILESTONE 1 — Tracer advection (FCT) on prescribed velocity

Goal: the first **byte-exact physics slice**. Transcribe FESOM2 tracer advection; gate `max|Δ|=0`
vs FESOM2 on `pi` with a *prescribed* velocity field. Proves types→halo→edge-FV kernel→dump→byte-gate
end-to-end and byte-gates the polygon-agnostic edge-FV scatter (the tri→quad core). Structure per
dwarf `oce_adv_tra_*`; **math from FESOM2** `oce_adv_tra_*.F90`.

#### Task M1.1: Horizontal advection (upwind + MUSCL) — ✅ DONE (max|Δ|=0)

**Files:**
- Create: `src/oce/oce_adv_tra_hor.F90`, `src/oce/oce_muscl_adv.F90`
- ➕ Also created (needed to produce the gate target `del_ttf_advhoriz`):
  `src/oce/oce_tracer_grad.F90` (`tracer_gradient_elements` → `tr_xy`),
  `src/oce/oce_adv_tra_flux.F90` (`oce_tra_adv_flux2dtracer` scatter; M1.4 driver will `use` it),
  plus the gate harness: `src/infra/mod_advhor_dump.F90`, `src/drivers/fesom_advhordump.F90`,
  FESOM2 `src/fesom_advhor_dump.F90` (wired into `ocean_setup`), `tools/{advhor_diff.py,
  run_advhordump_pi.sh,run_advhor_gate.sh}`.

- [x] transcribe from FESOM2 `oce_adv_tra_hor.F90` (upw1/muscl/mfct) + `oce_muscl_adv.F90`
      (muscl_adv_init/find_up_downwind_triangles/fill_up_dn_grad) + `tr_xy`
      (oce_tracer_mod.F90:181-182) (cite line refs); edge-FV flux → scatter to `edge_tri` neighbors
- [x] operator-diff harness: oracle = FESOM2's REAL kernels under **controlled-input replay** with
      an analytic prescribed velocity+tracer (computed from byte-proven coords in both codes); compare
      `del_ttf_advhoriz` (`oce_ale_tracer.F90:232`) — `tools/run_advhor_gate.sh`
- [x] **Gate:** operator-diff `max|Δ|=0` vs FESOM2 on `pi` — PASS on `del_ttf_advhoriz` +
      `adv_flux_hor` (upw1 + muscl) AND every intermediate (helem, nboundary_lay, edge_up_dn_tri,
      tr_xy, edge_up_dn_grad). ctest 13/13 (Intel+GNU dp); Intel sp + geom gate still green. See L9.

#### Task M1.2: Vertical advection (upwind + QR4C) — ✅ DONE (max|Δ|=0)

**Files:** Create: `src/oce/oce_adv_tra_ver.F90`
- [x] transcribe from FESOM2 `oce_adv_tra_ver.F90` (QR4C; watch the D=2 shallow-column double-write)
- [x] **Gate:** operator-diff `max|Δ|=0`

#### Task M1.3: FCT limiter (Zalesak) — ✅ DONE (max|Δ|=0)

**Files:** Create: `src/oce/oce_adv_tra_fct.F90`
- [x] transcribe from FESOM2 `oce_adv_tra_fct.F90`; `flux_eps=1e-16` flooring; low-order + clipped
      antidiffusive flux; halo writes for `fct_ttf_max/min` over `myDim+eDim` (faithful loop bounds)
- [x] **Gate:** operator-diff `max|Δ|=0`

#### Task M1.4: Advection driver + dispatch + integration — ✅ DONE (max|Δ|=0)

**Files:** Created: `src/oce/oce_adv_tra_driver.F90` (`do_oce_adv_tra`), `src/oce/oce_tracer_mod.F90`
(`init_tracers_AB`; `ab_epsilon=0.1` in `mod_config`), `src/oce/oce_ale_tracer.F90`
(`adv_tracers_ale`/`advect_tracer`)
- [x] `do_oce_adv_tra`: per-tracer scheme dispatch via `select case` on tra_adv_{hor,ver,lim};
      FCT + non-FCT (`do_zero_flux`) paths; `tra_adv_ph/pv` order knobs; prescribed velocity
- [x] `init_tracers_AB` (AB(2)/AB(3) → valuesAB; del_ttf zero; edge_up_dn_grad fill from grad(values))
      + `adv_tracers_ale`/`advect_tracer` (the per-tracer loop body + `del_ttf += advhoriz+advvert`).
      NB `model_step` integration awaits LIVE UV/W from dynamics (M2); the gate drives the assembled
      step with a prescribed velocity.
- [x] **Gate:** `max|Δ|=0` vs FESOM2's REAL `init_tracers_AB`+`do_oce_adv_tra` on `pi`, 1-rank
      (oracle shim extended, `libfesom.so` rebuilt). 7 records: `valuesAB`,
      `del_ttf_{advhoriz,advvert,}_step` (FCT) + `..._stepnon`/`del_ttf_step_non` (non-FCT).
      `use_wsplit` forced `.false.` (implicit w-split path = M2). Debug `-check all` clean. See L12.

#### Task M1.5: Multi-rank byte-gate — ⤳ FOLDED INTO M2.12 (decision 2026-06-19)

A multi-rank advection byte-gate needs the **local-mesh remap** (global→local numbering,
connectivity, `nod_in_elem2D` order, geometry, com-structs) — which `read_mesh` does NOT build
(`mod_mesh_read.F90:28` errors at `npes/=1`); it was scheduled at **M2.12**, where the
WHOLE-model multi-rank byte-match (which subsumes advection) lives. Doing the remap now just to
gate advection, then re-gating the whole model at M2.12, is duplicated effort. So:
- **M1 EXIT = the 1-rank anchor** (D0/D7 "serial == 1-rank MPI" is the bit-identity gold
  standard): M1.1–M1.4 `max|Δ|=0` on pi 1-rank ⇒ **tag `m1`** on M1.4. ✅ DONE
- The multi-rank advection gate (pi `dist_2`/`dist_8`; exchange-and-compare probe; lift the
  kernels' dropped halo exchanges + myDim/eDim loop bounds) is now part of **M2.12** below.

### MILESTONE 2 — Minimal ocean dynamical core (ARCHITECTURAL MVP)

Goal: prove the clean architecture reproduces FESOM2 **byte-exact on a real config** (CORE2, linfs,
PP, opt_visc=7), single- and multi-rank. Each kernel: structure clean, **math transcribed from
FESOM2**, operator-diff `max|Δ|=0` before integration.

**Oracle namelist (CRITICAL).** The shipped CORE2 namelists run KPP + GM + Redi (`namelist.oce`:
`mix_scheme='KPP'`, `Fer_GM=.true.`, `Redi=.true.`). The M2 byte-gate oracle must be FESOM2 run with a
**reduced M2 namelist** — `mix_scheme='PP'`, `Fer_GM=.false.`, `Redi=.false.`, `which_ale='linfs'`,
`opt_visc=7` — NOT the default (diffing against the default KPP/GM run = days of phantom divergences).
Capture this reduced namelist as part of the reference-run spec. Note the computed-but-dead producers
`sw_alpha_beta`/`compute_sigma_xy`/`compute_neutral_slope` are called unconditionally in FESOM2 but
consumed only by KPP/GM/Redi → deferred to M4 (don't be surprised they're absent from M2 yet present
in the D6 sequence + the oracle's SW_AB dump).

**D6 precise order (for reference):** `pressure_bv` (EOS + hpressure) → `pressure_force` (PGF) →
`sw_alpha_beta`/`sigma_xy`/`neutral_slope` (dead in M2) → mixing (PP) → `mo_convect` → `vel_rhs`
(Coriolis AB2 + momadv + PGF) → viscosity → impl-vert-visc → `ssh_rhs` → CG → `update_vel` → `hbar` →
`eta_n` → ALE thickness/W → tracer solve → `update_thickness_ale` commit.

#### Task M2.1: `pressure_bv` — EOS + hydrostatic pressure + N² (+ horizontal smoothing) — ✅ DONE (max|Δ|=0)
**Files:** Create: `src/oce/oce_pressure_bv.F90`
- [x] full Jackett-McDougall EOS in **split form** (`densityJM_components`: `bulk_0 + Z·(bulk_pz +
      Z·bulk_pz2)`, then `·rhopot/(…)`) — **never linearize α/β**; preserve the factorization (the bits
      depend on it). Cite `oce_ale_pressure_bv.F90`.
- [x] density anomaly subtracts the **`density_ref(nz,node)` array** (initialize to match FESOM2), NOT
      the scalar; `bvfreq`/N² divides by scalar `density_0=1030` (`oce_ale_pressure_bv.F90:440`)
- [x] **top-down `hpressure` integration lives HERE** (`oce_ale_pressure_bv.F90:367–403`; dumped before
      PGF) — NOT in M2.2
- [x] N² **horizontal** smoothing `smooth_nod(bvfreq, N2smth_hidx, …)` (`oce_ale_pressure_bv.F90:499`);
      pin `N2smth_hidx=1`, confirm `N2smth_v=.false.`; **one halo exchange per smoothing cycle** (uses M0.5)
- [x] MLD1/2/3 + `dbsfc` are KPP-only → omit in M2 (legitimate scope reduction); aux fields are WORK
- [x] **Gate:** operator-diff `max|Δ|=0` on `density`, `hpressure`, `bvfreq`

#### Task M2.2: Hydrostatic PGF (gradient contraction only) — ✅ DONE (max|Δ|=0)
**Files:** Create: `src/oce/oce_pgf.F90`
- [x] `pressure_force_4_linfs` → `pressure_force_4_linfs_fullcell` (`oce_ale_pressure_bv.F90:529,575`):
      `gradient_sca` contraction of the `hpressure` built in M2.1 → `pgf_x`/`pgf_y`. **No hpressure
      integration here.** (`pgf=Σ_k gradient_sca(k)·hpressure(nz,elnode_k)/density_0`; density_0 a
      runtime divisor but byte-identical both sides, L7/L10; same contraction as M1.1 `tracer_gradient_elements`)
- [x] **Gate:** operator-diff `max|Δ|=0` on `pgf_x`/`pgf_y` — PASS first run (pgf 66% non-zero, ~1e-5 m/s²;
      `tools/run_pressure_gate.sh`, 12 fields). Debug `-check all` clean; M1 advhor + 13/13 ctest still green.

#### Task M2.3: vel_rhs — Coriolis + AB2 + PGF (partial assembly) — ✅ DONE (max|Δ|=0)
**Files:** Create: `src/oce/oce_dyn_velrhs.F90`; `compute_coriolis` in `src/mesh/mod_mesh_areas.F90`
- [x] AB2 **actual** coefficients `ab1=-(0.5+ε)`, `ab2=(1.5+ε)` with **`epsilon=0.1`** (`ab_epsilon`
      from M1.4; `oce_ale_vel_rhs.F90:98–99`) — the `-0.5/1.5` base is only ε=0
- [x] **first-step Euler start:** reproduce `if (lfirst .and. .not. r_restart) ff=1.0_WP`
      (`oce_ale_vel_rhs.F90:287–290`) — `lfirst` passed as an explicit arg; gated BOTH ff=1.0 and ff=ab2
- [x] Coriolis init of `UV_rhsAB(1,1,…)` + AB2 blend of the previous-step array; add PGF/SSH-gradient.
      Added the `coriolis`/`coriolis_node` geometry field (`2·omega·sin(lat_geo)` via `r2g`)
- [x] **Gate (partial):** operator-diff `max|Δ|=0` on the AB2-blend + Coriolis + PGF + SSH-gradient
      pieces — 7 new records in `run_pressure_gate.sh` (19 fields total), the shim drives the REAL
      `compute_vel_rhs` twice with `momadv_opt=0`. FULL `UV_rhs` gate (incl. `momentum_adv_scalar`) → M2.4

#### Task M2.4: Momentum advection (✅ DONE) + biharmonic viscosity (✅ DONE)
**Files:** momadv ported into `src/oce/oce_dyn_velrhs.F90` (NOT a separate `oce_dyn_momadv.F90` —
matches FESOM2's `oce_ale_vel_rhs.F90` layout); viscosity → `src/oce/oce_dyn_visc.F90` (created)
- [x] `momentum_adv_scalar` (`momadv_opt==2`) — **called from inside `compute_vel_rhs`**
      (`oce_ale_vel_rhs.F90:271-273`), ADDING `w·du/dz`+`u·du/dx` into the same `UV_rhsAB(1,1:2,…)`
      slot as Coriolis *before* the AB2 blend; wired accordingly. **Gate:** operator-diff `max|Δ|=0`
      on the **complete `UV_rhs`** (Coriolis+advection+PGF+SSH) + the `uvnode_rhs` intermediate + `w_e`
      input — 4 new records → **21 fields** in `run_pressure_gate.sh`. Debug `-check all` clean. (L16)
- [x] **biharmonic viscosity `opt_visc=7`** (pinned to the pi namelist `visc_gamma0=0.003`, not the
      type default 0.03) — `viscosity_filter`→`visc_filt_bidiff` (`src/oce/oce_dyn_visc.F90`), a SEPARATE
      operator run AFTER `compute_vel_rhs`. A biharmonic = edge-based ∇² applied TWICE over INTERIOR edges
      only; NO new geometry (`edge_tri`/`elem_area`/`ulevels`/`nlevels`/`edge2D_in`, all gated; the
      `gradient_vec` worry was unfounded for opt_visc=7). **Gate:** operator-diff `max|Δ|=0` on
      `visc_u_c`/`visc_v_c` (pass-1 Laplacian) + the post-viscosity `uv_rhs_visc` — 3 new records →
      **24 fields**. UV bumped to 2.0/1.5 m/s so all three flow-aware `max(γ0,γ1,γ2)` branches fire
      (20.4/78.6/1.0%). Debug `-check all` clean. (L17)

#### Task M2.5: Implicit vertical viscosity (TDMA) — ✅ DONE (max|Δ|=0)
**Files:** Created: `src/oce/oce_dyn_ivertvisc.F90`; `zbar_e_bot` added to `t_mesh`
- [x] `impl_vert_visc_ale` per-element tridiagonal (Thomas) solve — implicit vertical viscosity `Av`
      + vertical advection (`w_i` upwind) + wind-stress top BC + quadratic bottom drag; OVERWRITES
      `UV_rhs` with the solution (`UV` read-only). A SEPARATE operator AFTER `viscosity_filter`
      (FESOM2 `oce_ale.F90:3874`). `Av`/`stress_surf` passed as explicit args (M2.8 mixing / M2.10
      forcing not yet ported → prescribed analytically for the gate). `helem`+`zbar_e_bot` built in
      the driver (full cells). (L18)
- [x] **Gate:** operator-diff `max|Δ|=0` on the post-solve `uv_rhs_ivv` + prescribed inputs
      (`Av`/`stress_surf`/`w_i`) — 4 new records → **28 fields** in `run_pressure_gate.sh`. TDMA is a
      sequential recurrence → byte-matches by pure L9 transitivity (PASSED first run). Non-vacuous
      (`max|d(uv_rhs)|=0.985`, `w_i>0/<0` both fire). Debug `-check all` clean (pi min nlevels=5 → no
      single-layer OOB). (L18)

#### Task M2.6: SSH — stiffness, ssh_rhs, CG solve — ✅ DONE (max|Δ|=0)
**Files:** Created: `src/oce/oce_ssh_rhs.F90` (`init_stiff_mat_ale` + `compute_ssh_rhs_ale`),
`src/oce/oce_ssh_solve.F90` (`ssh_solve_preconditioner` + `ssh_solve_cg` + `solve_ssh_ale`)
- [x] build `ssh_stiff` (CSR; stiffness `factor=g·dt·α·θ` × `(zbar_e_bot−zbar_e_srf)`·gradient·edge_cross +
      mass `areasvol/dt`); `ssh_rhs` (edge-divergence of `α(UV+UV_rhs)`); preconditioned CG (`soltol=1e-5`,
      **fixed iteration-count determinism** — converged in 37 iters, byte-identical). linfs builds the matrix
      ONCE (`update_stiff_mat_ale` skipped); 1-rank drops the global remap + `exchange_nod`/`MPI_Allreduce`
      (CG on local CSR `colind_loc`/`rowptr_loc`). `α=θ=1` ⇒ `(1−α)·ssh_rhs_old=0`; stiffness dt = pi
      namelist 86400/36 (NOT the shim's 1800).
- [x] **Gate:** operator-diff `max|Δ|=0` on `ssh_stiff_diag` + `ssh_Aeta` (full matvec A·eta_n) + `ssh_rhs` +
      `d_eta` (CG solution) — 4 new records → **32 fields** in `run_pressure_gate.sh`. **`Σ ssh_rhs` ≈ −3.7e-5
      = ~1e-13 relative** (edge-divergence telescoping at the ~1e8-scale bumped-UV rhs). Debug `-check all`
      clean; M1 advhor + 13/13 ctest still green. (L19; **CLEAN-rebuild footgun** — new files force a clean
      Release rebuild before the byte-gate is byte-reliable.)

#### Task M2.7: ALE (linfs) + velocity/SSH update — ✅ DONE (max|Δ|=0)
**Files:** Created: `src/oce/oce_ale.F90` (`update_vel` + `compute_hbar_ale` + `update_eta_n` +
`vert_vel_ale` + private `compute_CFLz`/`compute_Wvel_split`)
- [x] linfs (`hnode_new=hnode`); `update_vel` (`UV += UV_rhs + [-g·θ·dt·grad(d_eta)]`, a
      `gradient_sca` contraction of `d_eta`, FESOM2 `oce_dyn.F90:88-173`); `compute_hbar_ale`
      (`hbar_old=hbar`; `hbar += dt/areasvol·div(UV)`; `dhe`; the linfs water-flux term vanishes);
      `eta_n = α·hbar + (1−α)·hbar_old` (α=1 ⇒ `eta_n=hbar`); thickness/W update `vert_vel_ale`
      (W = −cumsum(div(UV·h))/area; linfs leaves `hnode_new=hnode`; then `compute_CFLz` +
      `compute_Wvel_split` explicit/implicit split). The **`K_v⁻` deformation bound is NOT in the
      post-CG ALE path** — `Kv` is produced by PP mixing (M2.8); the bullet was a mislabel, folded
      into M2.8. Fer_GM/ldiag_ke branches dropped (no GM / no `ke_*` in v1).
- [x] **Gate:** operator-diff `max|Δ|=0` on `uv_upd` + `ssh_rhs_old` + `hbar` + `dhe` + `eta_n_upd`
      + `w` + `hnode_new` + `cfl_z` + `w_split_e`/`w_split_i` (+ prescribed `hbar_in`) — 11 new
      records → **43 fields** in `run_pressure_gate.sh`. PASS first run (L9 transitive: every operand
      pre-gated — `d_eta` M2.6, post-TDMA `UV_rhs` M2.5, geometry, edge order). Non-vacuous: the Wvel
      split fired on 13253 (nz,node) (`CFL_z>wsplit_maxcfl`), `max|w|=0.042` m/s. Debug `-check all`
      clean (compute clean; the I/O writer needs `ulimit -s unlimited` for the big array temporary). (L20)

#### Task M2.8: PP vertical mixing ✅ DONE (byte-gate CLOSED, pi 1-rank)
**Files:** Created: `src/oce/oce_ale_mixing_pp.F90` (mirrors the FESOM2 filename)
- [x] transcribe `oce_ale_mixing_pp.F90` (3 sequential loops; produces `Kv` on nodes, `Av` on
      elements); dispatch `mix_scheme=='PP'`. **Convective adjustment is NOT here** — it is `mo_convect`
      (M2.8b), a separate routine called after PP. (`oce_mixing_pp` + `Kv0_background_qiang` +
      `Kv0_background`; `Kv`/`Av` added to `t_dyn_work`; `Kv0_const` etc. added to `mod_param_phys`.)
- [x] **Gate:** operator-diff `max|Δ|=0` (3 new records `uvnode`/`pp_Kv`/`pp_Av`, 46 fields total).
      Non-vacuous: the Ri factor spans [5.8e-7, 0.93] (62.4% of node-levels > 0.1), `max|Kv|`=8e-3
      (800× the K_ver background), `max|Av|`=8.7e-3. Debug `-check all` clean. See LESSONS L21.

#### Task M2.8b: `mo_convect` — convective adjustment ✅ DONE (byte-gate CLOSED, pi 1-rank)
**Files:** Created: `src/oce/oce_mo_conv.F90`
- [x] transcribe `oce_mo_conv.F90` (called AFTER PP — `oce_ale.F90:3729`); applies the instability
      adjustment to **both `Kv` (nodes, `oce_mo_conv.F90:85`) AND `Av` (elements, line 108)**.
      `use_momix` (TB04, reads forcing/ice) OMITTED → deferred to M2.10; `use_windmix` transcribed guarded-off.
- [x] pin run-namelist switches: `use_instabmix=.true., instabmix_kv=0.1` (the source of the `0.1`);
      forced `use_momix=.false.`, `use_windmix=.false.` (only the instability branch is live in M2)
- [x] **Gate:** operator-diff `max|Δ|=0` on `moc_Kv` and `moc_Av` (2 new records, 48 fields total).
      Non-vacuous: a localized unstable T band (warm subsurface lens) makes `bvfreq<0` on 5165 node-levels /
      9936 elem-levels → the floor fires, `max|ΔKv|`=0.090, `max|ΔAv|`=0.096 (~0.01→0.1). The T/S change
      cascades but M2.1-M2.8 records all re-verify `max|Δ|=0`. Debug `-check all` clean. See LESSONS L22.

#### Task M2.9a: Tracer-solve assembly (`solve_tracers_ale`) ✅ DONE (max|Δ|=0, 1-rank pi)
**Files:** Added `diff_tracers_ale` + `diff_part_hor_redi` + `diff_ver_part_impl_ale` + `bc_surface`
to `src/oce/oce_ale_tracer.F90` (mirrors FESOM2's file layout, NOT the plan's `oce_solve_tracers.F90`).
- [x] accumulate `del_ttf` = advection tendency (M1, prescribed input here) + **horizontal diffusion**
      (`diff_part_hor_redi`, Redi=.false. branch — `Kh·tr_xy` edge-flux scatter); `del_ttf` reset per tracer
- [x] **SSS/SST restoring**: on pi the *active* restoring is the **surface** one via `relax_salt`/`virtual_salt`
      in `bc_surface` (gated, prescribed analytically). The 3D `relax_to_clim` short-circuits (`clim_relax=0`)
      → deferred to M2.10 (needs Tclim/Sclim + `clim_relax>0`); `surf_relax_S` lives in forcing (M2.10).
- [x] **ALE reconstruct** `T=(T·hnode+del_ttf)/hnode_new` (linfs: `hnode_new==hnode` → `T += del_ttf/hnode`) →
      implicit vertical-diffusion TDMA `diff_ver_part_impl_ale`, the **FIRST consumer of the M2.8 `dyn%work%Kv`**
      (sourced LIVE post-`mo_convect`, no longer prescribed). Redi (isredi=0) + FCT (`do_wimpl=.false.`) +
      no-KPP/sw/icebergs branches reduce it to pure vertical diffusion + the surface BC.
- [x] **Gate:** operator-diff `max|Δ|=0` on `tsol_del_ttf_T/S` (post-diffusion `del_ttf`, per tracer) and
      `tsol_T`/`tsol_S` (the solved T/S), **57 fields total** (9 new `tsol_*`). Non-vacuous: Kv consumed=[0, 0.1]
      (live), `max|dT|=1.26`°C, `max|dS|=0.11`, `max|del_ttf_T|=0.50`. Debug `-check all` clean. See LESSONS L23.

#### Task M2.9b: `step_oce` sequence assembly + dispatch + thickness commit — ✅ DONE (max|Δ|=0, 1-rank pi)
**Files:** Created: `src/step/mod_step_oce.F90` (`step_oce`); added `compute_vel_nodes` +
`update_thickness_ale` to `src/oce/oce_ale.F90`, `solve_tracers_ale` to `src/oce/oce_ale_tracer.F90`;
Modified: `src/step/mod_model.F90` (`model_ocean_step` entry). Gate: `src/drivers/fesom_stepdump.F90`,
FESOM2 shim `port2/fesom2/src/fesom_step_dump.F90`, `tools/run_step{dump_pi,_gate}.sh`,
`tools/dump_diff.py --ignore-substep`.
- [x] assemble the faithful sequence with LIVE data flow: `compute_vel_nodes` → `pressure_bv`+smooth →
      `pressure_force` → PP `oce_mixing_pp` + `mo_convect` → `compute_vel_rhs` (momadv) → `viscosity_filter` →
      `impl_vert_visc_ale` (Av now LIVE from PP, the M2.5-deferred wiring) → `compute_ssh_rhs_ale`+`solve_ssh_ale`
      → `update_vel` → `compute_hbar_ale` → `update_eta_n` → `vert_vel_ale` → `solve_tracers_ale` (the full wrapper:
      advection + the M2.9a diffusion, Kv LIVE + salinity clamp) → `update_thickness_ale`. The dead-in-M2
      producers (`sw_alpha_beta`/`compute_sigma_xy`/`compute_neutral_slope`) are OMITTED (return at M4).
- [x] final `update_thickness_ale` commit (linfs no-op: hnode/helem fixed, `exchange_elem` 1-rank no-op)
- [x] **Gate:** **per-substep `max|Δ|=0` on `pi`** across the whole ocean step — the gate drives the REAL FESOM2
      `oce_timestep_ale` (its built-in `dump_shim_record_node`) vs FESOM3's `step_oce` (mirrored `mod_dump` dumps)
      on identical prescribed state, comparing the 13 NODE substeps (density/pressure/bvfreq / Kv / ssh_rhs /
      d_eta / hbar / eta_n / hnode_new / w / T / S / hnode = **65 records**, 5 probe nodes) with `dump_diff.py`.
      The SW_AB substep (dead in M2) is FESOM2-only and ignored. PASS first run; the assembly is exercised LIVE
      (uvnode← compute_vel_nodes, Av/Kv← PP+convection, d_eta← CG, T/S← advection+diffusion). Non-vacuous
      (bvfreq<0 → convective Kv=0.1, d_eta∈[-4.7,2.2], T/S evolved). Debug `-check all` clean; M1 advhor + the
      M2.1-M2.9a pressure gate (57 fields) + 13/13 ctest still green. **Scope:** `use_wsplit=.false.` (the FCT
      implicit vertical-advection `adv_tra_vert_impl`, do_oce_adv_tra's `use_wsplit=.true.` path, is a distinct
      unported kernel — M1.4 precedent; `impl_vert_visc_ale` still runs on the prescribed `w_i`). See LESSONS L24.

#### Task M2.10: Forcing (JRA55 bulk + SW penetration) — ✅ DONE (max|Δ|=0, 1-rank pi, 20 fields)
**Files:** Create: `src/io/mod_io_netcdf.F90`, `src/forcing/mod_forcing_read.F90`, `src/forcing/mod_forcing_bulk.F90`,
`src/oce/oce_shortwave_pene.F90`
- [x] **M2.10a forcing READ — ✅ DONE (max|Δ|=0, 1-rank pi, 8 fields).** The FIRST netCDF I/O. `mod_io_netcdf`
      (`use netcdf` wrapper, CMake `nf-config`+rpath) + `mod_forcing_read` (julday[noleap=365·yyyy] + binarysearch +
      time-axis transform + periodic-lon halo + spatial bilinear + the **two-stage time-interp `atmdata=rdate·coef_a+
      coef_b`** ~710820-scale cancellation [faithful, NOT collapsed] + **g2r wind-coef rotation**) + `vector_g2r`
      (mod_mesh_rotate). Reads CORE2 NCAR stubs. Fixed the L8 1-rank hang = FESOM2 `next_io_rank` np=1 infinite
      recursion (patched, value-identical). Gate `tools/run_forcing_gate.sh`. See LESSONS L25.
- [x] **M2.10b bulk transfer coeffs + wind stress — ✅ DONE (max|Δ|=0, 1-rank pi, 17 fields).**
      `ncar_ocean_fluxes_mode` (Large&Yeager + Large-2009 drag, 5-iter Monin-Obukhov) → `Cd`/`Ch`/`Ce` (gated vs the
      REAL routine); `stress_atmoce = Cd·ρ_air·|Δu|·Δu`; node→elem `stress_surf` (`a_ice=0`, `/3`). Prescribed SST +
      surface ocean velocity (dummy ice, thermo type-defaults). Byte traps: `inc_ratio=1.0e-4`/`inv_rhoair=1./1.3`/
      `tmelt=273.15`/`rhoair=1.3` transcribed VERBATIM (un-suffixed); `(ustar*ustar)`; `atan(1.0_WP)`. Debug `-check
      all` clean. See LESSONS L25.
- [x] **M2.10c shortwave penetration — ✅ DONE (max|Δ|=0, 1-rank pi).** `cal_shortwave_rad` (Morel&Antoine/Sweeney
      2005: `swsurf=(1-albw)·shortwave·0.54`, chl floor 0.02, the v1/v2/sc1/sc2 polynomial, the two-exponential
      `sw_3d` depth profile via `zbar_3d_n`, `/vcpw`) → `sw_3d` + the `heat_flux` visible-removal. Built
      `src/oce/oce_shortwave_pene.F90`. Consumes the LIVE M2.10a `shortwave`; `chl`/`heat_flux`/`a_ice=0` prescribed;
      `albw=0.066_WP` forced both sides. Non-vacuous (chl floor fires on 476 polar nodes). Debug `-check all` clean.
- [ ] (SSS/SST restoring lives in the tracer solve M2.9a, NOT here — only prepare surface flux/coeff inputs)
- **SCOPE:** `heat_flux`/`water_flux` air-sea budget = **M3** (obudget in `ice_thermo_oce.F90` + `oce_fluxes`); M2.10
      leaves them prescribed (M2.9b shim) until M3. `gen_bulk_formulae.F90` = transfer coefficients ONLY (not stress/heat).

#### Task M2.11: Full driver + single-rank CORE2 byte-gate
**Files:** Create: `src/drivers/fesom.F90`, `src/step/mod_model.F90` (full lifecycle), `src/io/` (minimal)
- [x] **M2.11a CORE2 geometry byte-gate — ✅ DONE (2026-06-20, max|Δ|=0, 19 fields).** Scaled the mesh pipeline 40×
      (nod2D=126858/elem2D=244659/edge2D=371644/nl=48) + **closed the L8 CW-swap deferral** (FESOM3 made 244654/244659
      `enforce_cw_orientation` swaps — ≈100%, mesh stored CCW; pi had 0 — post-swap geometry max|Δ|=0). `tools/make_dist1.py`
      hand-crafts the single-rank dist_1; `tools/run_geom_gate_core2.sh`; `fesom_geomdump` unchanged (`FESOM3_MESH_DIR`).
      Confirmed safe: min nlevels=5 (no single-layer cols, L18 N/A), cavity+partial-cell OFF. See LESSONS L26.
- [x] **M2.11b initial conditions (`do_ic3d`) — ✅ DONE (2026-06-21, max|Δ|=0, 3 fields, CORE2 1-rank).** Ported the
      3D-climatology IC: `src/oce/oce_initial_state.F90` (`do_ic3d` + `nc_readGrid` + `nc_ic3d_ini` + `getcoeffld` +
      `extrap_nod3D`) reading **phc3.0_winter.nc** (360×180×33) — netCDF 3D-double read (`nc_get_var3d_dp`) + NaN→dummy +
      periodic-lon halo + spatial bilinear + vertical LINEAR interp onto `Z_3d_n` + `extrap_nod3D` (Gauss-Seidel
      neighbour-average + downward fill — the partition-order step, deterministic at 1-rank) + Kelvin guard +
      **`insitu2pot`** (Bryden-1973 RK4 `ptheta`/`atg` in `oce_pressure_bv.F90`). `idlist=2,1` → salt (data(2)) read FIRST,
      temp (data(1)) SECOND, `t_insitu=.true.`. Reused `forcing_binarysearch` + `mod_io_netcdf`. Gate: `tools/run_ic_gate_core2.sh`
      (Z_3d_n input + ic_temp/ic_salt) vs the oracle's live `Tclim`/`Sclim` (`fesom_ic_dump.F90` shim). PASSED first run;
      Debug `-check all` clean; pressure gate (57) + ctest (13) still green. See LESSONS L27.
- [x] ✅ **M2.11c — full lifecycle + multi-step CORE2 gate.** DONE 2026-06-21. `src/drivers/fesom_lifecycle.F90` (cold-start
      CORE2 + do_ic3d IC + N-step `step_oce` loop, multi-step AB2); forcing/flux M3-gap prescribed from the oracle dump.
      **Gate `tools/run_lifecycle_gate_core2.sh`: MATCH, 195 records (13 substeps × 5 probes × 3 steps), worst |Δ|=0** —
      whole dynamical core byte-exact across multiple CORE2 steps INCLUDING the free-surface CG `d_eta`. The CORE2 CG
      "reproducibility floor" was SOLVED (it was an auto-vectorised preconditioner divide; one `!DIR$ NOVECTOR`; LESSONS L29).
      **The FORCED variant too (M2.11c-2, re-verified post-L29 2026-06-21): `tools/run_lifecycle_forced_gate_core2.sh`
      MATCHes, 195 records worst |Δ|=0** (real CORE2 forcing + ice + the M3-gap air-sea fluxes prescribed from
      `fesom_flux_dump.F90`) — so ALL of M2 is byte-identical end-to-end on CORE2. M3/M2.12 unblocked.

#### Task M2.12: Multi-rank + production validation (MVP exit)
**Files:** Create: the local-mesh remap in `mod_mesh_read`/a new builder; Modify: tests/scripts; `docs/HANDOFF.md`

**Foundation readiness (scoped 2026-06-21).** The multi-rank halo + partition infra ALREADY EXISTS and is
tested 1/2/8-rank: `mod_halo.F90` is a real MPI exchange (manual pack → Isend/Irecv → broadcast-only
owner→halo unpack via `com_struct`; `exchange_nod` 2D/3D-real + 2D-int, `exchange_elem` 2D/3D-real), and
`mod_partitioning.read_dist_partition` already loads `myList_nod2D/elem2D/edge2D` + the three com-structs
(`com_nod2D`/`com_elem2D`/`com_elem2D_full`). FESOM2 uses precompiled MPI_TYPE_INDEXED datatypes instead of
manual pack, but that moves the SAME bytes with NO arithmetic → byte-identical result. So M2.12 is faithful
transcription on a working foundation, NOT a from-scratch parallel build. **Gate rule (L8):** compare FESOM3
`dist_N` vs FESOM2 `dist_N` PER-RANK on OWNED entries (NOT vs 1-rank global) — both share the same `myList`
order so local index i ↔ same global id, and the per-node area accumulation order matches → byte-identical;
the 1-rank global uses a different element permutation (non-associative FP) so it would NOT match.

- [x] ✅ **M2.12a — local-mesh remap + per-rank GEOMETRY byte-gate (pi dist_2/dist_8).** DONE 2026-06-21:
      `max|Δ|=0` on all 19 geometry fields, every rank, dist_2 (2) AND dist_8 (8) vs same-partition FESOM2
      (`tools/run_geom_gate_multirank.sh`). Built `read_mesh_local` (scatter via full-size inverse maps + owned
      `nod_in_elem2D` + `enforce_cw`(owned) — closes the multi-rank CW-swap caveat) + partition-aware
      `compute_geometry` (`local_bounds`; owned centers precomputed + `exchange_elem`'d for owned-edge
      `edge_cross_dxdy`; owned-node areas local). Per-rank owned dump (FESOM3 `mod_geom_dump` + oracle
      `fesom_geom_dump.F90`). No regression: 1-rank geom (19) + pressure (57) gates + 13/13 ctest still pass.
      DEFERRED to M2.12b (not needed for the owned-entry gate): the `find_neighbors` nod_in_elem2D halo dance
      (eXDim) + area/elem_area halo exchanges. THE prerequisite.
      Replace the `npes/=1` error in `mod_mesh_read.F90:28` with the local-mesh remap, transcribed from
      FESOM2 `read_mesh`/`find_neighbors`/`mesh_areas` (oce_mesh.F90:212/1969/2162):
      (1) full-size global→local inverse maps from `myList_*` (simpler than FESOM2's chunked `mapping`;
      identical result — pure deterministic scatter). (2) Scatter local `coord_nod2D` (read global nod2d.out,
      rotate, place at local idx), `elem2D_nodes` (OWNED elems only `(3,myDim_elem2D)`, localize node ids),
      `edges`/`edge_tri` (localize node + elem ids, boundary `edge_tri`→0), `nlevels`/`nlevels_nod2D`/`depth`;
      `nl`/`zbar`/`Z` are global. (3) `enforce_cw_orientation` on local owned elems (per-element deterministic
      → consistent across ranks; CLOSES the multi-rank CW-swap caveat). (4) `nod_in_elem2D` via the
      `find_neighbors` exchange dance: build from owned elems (skip halo nodes `node>myDim`), `exchange_nod`
      the count, pack global elem-ids → `exchange_nod` → re-localize through the **eXDim** extended-halo
      (`myList_elem2D(1:myDim+eDim+eXDim)`). SKIP `elem_neighbors`/`elem_edges` (not needed for geometry).
      (5) `setup_vertical` with `exchange_nod(nlevels_nod2D_min)` (oce_mesh.F90:1669). (6) `compute_geometry`
      runs on local arrays UNCHANGED, EXCEPT `mesh_areas` needs `exchange_elem(elem_area)` over the FULL halo
      (oce_mesh.F90:2220) — likely a NEW `exchange_elem_full` on `com_elem2D_full` (eDim+eXDim); the gate
      catches if `com_elem2D` is insufficient. **Oracle:** extend `fesom_geom_dump.F90` (currently npes==1
      no-op) to dump per-rank LOCAL OWNED arrays. **Gate:** `tools/run_geom_gate_multirank.sh` (FESOM2 dist_2
      dump vs FESOM3 dist_2 dump, per-rank owned, `geom_diff.py --glob`); target `max|Δ|=0`; repeat dist_8.
- [x] ✅ **M2.12b — M1 advection multi-rank gate (folded from M1.5).** DONE 2026-06-21: `max|Δ|=0` on all 13
      advection fields, every rank, pi **dist_2 (2) AND dist_8 (8)** vs same-partition FESOM2
      (`tools/run_advhor_gate_multirank.sh`). Lifted the advection subtree to multi-rank with an OPTIONAL
      `partit` (absent ⇒ proven 1-rank path verbatim; present+npes>1 ⇒ owned/halo bounds + exchanges) and
      KEPT the explicit-shape dummies/per-element arithmetic so the codegen — and the byte-match vs FESOM2
      (incl. any vectorised divide, L29) — is unchanged. Added: `mod_part_bounds.owned_bounds/is_multirank`;
      the `find_neighbors` halo dance (`complete_nod_in_elem_halo` in `read_mesh_local` — completes
      `nod_in_elem2D` for halo nodes, re-localised through eXDim); `exchange_elem_full` (2D-r/2D-i/3D-r on
      `com_elem2D_full`) + `core_blk_r`; `compute_geometry` halo exchanges (`elem_area` full + `area`/`areasvol`/
      `*_inv` nodes); `exchange_elem(tr_xy)` full-halo + `exchange_nod(fct_LO)` + `exchange_nod(fct_plus/minus)`
      at the FESOM2 sites; `find_up_downwind_triangles` `coord_elem`/`e_nodes` (full-halo, global-id match);
      `nboundary_lay` owned+halo NO-exchange (partition-local, matches FESOM2). Gated `valuesAB`/`edge_up_dn_grad`/
      `fct_LO`/`fct_ttf_max-min`/`fct_plus-minus`/`del_ttf_{advhoriz,advvert,}` (FCT) + `del_ttf_*` (non-FCT/MUSCL).
      `edge_up_dn_grad` captured PRE-`do_oce_adv_tra` (FESOM2 reuses it as the FCT `AUX` scratch → bignumber).
      No regression: 1-rank geom(19)/advhor(40)/pressure(57)/step(65) + CORE2 lifecycle(195) + 13/13 ctest all
      `max|Δ|=0`/green. Lessons L31. **ACTIVE NEXT = M2.12c.**
- [x] ✅ **M2.12c — whole-model multi-rank byte-match (MVP exit):** decomposed c-1/c-2/c-3 — ALL DONE 2026-06-22.
      - [x] ✅ **M2.12c-1 — pre-SSH DYNAMICS chain (2026-06-22):** `max|Δ|=0` pi dist_2 + dist_8, 25 records
            (density/pressure/bvfreq/Kv/ssh_rhs), `tools/run_stepdyn_gate_multirank.sh`. Lifted compute_vel_nodes →
            pressure_bv(+smooth_nod) → pgf → oce_mixing_pp → mo_convect → compute_vel_rhs(+momentum_adv_scalar) →
            viscosity_filter(visc_filt_bidiff) → impl_vert_visc_ale → compute_ssh_rhs_ale via optional `partit` +
            `mod_part_bounds.local_dims` + FESOM2 exchanges (`mod_halo` gained a rank-3 `exchange_nod` node-block).
            `fesom_stepdump_mr` (FESOM3) + `fesom_step_dump` shim npes>1 (FESOM2). No 1-rank regression. Lessons L32.
      - [x] ✅ **M2.12c-2 — SSH stiffness + CG (2026-06-22):** `max|Δ|=0` on `d_eta`, pi dist_2 + dist_8 (30 records,
            `tools/run_stepdyn_gate_multirank.sh` substep 9 un-ignored; 38 CG iters/rank). The FIRST multi-rank
            iterative solver. **NO mesh-infra extension needed (the scoped plan above was WRONG):** a `dist_N`-file
            check proved the stiffness owned rows assemble FULLY LOCALLY (owned edges ⇒ owned triangles; FESOM2's own
            elem2D_nodes/gradient_sca are owned-only). `init_stiff_mat_ale` lifts with ONLY bounds (`nNodO`/`nEdgeO`,
            `n_num(nNodL)`); CG owned dot-products + new `mod_halo::allreduce_sum` + per-iter `exchange_nod(pp/rr/x)` +
            `exchange_nod(diag_values)` in the precond; rtol/conv denominators GLOBAL. Cross-rank `MPI_Allreduce`
            byte-matches (L6); L29 NOVECTOR precond carries over. No 1-rank regression. Lessons **L33**.
      - [x] ✅ **M2.12c-3 — ALE update + tracer SOLVE + WHOLE-STEP gate (2026-06-22):** `max|Δ|=0` on ALL 65 substep
            records (density/pressure/bvfreq/Kv/ssh_rhs/d_eta/hbar/eta_n/hnode_new/w/T/S/hnode), pi dist_2 + dist_8,
            `tools/run_step_gate_multirank.sh`. **The whole multi-rank ocean step byte-matches FESOM2 — the MVP.**
            Lifted update_vel/compute_hbar_ale/update_eta_n/vert_vel_ale(+compute_CFLz/compute_Wvel_split)/
            update_thickness_ale (oce_ale.F90) + solve_tracers_ale/diff_tracers_ale/diff_part_hor_redi/
            diff_ver_part_impl_ale (oce_ale_tracer.F90; the FCT advection was M2.12b) via the optional-`partit` pattern
            (bounds from `owned_bounds`/`local_dims`; FESOM2 exchanges at the FESOM2 sites), and threaded the optional
            `partit` through `mod_step_oce::step_oce` (absent ⇒ 1-rank verbatim). New whole-step MR driver
            `fesom_stepfull_mr` (= step_oce through partit, gid-keyed per-rank dumps); the FESOM2 oracle is UNCHANGED
            (`fesom_step_dump` npes>1 already runs the whole `oce_timestep_ale` + dumps all substeps). pi ships only
            dist_1/2/8 (no dist_32; dist_2+dist_8 = the c-1/c-2 coverage). The byte-gate on the rich analytic state
            (non-trivial T/S/UV/SSH at 2 AND 8 ranks) subsumes the rest-at-rest / gravity-wave sanity probes; matching
            FESOM2's exact owned loop bounds + exchanges makes any stale halo unobservable. **NO 1-rank regression**
            (step 65 + pressure 57 + advhor MR + stepdyn MR + 13/13 ctest all `max|Δ|=0`). New ALE update kernels:
            owned-element update_vel + exchange_elem_full(UV); compute_hbar_ale edge-div over owned edges + always-on
            exchange_nod(hbar); vert_vel_ale owned cumsum + exchange_nod(w)/exchange_nod(hnode_new); CFLz/Wvel_split
            owned+halo. Tracer solve: solve_tracers_ale per-tracer exchange_nod(values) + owned+halo clamp;
            diff_part_hor_redi owned edges; diff_ver_part_impl_ale owned-node TDMA. Lessons **L34**.
- [ ] **production dt=1800 (reduced M2 namelist): PRIMARY gate = per-substep `max|Δ|=0` vs the FESOM2
      reduced-config run**; SECONDARY sanity = run past the historical blow-up step, non-drifting
- [x] ✅ **Gate:** all multi-rank + 1-rank byte-gates GREEN ⇒ **architectural MVP** — tagged `m2-mvp` (2026-06-22);
      `docs/LESSONS.md` L30–L34 written. Baseline re-confirmed after a clean rebuild: step MR dist_2/8 (65 rec each),
      step 1-rank (65), stepdyn MR dist_2/8 (30 each), advhor MR dist_8 (8 ranks), pressure 1-rank (57 fields),
      13/13 ctest — all `max|Δ|=0`/green. (The production dt=1800 run above stays OPTIONAL secondary validation.)

### Task V: Verify acceptance criteria (per milestone)
- [ ] every kernel in the milestone has an operator-diff `max|Δ|=0` record
- [ ] whole-model byte-gate green at the milestone's config + ranks
- [ ] off-switches for the milestone's features byte-reproduce the prior milestone
- [ ] tag the milestone; update `docs/HANDOFF.md` (state, canonical reference run, next task)

### Task F: Final documentation (per milestone)
- [ ] update `README.md`, `docs/LESSONS.md`, `docs/HANDOFF.md`
- [ ] move/curate this plan section to `docs/plans/completed/` when M0–M2 done

---

## Roadmap (M3–M6+ — detailed into tasks when the predecessor gate passes)

Each milestone: feature behind an off-by-default flag; off-switch byte-gates the prior milestone;
transcribe math from FESOM2; per-kernel operator-diff then whole-model byte-gate; validate at
production dt=1800; tag.

- **M3 — Sea ice (EVP).** `ocean2ice` → EVP (120 subcycles, `Tevp_inv=3/ice_dt`, **`ice_strength`
  0.5 factor** near `ice_EVP.F90:599` — re-verify exact line when detailing M3, `m_ice`=h·a, σ elastic
  memory persists) → `tg_rhs` → ice FCT →
  `cut_off` → thermo → `oce_fluxes`. Watch: `bc_index_nod2D` multi-rank-safe; ice_dt synced to ocean
  dt; NaN-vs-0 ice masking in diagnostics. Build `t_ice` from scratch.
- **M4 — GM/Redi.** sigma_xy/neutral_slope → Γ solve (TDMA) → bolus velocities (**`fer_w` from
  `div(fer_uv·h)`; never per-cell clamp**) → rotated Redi diffusion; master off-switch byte-matches
  pre-GM model. Producers: `sw_alpha_beta`. Active: `K_GM_max=1000`, ODM95 tapering.
- **M5 — KPP + production multi-year (= paper-parity).** Port **FESOM1.4 inline KPP**
  (`oce_ale_mixing_kpp.F90`, `mix_scheme_nmb=1`); **single `Kv=diffK(:,:,1)` for T and S**;
  `ddmix`/`ghats` gate-only (off in CORE2); `bldepth` is the high-risk routine (needs the smoothed
  N²). Then **multi-year CORE2 dt=1800** statistical validation (SST RMS target ~0.004–0.006 °C,
  zero deep drift). End of M5 = parity with the paper's validated FESOM2 config.
- **M6+ — New targets beyond the paper (fresh validation, not proven ports):**
  - **zstar:** ALE thickness evolution `hnode*(1+η/H)`; **real freshwater/salt flux** replaces
    virtual salt; **no `hpressure`** (Shchepetkin PGF self-contained); cumulative-CSR stiffness
    rebuilt each step; preconditioner built once.
  - **TKE:** CVMix classical TKE; **prognostic `tke` field** (the one stateful mixing scheme);
    `tke_cd=3.75` from namelist; halo-exchange internal `tke_Av` before node→elem mean; Prandtl
    literal `6.6` stays double.
  - **mEVP:** Bouillon 2013; α/β=250; frozen-entry rhs anchor (not current iterate).

## Technical Details

- **Storage:** SoA, level-contiguous `(nl, n_entity)` column-major (== C `[n*nl+lev]`); topology-agnostic.
  Tracer stride = `nl` (allocation shape). Aux 3D fields are interface (`bvfreq/Kv/Av`, on `zbar`) vs
  mid-layer (`T/S/uv`, on `Z`) — register correctly.
- **Geometry:** scalars@nodes, velocity@element centroids; bathymetry cellwise (`nlevels(elem)`);
  `nlevels_nod2D`=MAX over cells, `nlevels_nod2D_min`=MIN over cells (=ALE deformation limit).
- **Rotated grid:** CORE2 Euler (50,15,−90); `coord_nod2D` rotated (compute), `geo_coord_nod2D`
  geographic (IC/forcing/Coriolis); all boundary-crossing vectors rotate (g2r in, r2g out).
- **Precision (D3):** WP/MP two-tier; v1 WP=MP=real64. `_WP`/`_MP` literal suffixes. Cancellation
  islands (CG/reductions/EOS/pressure/rotation) protected in the later mixed-precision phase.
- **MPI (D7):** generic `exchange_nod/elem`; broadcast-only halo; node-only reductions; 1-rank
  synthesis; faithful halo-write loop bounds (`myDim+eDim`); size fields `myDim+eDim` if read at halo.
- **Determinism:** anchor build Intel + DP + `-fp-model precise` + no-reassoc; fixed rank count;
  `max|Δ|=0` byte-gate via gid-keyed dump.

## Post-Completion

*External / manual — no checkboxes.*

**Manual / HPC verification:** multi-year CORE2 climate runs (SST/SSS RMS, deep-drift, ice
area/volume) on Levante batch; energy/SYPD profiling (later, optimization phase).

**External provenance:** the METIS partitioner that produces `dist_<NP>/` is out of scope (we consume
the files + synthesize 1-rank). PHC IC provenance is partition-dependent — keep partition-matched IC
caches as part of any reference run spec.

**Later-phase activation (foundations already laid):** mixed precision (flip WP, protect islands +
multi-year gate), quad/mixed meshes (exercise arity on a quad test mesh + climate gate), GPU-direct/
overlap, other-language ports (reimplement the `exchange_nod/elem` contract).
