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

#### Task M2.1: `pressure_bv` — EOS + hydrostatic pressure + N² (+ horizontal smoothing)
**Files:** Create: `src/oce/oce_pressure_bv.F90`
- [ ] full Jackett-McDougall EOS in **split form** (`densityJM_components`: `bulk_0 + Z·(bulk_pz +
      Z·bulk_pz2)`, then `·rhopot/(…)`) — **never linearize α/β**; preserve the factorization (the bits
      depend on it). Cite `oce_ale_pressure_bv.F90`.
- [ ] density anomaly subtracts the **`density_ref(nz,node)` array** (initialize to match FESOM2), NOT
      the scalar; `bvfreq`/N² divides by scalar `density_0=1030` (`oce_ale_pressure_bv.F90:440`)
- [ ] **top-down `hpressure` integration lives HERE** (`oce_ale_pressure_bv.F90:367–403`; dumped before
      PGF) — NOT in M2.2
- [ ] N² **horizontal** smoothing `smooth_nod(bvfreq, N2smth_hidx, …)` (`oce_ale_pressure_bv.F90:499`);
      pin `N2smth_hidx=1`, confirm `N2smth_v=.false.`; **one halo exchange per smoothing cycle** (uses M0.5)
- [ ] MLD1/2/3 + `dbsfc` are KPP-only → omit in M2 (legitimate scope reduction); aux fields are WORK
- [ ] **Gate:** operator-diff `max|Δ|=0` on `density`, `hpressure`, `bvfreq`

#### Task M2.2: Hydrostatic PGF (gradient contraction only)
**Files:** Create: `src/oce/oce_pgf.F90`
- [ ] `pressure_force_4_linfs` → `pressure_force_4_linfs_fullcell` (`oce_ale_pressure_bv.F90:529`):
      `gradient_sca` contraction of the `hpressure` built in M2.1 → `pgf_x`/`pgf_y`. **No hpressure
      integration here.**
- [ ] **Gate:** operator-diff `max|Δ|=0` on `pgf_x`/`pgf_y`

#### Task M2.3: vel_rhs — Coriolis + AB2 + PGF (partial assembly)
**Files:** Create: `src/oce/oce_dyn_velrhs.F90`
- [ ] AB2 **actual** coefficients `ab1=-(0.5+ε)`, `ab2=(1.5+ε)` with **`epsilon=0.1`** (the dt=1800
      trap; cite `o_PARAM` / `oce_ale_vel_rhs.F90:98–99`) — the `-0.5/1.5` base is only ε=0
- [ ] **first-step Euler start:** reproduce `if (lfirst .and. .not. r_restart) ff=1.0_WP`
      (`oce_ale_vel_rhs.F90:287–290`) — omitting it diverges the step-1 byte-gate
- [ ] Coriolis init of `UV_rhsAB(1,1,…)` + AB2 blend of the previous-step array; add PGF/SSH-gradient
- [ ] **Gate (partial):** operator-diff `max|Δ|=0` on the AB2-blend + Coriolis + PGF pieces. The FULL
      `UV_rhs` gate is at M2.4 (because `momentum_adv_scalar` adds into the same slot before the blend)

#### Task M2.4: Momentum advection + biharmonic viscosity
**Files:** Create: `src/oce/oce_dyn_momadv.F90`, `src/oce/oce_dyn_visc.F90`
- [ ] `momentum_adv_scalar` (`momadv_opt` per run namelist) — note it is **called from inside
      `compute_vel_rhs`** (`oce_ale_vel_rhs.F90:273`), adding into the same `UV_rhsAB(1,1,…)` slot as
      Coriolis *before* the AB2 blend; wire accordingly
- [ ] **biharmonic viscosity `opt_visc=7`** (pin to run namelist, not default 5)
- [ ] **Gate:** operator-diff `max|Δ|=0` on the **complete `UV_rhs`** (Coriolis+advection+PGF) and on viscosity

#### Task M2.5: Implicit vertical viscosity (TDMA)
**Files:** Create: `src/oce/oce_dyn_ivertvisc.F90`
- [ ] Thomas solver; wind stress top BC + bottom drag
- [ ] **Gate:** operator-diff `max|Δ|=0`

#### Task M2.6: SSH — stiffness, ssh_rhs, CG solve
**Files:** Create: `src/oce/oce_ssh_rhs.F90`, `src/oce/oce_ssh_solve.F90`
- [ ] build `ssh_stiff` (CSR; negative factor `-g·dt·α·hbar`); `ssh_rhs`; CG (`soltol=1e-5`,
      **fixed iteration-count determinism**); preconditioner; exchange preconditioner diagonal
- [ ] **Gate:** operator-diff `max|Δ|=0`; **`Σ ssh_rhs` over owned nodes telescopes to ~1e-13**

#### Task M2.7: ALE (linfs) + velocity/SSH update
**Files:** Create: `src/oce/oce_ale.F90`
- [ ] linfs (`hnode_new=hnode`); `update_vel`; `hbar`; `eta_n`; thickness/W update; `K_v⁻` deformation bound
- [ ] **Gate:** operator-diff `max|Δ|=0`

#### Task M2.8: PP vertical mixing
**Files:** Create: `src/oce/oce_mix_pp.F90`
- [ ] transcribe `oce_ale_mixing_pp.F90` (3 sequential loops; produces `Kv` on nodes, `Av` on
      elements); dispatch `mix_scheme=='PP'`. **Convective adjustment is NOT here** — it is `mo_convect`
      (M2.8b), a separate routine called after PP.
- [ ] **Gate:** operator-diff `max|Δ|=0`

#### Task M2.8b: `mo_convect` — convective adjustment
**Files:** Create: `src/oce/oce_mo_conv.F90`
- [ ] transcribe `oce_mo_conv.F90` (called AFTER PP — `oce_ale.F90:3729`); applies the instability
      adjustment to **both `Kv` (nodes, `oce_mo_conv.F90:85`) AND `Av` (elements, line 108)**
- [ ] pin run-namelist switches: `use_instabmix=.true., instabmix_kv=0.1` (the source of the `0.1`);
      confirm `use_momix=.false.`, `use_windmix=.false.` (only the instability branch is live in M2)
- [ ] **Gate:** operator-diff `max|Δ|=0` on `Kv` and `Av`

#### Task M2.9a: Tracer-solve assembly (`solve_tracers_ale`)
**Files:** Create: `src/oce/oce_solve_tracers.F90`
- [ ] accumulate `del_ttf` = horizontal adv + explicit-vertical adv + diffusion (`del_ttf_advhoriz/advvert`,
      `oce_ale_tracer.F90:232,240`); `del_ttf` zeroed in init
- [ ] **SSS/SST restoring** `relax_to_clim`/`relax_2_tsurf` (`oce_ale_tracer.F90:255–264`) — it modifies
      T/S *inside* the solve; pin `surf_relax_S=1.929e-06`, `balance_salt_water=.true.` + virtual-salt balancing
- [ ] **ALE reconstruct** `T=(T·hnode+del_ttf)/hnode_new` (`oce_ale_tracer.F90:468–471`) → implicit
      vertical-diffusion TDMA `diff_ver_part_impl_ale` (line 482)
- [ ] **Gate:** operator-diff `max|Δ|=0` on `del_ttf` and the updated T/S

#### Task M2.9b: `step_oce` sequence assembly + dispatch + thickness commit
**Files:** Create: `src/step/mod_step_oce.F90`; Modify: `src/step/mod_model.F90`
- [ ] assemble the faithful D6 sequence with feature dispatch (`select case`) + byte-identical off-switches
- [ ] final `update_thickness_ale` commit
- [ ] **Gate:** **per-substep `max|Δ|=0` on `pi`** across the whole ocean step

#### Task M2.10: Forcing (JRA55 bulk + SW penetration)
**Files:** Create: `src/forcing/mod_forcing_bulk.F90`, `src/forcing/mod_forcing_read.F90`
- [ ] JRA55-do bulk formulae; **bilinear time-interp association order faithful** (the 2.4e6
      cancellation); g2r wind rotation on input; shortwave penetration
- [ ] (SSS/SST restoring lives in the tracer solve M2.9a, NOT here — only prepare surface flux/coeff inputs)
- [ ] **Gate:** operator-diff `max|Δ|=0` on forcing fields

#### Task M2.11: Full driver + single-rank CORE2 byte-gate
**Files:** Create: `src/drivers/fesom.F90`, `src/step/mod_model.F90` (full lifecycle), `src/io/` (minimal)
- [ ] `fesom` driver: read CORE2 mesh + **PHC IC (partition-matched)** + JRA55; `model_init/step/finalize`
- [ ] PHC IC: in-situ→potential T conversion; `extrap_nod3D` fill (accept partition-order dependence)
- [ ] **Gate:** **per-substep `max|Δ|=0` vs FESOM2 on single-rank CORE2**

#### Task M2.12: Multi-rank + production validation (MVP exit)
**Files:** Create: the local-mesh remap in `mod_mesh_read`/a new builder; Modify: tests/scripts; `docs/HANDOFF.md`
- [ ] **local-mesh remap** (prerequisite for ALL multi-rank): from the global mesh + `dist_<NP>/`
      (myList, com-structs already read), build the per-rank LOCAL mesh — remap global→local node/
      elem/edge ids, local `elem2D_nodes`/`edges`/`edge_tri`, local `nod_in_elem2D` IN LOCAL ELEMENT
      ORDER (L8), local geometry via `compute_geometry`. Replace the `npes/=1` error in
      `mod_mesh_read.F90`. Byte-gate the local geometry per rank vs a multi-rank FESOM2 geom dump.
- [ ] **M1 advection multi-rank gate (folded from M1.5):** lift the M1.1–M1.4 kernels/driver to
      FESOM2's multi-rank structure — exchange_nod/elem of tr_xy/edge_up_dn_grad/fct_LO/fct_plus_minus/
      del_ttf at the FESOM2 sites; scatter over `myDim_edge2D`/`myDim_nod2D`; then byte-gate
      `del_ttf` on pi `dist_2`/`dist_8` vs a multi-rank FESOM2 advection reference (post-exchange
      OWNED values, same partition); exchange-and-compare probe clean.
- [ ] rest-stays-at-rest; SSH gravity wave at expected speed
- [ ] **1/8/32-rank byte-match** T/S at step 200 (accept only `extrap_nod3D` diff); halo-identity clean
- [ ] **production dt=1800 (reduced M2 namelist): PRIMARY gate = per-substep `max|Δ|=0` vs the FESOM2
      reduced-config run**; SECONDARY sanity = run past the historical blow-up step, non-drifting
- [ ] **Gate:** all green ⇒ **architectural MVP**; tag `m2-mvp`; write `docs/LESSONS.md` updates

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
