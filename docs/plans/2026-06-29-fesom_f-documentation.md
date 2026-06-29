# FESOM_F User & Developer Documentation (README + IMPLEMENTATION + TESTING)

## Overview

Produce three repository-root documentation files that describe **FESOM_F** — the
clean-architecture Fortran reimplementation of FESOM2 (ocean + sea ice) currently living in
this repo under the working name `fesom3`:

1. **`README.md`** — REPLACE the existing stale one (it still says "M0 complete, next step M1").
   A concise, code-block-heavy guide to compiling on Levante and running with different meshes.
2. **`IMPLEMENTATION.md`** — NEW. The architecture and design-decisions document: what differs
   from FESOM2 and why. Short overview, then a comprehensive reference.
3. **`TESTING.md`** — NEW. The validation approach: how FESOM_F is byte-compared against FESOM2,
   the test catalog, results, and scalability. Short overview, then a comprehensive reference.

**Problem it solves:** the only current entry points are the stale README and the extremely dense
internal logs (`docs/HANDOFF.md`, `docs/LESSONS.md`). A collaborator or future-self cannot currently
(a) build & run from the README, (b) understand the architecture without reverse-engineering the
source, or (c) understand the validation story without parsing 1400 lines of milestone shorthand.

**This came out of a completed brainstorm. All design decisions are settled (see Context). Do not
re-litigate them.**

## Context (from discovery)

### Decisions already made (brainstorm)
- **Audience:** internal team / future-self, BUT keep jargon limited — define byte-gate terms on
  first use; do not use milestone-code shorthand (e.g. "M2.12c-3") as load-bearing prose references.
- **Naming:** brand the project **FESOM_F** in all prose. Document concrete identifiers with their
  REAL names (they still say `fesom3`): env-var prefix is literally `FESOM3_`; repo dir is `fesom3`;
  build dirs `build_<compiler>_<precision>`; git tags `m0…m10`; drivers `fesom_*`. Add ONE line in
  README explaining the `FESOM_F` ↔ `fesom3` legacy naming. **Do NOT rename code or env vars.**
- **Local/absolute Levante paths are OK** (internal audience): `/pool/data/AWICM/...`, the FESOM2
  oracle at `/home/a/a270088/port2/fesom2`, etc.
- **Style:** README concise; IMPLEMENTATION.md & TESTING.md = short focused OVERVIEW then
  comprehensive reference.
- **Existing docs:** replace README; KEEP `docs/HANDOFF.md`, `docs/LESSONS.md`, `docs/plans/` as
  internal logs; the three new docs cross-link to them and to each other.
- **Run config presentation:** document the `FESOM3_*` environment-variable reality (what the
  production drivers use today) + a short callout that the consolidated single-file `namelist.config`
  workflow is the intended future UX — being precise that the `mod_config` write-once *module* exists
  now, but the lifecycle drivers don't read `namelist.config` yet. POINT to where env vars are defined
  rather than only hand-listing them.

### Project state (grounded)
- `main` was fast-forwarded to commit `4991e62` (= tag `m10` + 1 doc commit) at the start of this
  work, so the restart/checkpoint code is present on `main`.
- FESOM_F v1 is **feature-complete & byte-exact vs FESOM2**: milestones **M0–M10** complete & tagged
  `m0…m10`. **Next milestone is OPEN.**
  - M0 foundation · M1 advection · M2 dynamical core · M3 sea ice EVP/mEVP · M4 GM/Redi · M5 KPP ·
    M6 zstar ALE · M7 TKE · M8 long production runs (JRA55-do) · M9 Zarr output · M10 restart/checkpoint.
- **Validation contract:** `max|Δ|=0` (byte-exact) vs an instrumented FESOM2 v2.7.3 oracle
  (SHA `9271ae92`), 1-rank AND multi-rank up to 864 ranks, both EVP variants.

### Files/components involved (writing targets)
- Create/replace: `README.md`, `IMPLEMENTATION.md`, `TESTING.md` (repo root).

### Grounding files to READ (for accuracy)
- **Build/env:** `CMakeLists.txt`, `cmake/fesom_flags.cmake`, `configure.sh`, `env.sh`,
  `env/levante.dkrz.de/shell.intel`, `env/levante.dkrz.de/shell.gnu`, `config/namelist.io`.
- **Run:** `src/drivers/fesom_lifecycle_native_mr.F90` (env-var parsing block + restart vars),
  `tools/run_lifecycle_2yr_freerun_f3_dist864.sbatch`, `tools/run_bench_f3_scaling.sbatch`,
  `tools/run_lifecycle_jra55_gate_multirank.sh`.
- **Architecture:** full `src/` tree; `src/params/mod_precision.F90`, `mod_config.F90`;
  `src/types/mod_mesh.F90`, `mod_partit.F90`, `mod_dyn.F90`, `mod_tracer.F90`, `mod_ice.F90`;
  `src/infra/mod_halo.F90`, `mod_part_bounds.F90`; `src/step/mod_step_oce.F90`;
  `src/io/mod_io_zarr.F90`, `mod_io_means.F90`, `mod_io_decomp.F90`, `mod_io_restart.F90`,
  `mod_io_posix.F90`, `mod_io_coords.F90`.
- **Testing:** `test/CMakeLists.txt`, `tools/dump_diff.py`, `tools/zarr_diff.py`,
  `tools/run_oracle_pi.sh`, `tools/run_step_gate.sh`, `tools/run_lifecycle_gate_core2.sh`,
  `tools/run_output_gate.sh`, `tools/run_meshdiag_gate.sh`, `tools/run_restart_gate_core2.sh`,
  `tools/run_restart_gate_multirank.sh`, `docs/plans/2026-06-28-forcing-perf-investigation.md`,
  `docs/plans/completed/2026-06-28-restart-checkpoint.md`, `docs/HANDOFF.md`, `docs/LESSONS.md`.

## Development Approach

- **This is a documentation task, not a code task.** The standard plan template's "unit tests" are
  replaced by **fact-verification steps**: every command, path, flag string, count, env-var name,
  and numeric result printed in a doc MUST be checked against the actual repository before the task
  is considered done. Treat a wrong path or a stale number with the same severity as a failing test.
- **Verify, do not guess.** All "VERIFICATION ITEMS" (below) must be resolved by reading source or
  running a command — never paraphrased from memory or from this plan's grounded facts (which are a
  starting point, not authority).
- Complete each document fully (write + self-verify) before moving to the next.
- Keep prose readable and low-jargon; define a term the first time it appears.
- Do NOT modify source code, env vars, tags, or build dirs. Docs only.

## Testing Strategy (= verification strategy for docs)

- **Per-task fact-checks (required):** each writing task ends with verification checklist items that
  run the doc's own commands / `grep` the claimed identifiers / `ls` the claimed paths.
- **Cross-document consistency:** the three docs must agree on the status line (M0–M10, tags
  `m0…m10`), the validation contract wording, and the FESOM_F↔fesom3 naming note.
- **No-regression for the README:** the build command, a gate-script path, and the env-var pointer
  must all resolve against the real tree.
- **No build/compile required** (writing docs), but `ctest -N` and a few `grep`/`find` commands ARE
  run as verification.

### VERIFICATION ITEMS (resolve at write-time; do NOT guess)
- [ ] **Git tag/commit status line (LOAD-BEARING — repeated in all three docs):** `git tag -l 'm*'`,
  `git rev-parse --short HEAD`, `git log --oneline -5`, `git merge-base --is-ancestor m10 HEAD`.
  Confirm `m0…m10` exist and `m10` is in `main`'s history before writing the status line. (NOTE: this
  session already fast-forwarded `main` `e36fe94`→`4991e62`, so `m10` IS an ancestor of HEAD — the
  Task-0 check just re-confirms it; do NOT trust the session-start snapshot, which predates the FF.)
- [ ] **Exact ctest count + per-test np breakdown:** `cd build_intel_dp && ctest -N` (reading the CMake
  source suggests ~21 = 18 `add_fesom_test` + 3 `add_test`, matching HANDOFF's "21/21" — but confirm
  with the live `ctest -N`; discard the earlier "~22" guess).
- [ ] **Exact module count:** `find src -name '*.F90' | wc -l` and a per-subsystem breakdown.
- [ ] **Exact compiler flag strings:** read `cmake/fesom_flags.cmake` (Intel anchor, GNU, Levante
  `-march`) — quote verbatim, do not paraphrase.
- [ ] **Exact Levante module versions:** `grep` `env/levante.dkrz.de/shell.intel` + `shell.gnu` for the
  real module strings (intel-oneapi, openmpi, gcc, netcdf-c/fortran) — README §2 must match these.
- [ ] **Exact env-var names + restart vars:** read `src/drivers/fesom_lifecycle_native_mr.F90`
  (core block, physics toggles, output, and the M10 restart vars: `RestartInPath`/`RestartOutPath`/
  `r_restart`/cadence). Capture the real names and the line numbers to cite.
- [ ] **Exact pi mesh location** a user points `FESOM3_MESH_DIR` at for a login-node test (shipped in
  repo vs referenced from `port2`/`design_refs`).
- [ ] **netCDF requirement + liblz4 optional-detection** behavior — confirm from `CMakeLists.txt`.
- [ ] **Oracle SHA / version** (`9271ae92`, v2.7.3) — confirm against `tools/run_oracle_pi.sh` / repo.

## Progress Tracking
- Mark completed items `[x]` immediately when done.
- Add newly discovered tasks with ➕ prefix; blockers with ⚠️ prefix.
- Keep this plan in sync with what was actually written.

## What Goes Where
- **Implementation Steps** (`[ ]`): writing each doc + verifying its facts within this repo.
- **Post-Completion** (no checkboxes): a human read-through for tone, and the manual `ushow`/xarray
  display smoke that only a person can eyeball.

## Implementation Steps

### Task 0: Pre-flight fact-gathering (resolve all VERIFICATION ITEMS)

**Files:**
- Create: `docs/plans/notes-fesom_f-docs-facts.md` (scratch capture; delete or fold in at the end)

- [ ] **Git state:** `git tag -l 'm*'`, `git rev-parse --short HEAD`, `git log --oneline -5`, `git merge-base --is-ancestor m10 HEAD && echo "m10 in HEAD"`; reconcile the M0–M10 / `m10` / commit status line (FF to `4991e62` already done this session — just re-confirm).
- [ ] Run `find src -name '*.F90' | wc -l` and `for d in src/*/; do echo "$d $(ls $d/*.F90 2>/dev/null | wc -l)"; done`; record total + per-subsystem module counts.
- [ ] `cd build_intel_dp && ctest -N` (build first if absent: `./configure.sh --compiler intel --precision dp --build`); record the exact test list + count + np variants (expect 21; verify).
- [ ] Read `cmake/fesom_flags.cmake`; copy the exact Intel, GNU, and Levante flag strings verbatim.
- [ ] `grep` `env/levante.dkrz.de/shell.intel` + `shell.gnu` for the exact module versions (intel-oneapi/openmpi/gcc/netcdf) for README §2.
- [ ] Read `src/drivers/fesom_lifecycle_native_mr.F90` env-parsing region; list every `FESOM3_*` var actually read (incl. restart) with its default and line number.
- [ ] Read `CMakeLists.txt`; confirm netCDF-Fortran requirement + liblz4 optional auto-detection + output dirs (`bin/`, `module/`).
- [ ] Confirm the pi mesh path a login-node run uses, and the CORE2 pool path + node/elem counts; confirm oracle version/SHA from `tools/run_oracle_pi.sh`.
- [ ] Record all findings in the scratch notes file so Tasks 1–3 cite verified values (no guessing).

### Task 1: Write `README.md` (replace stale one)

**Files:**
- Modify (replace): `README.md`
- Read: `configure.sh`, `env.sh`, `env/levante.dkrz.de/shell.{intel,gnu}`, `CMakeLists.txt`, `cmake/fesom_flags.cmake`, `config/namelist.io`, `src/drivers/fesom_lifecycle_native_mr.F90`, `tools/run_lifecycle_2yr_freerun_f3_dist864.sbatch`, `tools/run_step_gate.sh`

- [ ] **§1 What FESOM_F is** — 1 paragraph + status line (v1 feature-complete & byte-exact vs FESOM2; M0–M10, tags `m0…m10`) + the ONE-LINE FESOM_F↔fesom3 legacy-naming note + links to IMPLEMENTATION.md / TESTING.md.
- [ ] **§2 Requirements** — Levante (DKRZ); Intel oneAPI 2022.0.1 + OpenMPI 4.1.2 (anchor) OR GCC 11.2 + OpenMPI 4.1.2 (portability); netCDF-Fortran; optional liblz4.
- [ ] **§3 Build** — `./configure.sh --compiler intel --precision dp --clean --build` → `build_intel_dp/bin/`; note env.sh auto-loads modules; GNU portability variant; brief (verified) flag explanation.
- [ ] **§4 Quick check** — `ctest` on a login node + one byte-gate (`tools/run_step_gate.sh`) as the "it reproduces FESOM2" smoke.
- [ ] **§5 Running a simulation** — `FESOM3_*` env-var workflow with `fesom_lifecycle_native_mr`; copy-pasteable `srun` CORE2 example (modeled on the sbatch); curated physics-toggle table; output via `FESOM3_OUTPUT` + `namelist.io`; **POINT TO** `src/drivers/fesom_lifecycle_native_mr.F90` (+ example `tools/*.sbatch` / `tools/run_lifecycle_*.sh`) as the authoritative env-var source; short callout that the single-file `namelist.config` workflow is the intended future config (the `mod_config` module exists, but the production drivers don't read it yet).
- [ ] **§6 Different meshes** — `FESOM3_MESH_DIR` + auto `dist_<N>/`; mesh-dir contents; pi (login-node) vs CORE2 (pool path + counts); note on partitioning a new mesh.
- [ ] **§7 Output & reading results** — Zarr stores (one per variable per period), `xarray.open_zarr` one-liner, `fesom.mesh.diag.zarr`.
- [ ] **§8 Repository layout** — short source tree (src/ subsystems, tools/, test/, config/, docs/).
- [ ] **§9 Documentation map** — IMPLEMENTATION.md, TESTING.md, `docs/HANDOFF.md`, `docs/LESSONS.md`, `docs/plans/`.
- [ ] **VERIFY:** run the §3 build command form against `configure.sh`; `ls` the bin path; `grep` each env var in the driver; `ls` `tools/run_step_gate.sh` + the mesh-dir files; confirm flag text matches `cmake/fesom_flags.cmake`.

### Task 2: Write `IMPLEMENTATION.md` (new)

**Files:**
- Create: `IMPLEMENTATION.md`
- Read: `src/params/mod_precision.F90`, `mod_config.F90`; `src/types/mod_{mesh,partit,dyn,tracer,ice}.F90`; `src/infra/mod_halo.F90`, `mod_part_bounds.F90`; `src/step/mod_step_oce.F90`; `src/io/mod_io_{zarr,means,decomp,restart,posix,coords}.F90`; `docs/LESSONS.md` (for L29/L51 refs)

- [ ] **Part I — Overview** — guiding principle ("first reproduce exactly, then optimize"); kernels transcribed line-for-line (byte-identical) vs plumbing redesigned; a FESOM2↔FESOM_F summary table; v1 scope & deferred.
- [ ] **§1 Goals & byte-exact validation contract.**
- [ ] **§2 Source organization** — 9 subsystems (params/types/infra/mesh/oce/ice/forcing/io/step) + drivers; module-count table (verified counts from Task 0).
- [ ] **§3 Data model** — derived types & dependency injection (`t_mesh`/`t_partit`/`t_dyn`/`t_tracer`/`t_ice`; no globals; prognostic-vs-work split; config-in-types). Quote one verified type/signature.
- [ ] **§4 Serial↔parallel** — optional-`partit` + `owned_bounds` pattern (one kernel byte-identical at 1 rank, correct at N); quote a real signature from `mod_step_oce.F90` or a kernel.
- [ ] **§5 Generic halo / MPI layer** — `exchange_nod`/`exchange_elem`/`allreduce`; `com_struct`; reduction determinism (verify names in `mod_halo.F90`).
- [ ] **§6 Precision scaffolding** — two-tier WP/MP; compile-flag switch; anchor WP=MP=8; mixed-precision future (quote `mod_precision.F90`).
- [ ] **§7 Mesh-arity generalization** — `MAX_NV=4`; triangle anchor, quad-ready (verify in `mod_mesh.F90`).
- [ ] **§8 Configuration** — the write-once `mod_config` module (exists today; reads `namelist.config` once then read-only) vs the honest reality that the production lifecycle drivers currently configure via `FESOM3_*` env vars + in-code defaults; the future single-file intent; point to where env vars are read.
- [ ] **§9 Timestep & lifecycle** — `step_oce` sequence; native-flux coupled driver.
- [ ] **§10 I/O & checkpointing** — hand-rolled Zarr v2; distributed chunk writers, NO rank-0 gather; partition-independent stores; `mod_io_decomp` + `decomp_gather`; restart `mod_io_restart` (atomic tmp→rename, `restart.latest`, keep-N prune, full state incl. ice + EVP `sigma`), `mod_io_posix` bind(C) shims; contrast FESOM2 rank-0-gather netCDF.
- [ ] **§11 Correctness decisions enabling byte-identity** — flush-to-zero (FTZ) match (L51); CG preconditioner `!DIR$ NOVECTOR` (L29); reduction order; the ~1-ULP redundant-element floor — cross-link LESSONS (verify L-numbers exist).
- [ ] **§12 Scope, limitations & roadmap** — done M0–M10; deferred (scientific/multi-decade validation, mixed precision WP→f32, performance/SYPD, beyond-v1 physics: IDEMIX/backscatter/BGC/cavity/icebergs/coupling).
- [ ] **VERIFY:** `grep` every type name, module name, `MAX_NV`, WP/MP definition, and L-number against the source; confirm the module-count table matches Task 0; confirm the §10 restart behaviors (atomic `tmp→rename`, `restart.latest`, keep-N prune, ice + EVP `sigma` in the field set) by reading `src/io/mod_io_restart.F90` rather than paraphrasing.

### Task 3: Write `TESTING.md` (new)

**Files:**
- Create: `TESTING.md`
- Read: `test/CMakeLists.txt`, `tools/dump_diff.py`, `tools/zarr_diff.py`, `tools/run_oracle_pi.sh`, `tools/run_step_gate.sh`, `tools/run_lifecycle_gate_core2.sh`, `tools/run_output_gate.sh`, `tools/run_meshdiag_gate.sh`, `tools/run_restart_gate_core2.sh`, `tools/run_restart_gate_multirank.sh`, `docs/plans/2026-06-28-forcing-perf-investigation.md`, `docs/plans/completed/2026-06-28-restart-checkpoint.md`, `docs/HANDOFF.md`

- [ ] **Part I — Overview** — central claim (`max|Δ|=0` byte-exact reproduction of FESOM2); two pillars (automated self-tests + byte-gates vs FESOM2 oracle); partition-independence; headline results.
- [ ] **§1 Philosophy** — byte-exact validation & reproducibility floors (what `max|Δ|=0` means; denormals/FTZ; reduction order; redundant-element ~1-ULP) — cross-link bit-identity-reality / LESSONS.
- [ ] **§2 The FESOM2 oracle** — v2.7.3, SHA `9271ae92` (verified), instrumented dump shims; `tools/run_oracle_pi.sh`; determinism check.
- [ ] **§3 Dump / gid-keyed comparison** — `mod_dump` binary format; substep enum; `tools/dump_diff.py`; histogram/threshold (read the tool to get the real record layout).
- [ ] **§4 Gate methodology** — 3-step recipe (oracle run → FESOM_F dump → diff) + one worked example from a real gate script.
- [ ] **§5 The ctest self-test suite** — enumerate from the verified `ctest -N` list; np 1/2/8; Intel+GNU; **use the verified count**.
- [ ] **§6 Byte-gate catalog** — table by subsystem (kernel / lifecycle / multi-rank / output / restart) — mesh / ranks / records / status (from `ls tools/run_*_gate*.sh` + HANDOFF).
- [ ] **§7 Multi-rank & partition independence** — `dist_1…864`; `zarr_diff.py --output-cmp`; restart cross-np `C2≡C8`.
- [ ] **§8 Output & restart validation** — `zarr_diff.py` modes; meshdiag vs FESOM2; restart split-vs-straight round-trip.
- [ ] **§9 Results** — milestone byte-exact table M0–M10; 2-year `dist_864` headline run (35040 steps / 1.1M records `max|Δ|=0` + FTZ fix); compiler/precision matrix.
- [ ] **§10 Scalability & performance** — two-point bench method (NSTEPS 100/600); forcing-perf 17.7→5.78 ms/step = 3.07× (F3 now beats F2); `dist_512` LOOP 44.49 ms/step ≈ 0.93× F2, ocean 28.66 ≈ 0.95× F2; "measure, don't guess"; ⚠️ note F2 s/run vs F3 ms/step conversion. Pull exact numbers from the forcing-perf doc.
- [ ] **§11 Reproducing the tests** — `ctest` + a gate on a login node.
- [ ] **VERIFY:** confirm ctest count/names against Task 0; `ls` every gate script named; confirm perf numbers + oracle SHA + restart results against the source docs; pin the 2-year headline numbers (steps `35040`, ~`1.1M` records) and the restart np=1 / np>1 ≤1-ULP / cross-np `C2≡C8` results to specific lines in `docs/HANDOFF.md` or `docs/plans/completed/2026-06-28-restart-checkpoint.md` (cite, don't paraphrase).

### Task 4: Cross-link, consistency & no-regression pass

**Files:**
- Modify: `README.md`, `IMPLEMENTATION.md`, `TESTING.md`

- [ ] Add reciprocal cross-links among the three docs and to `docs/HANDOFF.md` / `docs/LESSONS.md` / `docs/plans/`.
- [ ] Make the status line identical across all three (M0–M10, tags `m0…m10`, "next milestone OPEN").
- [ ] Ensure the FESOM_F↔fesom3 naming note appears once (README) and the other docs use FESOM_F consistently in prose while keeping real identifiers (`FESOM3_*`, `fesom_*`, tags) verbatim.
- [ ] Confirm no user-specific path is presented as portable (they're fine, but framed as "on Levante").
- [ ] **No-regression spot-check:** re-run the README build command form check, `ls` one gate script, and `grep` one env var — all must resolve.
- [ ] Delete or fold the Task 0 scratch notes file.

### Task 5: Verify acceptance criteria
- [ ] All three files exist at repo root; README replaced (no "M0 complete / next step M1" text remains).
- [ ] Every command, path, flag, count, env-var, and number in the docs was verified against the tree (no un-checked claim).
- [ ] Overview requirements satisfied: README builds+runs+mesh-setup; IMPLEMENTATION covers the 12 reference sections; TESTING covers the 11 reference sections.
- [ ] Cross-links resolve; status line consistent; naming consistent.

### Task 6: [Final] Commit & housekeeping
- [ ] (If user chose to commit) commit on `main` with a clear message.
- [ ] Move this plan to `docs/plans/completed/` (`mkdir -p docs/plans/completed`).
- [ ] Optionally update `MEMORY.md` pointer noting the docs now exist (user memory).

## Technical Details

- **Output format:** GitHub-flavored Markdown. Code blocks for every command. Tables for the
  FESOM2↔FESOM_F comparison, physics toggles, module counts, ctest list, gate catalog, results,
  and perf numbers.
- **Naming rule (mechanical):** prose noun = "FESOM_F"; verbatim identifiers stay as-is
  (`FESOM3_MESH_DIR`, `fesom_lifecycle_native_mr`, `build_intel_dp`, tags `m0…m10`,
  repo dir `fesom3`). One reconciling sentence in README §1.
- **Authority order when facts conflict:** source code > `cmake/*.cmake` / `CMakeLists.txt` /
  scripts > `docs/plans/*` (incl. forcing-perf & restart-checkpoint) > `docs/HANDOFF.md` >
  this plan's grounded facts > memory. Always prefer the more authoritative source.
- **Do not** modify code, env vars, tags, or build dirs.

## Post-Completion
*Items requiring manual intervention — no checkboxes, informational only*

**Manual verification:**
- A human read-through of all three docs for tone (low-jargon, academic-but-readable) and flow.
- The manual `ushow` / `xarray.open_zarr` display smoke for an output store (a person eyeballing a
  plot) — the automated `zarr_diff.py` round-trip is the proxy, but visual confirmation is human-only.
- Optionally render the Markdown (GitHub / VS Code preview) to confirm tables and code blocks lay out.
