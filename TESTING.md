# FESOM_F — Testing and Validation

This document describes how FESOM_F is tested, how it is compared against FESOM2, and what the
results are — correctness and scalability. For the architecture being tested, see
[IMPLEMENTATION.md](IMPLEMENTATION.md); for building and running, see [README.md](README.md).

> Naming, as elsewhere: the project is **FESOM_F** in prose; real identifiers (`FESOM3_*`, `fesom_*`,
> tags `m0…m10`) keep the repository's working name `fesom3`.

---

## Part I — Overview

FESOM_F is validated against one demanding standard: **`max|Δ| = 0`** — bit-for-bit identical results
to FESOM2. Because the numerical kernels are transcribed verbatim (see
[IMPLEMENTATION.md](IMPLEMENTATION.md)), the two codes should agree to the last bit, and the tests
assert exactly that. A non-zero difference is treated as a defect to be found and removed, not a
tolerance to be accepted.

Validation rests on two pillars:

1. **Automated self-tests** (`ctest`) — 21 cases over the infrastructure (types, partitioning, halo
   exchange, I/O decomposition, restart redistribution, …), run on 1, 2, and 8 ranks, on both the
   Intel and GNU toolchains.
2. **Byte-gates against the FESOM2 oracle** — an instrumented FESOM2 v2.7.3 and FESOM_F are run on
   identical input and their internal fields are compared substep by substep, demanding `max|Δ| = 0`.
   Gates exist for every subsystem, on single and multiple ranks, on a small test mesh and on the
   global CORE2 mesh.

A third property is checked throughout: **partition independence** — results, output, and restart
files do not depend on the number of MPI ranks.

**Headline results.** Milestones M0–M10 are byte-exact against FESOM2 on 1 rank and on up to 864
ranks, for both sea-ice solver variants. A full-physics run reproduces FESOM2 bit-for-bit over a
model year (1.1 million compared records), and a free-running two-year run is physically stable.
Restart resumes byte-exactly. On performance, FESOM_F is at or ahead of FESOM2 per timestep at scale.

---

## Part II — Reference

### 1. Philosophy: byte-exact validation and its floors

"Byte-exact" (`max|Δ| = 0`) means two independently compiled programs produce identical bits for
every compared value. This is achievable only because the kernels are transcribed verbatim and the
build and runtime are matched to FESOM2 (compiler flags, reduction order, floating-point mode — see
[IMPLEMENTATION.md](IMPLEMENTATION.md) §11). It is a stronger and more useful bar than "agrees to
N digits": it catches the class of subtle plumbing bugs that a tolerance would hide.

There are, however, a few honest **reproducibility floors** — places where exact equality is either
physically meaningless or inherited from FESOM2:

- **Denormals and flush-to-zero.** Bit-identity depends on the floating-point *runtime* mode, not
  just the source. A 2-year run once diverged because FESOM_F kept a denormal value FESOM2 flushed to
  zero; matching the mode (FTZ on) restored equality (LESSONS L51).
- **Reduction order.** Global sums are bit-identical only because FESOM_F mirrors FESOM2's exact MPI
  reduction; this is engineered, not assumed (LESSONS L6).
- **Redundantly-owned elements at multiple ranks.** FESOM2 assigns an element to a rank if *any* of
  its nodes is owned, so on CORE2/dist_2 about **562 elements** are "owned" by more than one rank and
  their copies can differ by ~1 ULP from order-dependent sums. This is inherent to FESOM2 (its own
  netCDF restart deduplicates them); FESOM_F's partition-independent stores deduplicate too. It is
  the one place the multi-rank bar is a tiny ε rather than exactly zero — see §7–8.

### 2. The FESOM2 oracle

The reference ("oracle") is **FESOM2 v2.7.3** (git SHA `9271ae92`), at
`/home/a/a270088/port2/fesom2`, instrumented with dump shims that write its internal fields at each
substep. It is run with `tools/run_oracle_pi.sh`, which sets up a run on the small **pi** mesh and
emits per-rank, gid-keyed dumps. The oracle is deterministic: a fresh run reproduces a committed
reference fixture bit-for-bit, which is what lets it serve as a fixed target.

### 3. The dump / gid-keyed comparison

Both codes write the same binary record (`src/infra/mod_dump.F90` and the FESOM2 shim), one record
per probed field at each substep:

```
int32 step | int32 substep_id | int32 probe_gid | int32 nlevels | char[24] field_name | float64 values[nlevels]
```

The key is **gid** — the *global* node or element id — so a value can be matched across runs
regardless of how the mesh was partitioned. `substep_id` names the point in the timestep
(`PRESSURE_BV`, `MIXING`, `SSH_SOLVE`, `TRACERS`, …). `tools/dump_diff.py` groups records by
`(step, substep, gid, field)`, reports the per-record `max|Δ|`, the first diverging substep, and a
magnitude histogram that separates real divergence from floating-point noise. Its options drive the
gates:

- `--glob` merges the per-rank dump files (multi-rank comparison keyed on gid),
- `--threshold T` sets the tolerance (default `0` — byte-exact),
- `--ignore-substep` skips a substep a given configuration does not compute,
- `--ignore-step` compares the *same physical state* reached via different per-segment step counts —
  the restart gate (straight run vs. checkpoint-and-resume).

### 4. Gate methodology

Every byte-gate follows the same three steps:

1. **Run the oracle** — instrumented FESOM2, dumping its fields.
2. **Run FESOM_F** on the *identical* input state, dumping the same fields.
3. **Compare** with `dump_diff.py`, requiring `max|Δ| = 0`.

For example, `tools/run_step_gate.sh` runs the real FESOM2 `oce_timestep_ale` and the FESOM_F
`step_oce` on identical single-rank pi state and compares every per-substep node dump (density,
pressure, `N²`, `Kv`, `ssh_rhs`, `d_eta`, `hbar`, `eta_n`, `w`, `T`, `S`, layer thickness). Because
each kernel consumes the previous kernel's live output, a clean `max|Δ| = 0` proves the whole
assembled step — data flow and dispatch — is byte-faithful, not just the kernels in isolation.

### 5. The `ctest` self-test suite

`test/CMakeLists.txt` registers **21 ctest cases** across 13 test programs, run on the Intel and GNU
builds:

| Test | Ranks | Checks |
|---|---|---|
| `test_params`, `test_types` | 1 | constants/precision parameters; derived-type definitions |
| `test_partit` | 1, 2, 8 | partition invariants (owned nodes partition exactly; elements/edges boundary-redundant) |
| `test_halo` | 1, 2, 8 | halo exchange identity; stale-halo corruption probe |
| `test_dump` | 1 | the gid-keyed dump binary round-trip |
| `test_mesh` | 1 | mesh infrastructure |
| `test_io_decomp` | 1, 2, 8 | output redistribution to writer subset, chunking, partial chunks |
| `test_io_decomp_gather` | 1, 2 | the inverse (restart-read) redistribution |
| `test_io_posix` | 1 | the `bind(C)` POSIX shims (rename/fsync/rmtree) |
| `test_vector_rotate` | 1 | vector rotation invertibility (output velocity to geographic) |
| `test_io_means` | 1 | output writer lifecycle + per-field cadence |
| `dump_diff_selftest` | — | the comparison tool itself (synthetic perturbation detection) |
| `fesom_analytic` | 1, 2 | end-to-end analytic driver runs clean |

These run in seconds on a login node:

```bash
cd build_intel_dp && ctest --output-on-failure
```

### 6. The byte-gate catalog

Beyond the self-tests, byte-gates cover each subsystem as it was built. They live in `tools/` as
`run_*_gate*.sh`; representative ones:

| Group | Example gate | Mesh / ranks | What it asserts |
|---|---|---|---|
| Kernels | `run_pressure_gate.sh`, `run_step_gate.sh`, `run_forcing_gate.sh` | pi, 1 | EOS/pressure, the whole assembled step, forcing read |
| CORE2 setup | `run_geom_gate_core2.sh`, `run_ic_gate_core2.sh` | CORE2, 1 | mesh geometry; 3-D initial conditions |
| Lifecycle | `run_lifecycle_gate_core2.sh`, `run_lifecycle_fullynative_gate_core2.sh`, `…_kpp_…`, `…_zstar_…` | CORE2, 1 | multi-step coupled runs per physics configuration |
| Multi-rank | `run_step_gate_multirank.sh`, `run_lifecycle_fullynative_gate_multirank.sh`, `run_lifecycle_jra55_gate_multirank.sh` | pi/CORE2, 2–864 | the same fields match per gid across rank counts |
| Output | `run_output_gate.sh`, `run_meshdiag_gate.sh` | pi, 1/2/8 | Zarr stores round-trip; mesh diagnostics match FESOM2 |
| Restart | `run_restart_gate_core2.sh`, `run_restart_gate_multirank.sh`, `run_restartroundtrip.sh` | CORE2, 1/2/8 | resume reproduces the straight-through trajectory |

Each gate is `max|Δ| = 0` at its stated anchor (the multi-rank element floor of §1/§7 aside). The
HANDOFF log lists the exact field/record counts per gate.

### 7. Multi-rank and partition independence

Correctness across ranks is tested in two complementary ways:

- **Per-gid byte-gates** at `dist_2`, `dist_8`, and on up to `dist_864` (CORE2). Because dumps are
  keyed on global id, the multi-rank result is compared against the single-rank (and the oracle)
  value for the same entity.
- **Partition independence**: the output gate's `zarr_diff.py --output-cmp` asserts that the Zarr
  stores written at np=1, np=2, and np=8 are byte-identical (`max|Δ| = 0`) for every variable. The
  restart path adds a cross-rank check: a checkpoint written at one rank count restores exactly at
  another (`C2 ≡ C8`).

The single exception is the ~1-ULP redundantly-owned element floor (§1): at np>1 with active forcing,
element-derived fields can differ by ≤1 ULP (e.g. temperature `2.22e-16`). The multi-rank gate
therefore asserts strict `max|Δ| = 0` at np=1 and admits a relative ε floor (`1e-12`, far above the
~`2e-16` residual and far below any real regression) at np>1. Fixing it would require per-rank
element storage and would *break* the partition-independence guarantee, so it is deliberately left as
FESOM2 leaves it.

### 8. Output and restart validation

**Output (Zarr).** `tools/zarr_diff.py` validates stores several ways: `--roundtrip` (values match a
known generator formula), `--lz4` (compressed read-back), `--meshdiag` (the mesh-diagnostics store
matches the FESOM2 oracle), `--output`/`--output-cmp` (per-variable correctness and
partition-independence), and `--frame` (vector rotation to geographic vs. native). The output gate
also sweeps the writer knobs (chunk sizes, compression, writer count) and asserts every variant is
value-identical to the default.

**Restart.** The restart gate runs the model straight through N steps, and separately runs it K
steps, checkpoints, starts a *fresh process*, reads the checkpoint, and runs the remaining N−K steps
— then compares the full state (including sea ice and the EVP stress tensor `sigma`) with
`dump_diff.py --ignore-step`. Results:

- **Same rank count:** byte-exact resume — `max|Δ| = 0` at np=1 for the entire state (FESOM_F
  serializes `sigma`, which FESOM2's own restart omits, so even the ice stress matches), and for the
  whole state at np>1 in the quiet case; with active forcing, the ≤1-ULP element floor of §1/§7
  applies.
- **Different rank count:** the state is *restored* exactly (partition-independent storage,
  `C2 ≡ C8`), after which the run continues on a physically valid — but, because element ownership
  differs, not bit-identical — trajectory.

### 9. Results

**Milestones (all byte-exact, `max|Δ| = 0`, against FESOM2 v2.7.3 at the stated anchor):**

| Milestone | Tag | Subsystem |
|---|---|---|
| M0 | `m0` | foundation: types, partitioning, halo, dump, mesh |
| M1 | `m1` | tracer advection (FCT) |
| M2 | `m2-mvp` | dynamical core + multi-rank whole-step |
| M3 | `m3` | sea ice EVP/mEVP + air–sea coupling |
| M4 | `m4` | GM/Redi |
| M5 | `m5` | KPP vertical mixing + shortwave penetration |
| M6 | `m6` | zstar vertical coordinate |
| M7 | `m7` | TKE vertical mixing |
| M8 | `m8` | multi-year production runs (CORE2, JRA55-do) |
| M9 | `m9` | Zarr output |
| M10 | `m10` | restart / checkpoint |

Each milestone is byte-exact on 1 rank and multi-rank, for both EVP variants.

**Full-physics headline run (M8).** The full default configuration (zstar + TKE + GM/Redi +
shortwave penetration + native sea ice, JRA55-do forcing from 1958) on CORE2 at **dist_864** (864
ranks) reproduces FESOM2 bit-for-bit over a **full model year — 1,138,800 compared records**,
spanning ~17,280 steps (3.4× past the day-107 mark where the flush-to-zero issue, §1/LESSONS L51,
once appeared). A free-running **two-year run (35,040 steps)** with no oracle is physically stable —
no NaN/Inf, bounded global peaks, a correct seasonal ice cycle — confirming the model is not merely
byte-faithful but well-behaved.

**Compiler / precision matrix.**

| Build | Status |
|---|---|
| Intel + double | **anchor** — byte-exact vs FESOM2; all gate results |
| GNU + double | portability — all 21 `ctest` cases pass; pure-data fields bit-identical to Intel (only rotation-derived coordinates differ at ~1e-13, expected libm differences) |
| single precision | builds (scaffolding, §IMPLEMENTATION 6); not a byte-gate target |

### 10. Scalability and performance

Performance is measured, not estimated, with a **two-point method**: run from cold start for 100 and
for 600 steps and take per-step cost as `(T₆₀₀ − T₁₀₀)/500`, which removes fixed setup cost. FESOM_F
reports per-component timings directly (`src/infra/mod_timer.F90`, in ms/step); FESOM2 reports
seconds-per-run, so a fair comparison must divide by the step count first.

The one component that was initially slower than FESOM2 was **surface forcing**: reading the
DEFLATE-compressed JRA55-do netCDF re-opened the file and read both time brackets on every rank at
every crossing — 15.74 ms/step vs FESOM2's 7.03 (2.24×, +8.7 ms/step). Switching to a persistent file
handle with a double buffer (`FESOM3_FORCING_PERSIST`, default on; `src/forcing/mod_forcing_read.F90`)
fixed it. Measured at **dist_512** (steady state, two-point):

| Component | FESOM2 | FESOM_F (before) | FESOM_F (after) |
|---|---:|---:|---:|
| forcing | 7.03 | 17.71 | **5.78** (3.07× faster than before; now below FESOM2) |
| ocean step | 30.24 | — | **28.66** (0.95×, FESOM_F faster) |
| **loop total** | 47.81 | 58.26 | **44.49** |

So after the fix FESOM_F's whole timestep at 512 ranks (44.49 ms/step) is below FESOM2's (47.81) —
its all-ranks forcing read beats FESOM2's rank-0-read-plus-broadcast on this machine. The full
investigation, including the ruled-out hypotheses, is in
`docs/plans/2026-06-28-forcing-perf-investigation.md`. (Systematic scaling/throughput work remains
future — see [IMPLEMENTATION.md](IMPLEMENTATION.md) §12.)

### 11. Reproducing the tests

On a Levante login node:

```bash
./configure.sh --compiler intel --precision dp --clean --build   # build the anchor
cd build_intel_dp && ctest --output-on-failure                   # 21 self-tests
cd .. && bash tools/run_step_gate.sh                             # one byte-gate vs FESOM2 → max|Δ| = 0
```

Multi-rank and CORE2 gates (and the dist_512/dist_864 SLURM jobs) are the `tools/run_*_gate*.sh` and
`tools/*.sbatch` scripts. The development log `docs/HANDOFF.md` records the exact command and
expected record count for each.

---

*See also:* [README.md](README.md) (build & run) · [IMPLEMENTATION.md](IMPLEMENTATION.md)
(architecture) · `docs/HANDOFF.md` (development log) · `docs/LESSONS.md` (reproducibility traps) ·
`docs/plans/`.
