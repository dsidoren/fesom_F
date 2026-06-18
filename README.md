# FESOM3

A clean-architecture **Fortran** reimplementation of the FESOM2 ocean + sea-ice model.

**Goal (v1):** byte-exactly reproduce a specific FESOM2 configuration (double precision, triangle
meshes, same Intel toolchain, fixed MPI ranks). The "clean architecture" lives **only in the
plumbing** — derived types, dependency injection, a generic MPI halo layer, write-once config, and
precision/mesh-arity *scaffolding* — while the numerical **kernels are transcribed line-for-line
from FESOM2**. Guiding rule: *"first reproduce exactly, then optimize."*

Later phases (foundations laid now, **not** exercised in v1): mixed precision, quad/mixed-polygon
meshes, performance, and ports to other languages.

Target physics: nonlinear EOS, hydrostatic PGF, AB2, CG SSH solve, MFCT advection, biharmonic
viscosity, **PP / KPP / TKE** vertical mixing, **GM/Redi**, **EVP + mEVP** sea ice, **linfs + zstar**.

## Status

**M0 (Foundation) complete** — buildable skeleton + validation harness + mesh/partition/halo
infrastructure, zero physics. 13/13 self-tests green on Intel and GNU (double precision).

- **M0.1** CMake build (anchor flags from FESOM2 v2.7.3) · **M0.2** params (precision/constants/
  config) · **M0.3** types (mesh/partit/dyn/tracer/ice) · **M0.4** partitioning + 1-rank synthesis ·
  **M0.5** generic halo exchange · **M0.6** gid-keyed dump harness + `dump_diff.py` · **M0.7** mesh
  geometry + CW + analytic generator · **M0.8** analytic driver + lifecycle.
- **Plan:** [`docs/plans/2026-06-18-fesom3-architecture.md`](docs/plans/2026-06-18-fesom3-architecture.md).
  **State/next:** [`docs/HANDOFF.md`](docs/HANDOFF.md). **Gotchas:** [`docs/LESSONS.md`](docs/LESSONS.md).
- **Design references** (not tracked): `design_refs/tracer_dwarf`, `design_refs/fesom3-design`.
- **FESOM2 oracle** (transcribe from / byte-gate against): `/home/a/a270088/port2/fesom2/src/`.

## Build

```bash
./configure.sh --compiler intel --precision dp --clean --build   # anchor
cd build_intel_dp && ctest --output-on-failure
```

## Next step

Milestone **M1** — byte-exact tracer advection (FCT) on a prescribed velocity, gated `max|Δ|=0`
vs instrumented FESOM2. First task requiring the FESOM2 reference oracle (also closes the deferred
M0.7 geometry byte-gate).
