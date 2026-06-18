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

Architecture brainstorm complete and plan written/reviewed; implementation not yet started.

- **Plan:** [`docs/plans/2026-06-18-fesom3-architecture.md`](docs/plans/2026-06-18-fesom3-architecture.md)
  — M0–M2 detailed into byte-gated tasks, M3–M6+ roadmap.
- **Design references** (not tracked here): `design_refs/tracer_dwarf` (Fortran structural template),
  `design_refs/fesom3-design` (FESOMx design docs); `paper/` (the porting-experience paper).
- **FESOM2 oracle** (transcribe from / byte-gate against): `/home/a/a270088/port2/fesom2/src/`.

## Next step

Begin milestone **M0** (foundation: build system, derived types, MPI/halo infrastructure, and the
validation harness) per the plan.
