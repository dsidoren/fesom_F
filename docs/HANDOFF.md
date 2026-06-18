# FESOM3 — Handoff (durable state across sessions)

Single source of truth for "where are we / what's next". Update at the end of every task.

## Where we are

- **Milestone:** M0 (Foundation) — **COMPLETE ✓** (tag `m0`). 13/13 ctest green on
  Intel dp + GNU dp; debug build (`-check all`) clean. Next: M1 (tracer advection).
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
- **Current task:** M1.1 — horizontal tracer advection (upwind + MUSCL).

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

## Build

```bash
./configure.sh --compiler intel --precision dp --clean --build   # anchor
./configure.sh --compiler gnu   --precision dp --clean --build   # portability
cd build_intel_dp && ctest --output-on-failure                   # self-tests
```
Anchor = Intel + DP + FESOM2-v2.7.3-exact flags (see docs/LESSONS.md L1). Build dirs:
`build_<compiler>_<precision>/`. Login-node runs of 1–8 ranks are fine for self-tests.

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
- For M1: add UV/Wvel/del_ttf_advhoriz dump calls to FESOM2 `oce_ale_tracer.F90` advection
  (node + NEW element shims) and rerun. For M2: use the reduced namelist (see below).

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

## Next task

M0.2 — `params/`: `mod_precision`, `mod_constants` (cite FESOM2 `oce_modules.F90` lines),
`mod_config` + `mod_param_phys` (namelist read-once), `hp_math_intrinsics`; `test_params`.

## Open notes / risks

- Byte-gates vs FESOM2 need the instrumented FESOM2 built + run with the reduced namelist +
  reference inputs. Self-tests (params/types/partit/halo/dump round-trips) run standalone now;
  FESOM2-oracle gates (M0.7 geometry, M1+ kernels) come online once a reference run is produced.
