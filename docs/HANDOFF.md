# FESOM3 — Handoff (durable state across sessions)

Single source of truth for "where are we / what's next". Update at the end of every task.

## Where we are

- **Milestone:** M0 (Foundation) — **COMPLETE ✓** (tag `m0`). 13/13 ctest green on
  Intel dp + GNU dp; debug build (`-check all`) clean.
- **M1 in progress.** **Geometry byte-gate CLOSED ✓ (on pi)** — FESOM3 mesh geometry is
  `max|Δ|=0` vs FESOM2 on pi (1-rank): elem_area, elem_cos, metric_factor, gradient_sca,
  edge_dxdy, edge_cross_dxdy, area/areasvol(+inv), coord_nod2D, elem2D_nodes, edges,
  edge_tri, all level arrays. Run it: `tools/run_geom_gate.sh`.
  ⚠️ **Caveat:** pi has 0/5839 CW swaps, so the `enforce_cw_orientation` vertex-reorder
  path is byte-identical-by-construction to FESOM2 `test_tri` but NOT yet empirically
  gated (soufflet=228/5700, CORE2 will have swaps). Confirm at the M2.11 CORE2 gate or on
  soufflet. See LESSONS L8.
- **M1.1 horizontal tracer advection byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs
  FESOM2 on every field: `del_ttf_advhoriz` (the gate target) AND `adv_flux_hor` for BOTH
  upwind (UPW1) and MUSCL, plus every intermediate (helem, nboundary_lay, edge_up_dn_tri,
  tr_xy, edge_up_dn_grad). Run it: `tools/run_advhor_gate.sh`. Built: `src/oce/`
  oce_tracer_grad (tr_xy), oce_muscl_adv (nboundary_lay/edge_up_dn_tri/edge_up_dn_grad),
  oce_adv_tra_hor (upw1/muscl/mfct), oce_adv_tra_flux (scatter).
- **M1.2 vertical tracer advection byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs
  FESOM2 on `del_ttf_advvert` (gate target) AND `adv_flux_ver` for BOTH upwind (UPW1) and
  QR4C, plus the vertical geometry it consumes (`zbar_3d_n`, `Z_3d_n`, `area`) and the
  prescribed `wvel`. Built `src/oce/oce_adv_tra_ver.F90` (adv_tra_ver_upw1, adv_tra_ver_qr4c).
  Same `tools/run_advhor_gate.sh` (extended). PASSED first run; debug `-check all` clean.
- **M1.3 FCT (Zalesak) limiter byte-gate CLOSED ✓ (on pi 1-rank).** `max|Δ|=0` vs FESOM2
  on the whole FCT path: `fct_LO` (low-order solution), the clipped `adv_flux_{hor,ver}_fct`,
  the limiter internals (`fct_ttf_max/min` bounds, `fct_plus/minus` factors, pre-clip
  `*_fct_ho` fluxes), the final `del_ttf_{advhoriz,advvert}_fct`, plus `hnode/hnode_new` and
  the standalone `adv_flux_hor_mfct`. Config = pi namelist (MFCT opth=0.0 / QR4C optv=1.0 /
  FCT). Built `src/oce/oce_adv_tra_fct.F90` (`oce_tra_adv_fct`); the `use_lo` branch of
  `oce_adv_tra_flux.F90` is now LIVE. Same `tools/run_advhor_gate.sh` (extended). PASSED first
  run; debug `-check all` clean. The limiter actively clipped (~41% of nodes, fct_plus/minus
  min=0; ~14% of fluxes changed), so b1/b2/b3 were genuinely exercised. See LESSONS L11.
  **Next: M1.4** (driver `do_oce_adv_tra` + dispatch + `model_step` integration).
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
  M0.7-geom ✓ (closed in M1) mesh geometry byte-matches FESOM2 on pi 1-rank.
  M1.1 ✓ horizontal tracer advection (upw1 + MUSCL): operator-diff `max|Δ|=0` on
  `del_ttf_advhoriz` + `adv_flux_hor` + all intermediates vs FESOM2 on pi 1-rank.
  M1.2 ✓ vertical tracer advection (upw1 + QR4C): operator-diff `max|Δ|=0` on
  `del_ttf_advvert` (oce_ale_tracer.F90:240) + `adv_flux_ver` + vertical geometry
  (zbar_3d_n/Z_3d_n/area) vs FESOM2 on pi 1-rank.
  M1.3 ✓ FCT (Zalesak) limiter (oce_adv_tra_fct.F90 → oce_tra_adv_fct): operator-diff
  `max|Δ|=0` on fct_LO + clipped adv_flux_{hor,ver}_fct + fct_ttf_max/min + fct_plus/minus
  + del_ttf_*_fct + hnode/hnode_new + standalone MFCT vs FESOM2 on pi 1-rank. Activated the
  `use_lo` hnode/hnode_new branch of oce_adv_tra_flux. Limiter clipped ~41% of nodes.
- **Current task:** M1.4 — tracer advection driver `do_oce_adv_tra` (dispatch on
  tra_adv_{hor,ver,lim}) + `init_tracers_AB` (AB interpolation → valuesAB, the tr_xy/
  edge_up_dn_grad fill) + wire into `model_step`. The FCT path orchestration is proven
  inline by the M1.3 gate; M1.4 lifts it into the real driver and gates the assembled step.

## Geometry byte-gate (CLOSED ✓) — the 1-rank FESOM2 oracle recipe

The first true byte-gate vs the live oracle. Reusable for ALL M1+ kernel gates.
- **FESOM2 must run 1-rank** (so per-node `area`/`nod_in_elem2D` accumulation order ==
  FESOM3 global order; multi-rank reorders → ULP diffs). METIS can't make `dist_1`, so
  it is HAND-CRAFTED in `port2/fesom2/tests/data/MESHES/pi/dist_1/` (mirrors
  `save_dist_mesh` for np=1: `rpart.out` = npes/counts/identity node→contiguous map;
  `my_list00000.out` = identity lists; `com_info00000.out` = empty com_structs with
  BLANK lines for the zero-size halo arrays, which the reader's zero-trip `read(*,*)`
  skips). FESOM2's reader accepts it (validated).
- **FESOM2 geom dump:** `src/fesom_geom_dump.F90` (NEW), wired at end of `mesh_setup`
  (oce_mesh.F90), env-gated `FESOM_GEOM_DUMP`, npes==1 only, writes full arrays and
  STOPS before forcing (1-rank forcing init hangs on login node — irrelevant, geometry
  is complete at mesh_setup). Rebuilt: `build/bin/fesom.x` + `build/lib64/libfesom.so`
  (the proven oracle `bin/fesom.x`+`lib64/` is UNTOUCHED; backups `.proven` exist).
- **FESOM3 side:** `src/infra/mod_geom_dump.F90` + driver `src/drivers/fesom_geomdump.F90`
  (same binary format) + `tools/geom_diff.py` (field-by-field max|Δ|). Rotation matches
  the pi namelist (alpha/beta/gamma=50/15/-90, cyclic 360, force_rotation).
- **Run:** `tools/run_geom_gate.sh` → PASS. Individually: `tools/run_geomdump_pi.sh`
  (FESOM2) then `fesom_geomdump` (FESOM3) then `geom_diff.py`.

## M1.1/M1.2/M1.3 advection byte-gate (CLOSED ✓) — the kernel-gate recipe (reusable for M1.4+)

Same 1-rank oracle as geometry, but the FESOM2 dump fires LATER (end of `ocean_setup`,
after `init_thickness_ale` builds `helem` + `muscl_adv_init` builds nboundary_lay/
edge_up_dn_tri), still BEFORE forcing (the 1-rank hang). The dump PRESCRIBES analytic
inputs and drives the REAL FESOM2 kernels, so it gates actual FESOM2 code:
- **FESOM2 side:** `src/fesom_advhor_dump.F90` (NEW), wired at end of `ocean_setup`
  (`oce_setup_step.F90`), env-gated `FESOM_ADVHOR_DUMP`, npes==1, STOPS after dumping.
  It sets `ttf`(nodes)/`vel`(elements) from an analytic formula of the rotated coords,
  then calls FESOM2's own `tracer_gradient_elements`→`fill_up_dn_grad`→`adv_tra_hor_upw1`/
  `_muscl`→`oce_tra_adv_flux2dtracer`. Built into `build/` (proven oracle `bin/*.proven`
  UNTOUCHED). dt=1800, num_ord=0.75 are pinned constants shared with FESOM3.
- **FESOM3 side:** `src/drivers/fesom_advhordump.F90` + `src/infra/mod_advhor_dump.F90`
  (FADVHDMP format, adds a 3-D-array writer). Same analytic prescription (byte-identical
  coords ⇒ byte-identical `ttf`/`vel`), same transcribed kernels.
- **Inputs that had to match:** `helem` (linfs/zstar agree at init since hbar=eta=0:
  `helem(nz,e)=zbar(nz)-zbar(nz+1)`), `areasvol` (geom-proven), and — crucially —
  `nod_in_elem2D` ORDERING, which the area gate already proved transitively (see L9).
- **Run:** `tools/run_advhor_gate.sh` → PASS (now gates BOTH horizontal M1.1 and vertical
  M1.2 fields). Individually: `tools/run_advhordump_pi.sh` (FESOM2) then `fesom_advhordump`
  (FESOM3) then `advhor_diff.py`.
- **M1.2 extension (DONE):** both the shim and the FESOM3 driver prescribe an analytic
  `wvel(nz,n)=1e-4·sin(2·lon)·cos(lat)·cos(0.3·nz)` (varies sign in space+depth → both
  upwind branches + QR4C interior exercised), call `adv_tra_ver_upw1`/`_qr4c` → scatter
  with a ZERO horizontal flux (so only the vertical branch fills `del_ttf_advvert`), and
  dump `wvel/zbar_3d_n/Z_3d_n/area/adv_flux_ver_{upw1,qr4c}/del_ttf_advvert_{upw1,qr4c}`.
  The FESOM2 shim uses the live `mesh%Z_3d_n/zbar_3d_n` (built by `init_ale` earlier in
  ocean_setup); the FESOM3 driver builds them from the init_ale formula (LESSONS L10).
- **M1.3 extension (DONE):** both the shim and the FESOM3 driver prescribe a second tracer
  `ttfAB` (sharp/large oscillatory, to force clipping), build `hnode/hnode_new`, and run the
  FCT branch inline — LO upw1(ttf) horiz+vert → `fct_LO`; HO MFCT(ttfAB)@opth + QR4C(ttfAB)@optv
  → antidiffusive `adv_flux_*`; `oce_tra_adv_fct` clips; `oce_tra_adv_flux2dtracer(use_lo=.TRUE.)`
  scatters. The shim drives the REAL FESOM2 `oce_tra_adv_fct` (passing a SEPARATE `aux_scratch`,
  not `edge_up_dn_grad`, so the latter survives for its dump). 15 FCT fields gated `max|Δ|=0`.
  Key facts in LESSONS L11: pi config = MFCT(opth=0)/QR4C(optv=1)/FCT; `edge_up_dn_grad`=grad(ttf)
  not grad(ttfAB); the a2 bignumber bottom-layer fill; the AUX cavity caveat.

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

M1.4 — tracer advection driver + step integration. Create `src/oce/oce_adv_tra_driver.F90`
(`do_oce_adv_tra`) transcribing FESOM2 `oce_adv_tra_driver.F90:46-419`: dispatch on
`tra_adv_hor` (UPW1/MUSCL/MFCT), `tra_adv_ver` (UPW1/QR4C/CDIFF/PPM) and `tra_adv_lim` (FCT
or not) — the exact branch the M1.3 gate already proved inline (LO build → `fct_LO` →
HO antidiffusive → `oce_tra_adv_fct` → `oce_tra_adv_flux2dtracer(use_lo)`). Add
`init_tracers_AB` (oce_tracer_mod.F90: AB interpolation `values→valuesAB`, the
`tracer_gradient_elements(values)`/`fill_up_dn_grad` fill — NB gradient from `values`, NOT
`valuesAB`, see L11) and wire `do_oce_adv_tra` into `model_step`/`adv_tracers_ale`. All
kernels (UPW1/MUSCL/MFCT/QR4C + scatter + FCT limiter + tr_xy/MUSCL setup + vertical geometry
+ hnode/hnode_new) are byte-proven (M1.1+M1.2+M1.3); M1.4 assembles them under the real driver
and gates the per-tracer `del_ttf` over a real (prescribed-velocity) step. Then M1.5 (multi-
rank, tag m1). Non-FCT branches (`do_zero_flux=.true.`) and the `tra_adv_ph/pv` per-tracer
order knobs still need a gate (the M1.3 gate fixed opth=0.0/optv=1.0); cover them at M1.4.

## Open notes / risks

- Byte-gates vs FESOM2 need the instrumented FESOM2 built + run with the reduced namelist +
  reference inputs. Self-tests (params/types/partit/halo/dump round-trips) run standalone now;
  FESOM2-oracle gates (M0.7 geometry, M1+ kernels) come online once a reference run is produced.
