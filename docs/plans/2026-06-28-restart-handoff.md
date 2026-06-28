# Restart / checkpoint — milestone handoff (brainstorm-prep)

**Status:** NEXT milestone after M9 (Zarr output, tagged `m9` 2026-06-28). The M9 brainstorm deliberately
split restart OUT (output-first). This doc is the **brainstorm-prep**, not a plan: it grounds the milestone in
the FESOM2 oracle + the current FESOM3 state and surfaces the **design forks for the user to own** — same pattern
as M9 (architecture forks get surfaced and decided WITH the user before any code; see [[project-m9-brainstorm-first]]
and the memory `feedback`/`stay-close-to-fortran`). **Next session: run the brainstorm on §"Design forks", THEN
write the plan, THEN implement.** Do NOT pre-decide the forks.

---

## 1. What restart must do

- **Write** a checkpoint (end-of-run and/or periodic) holding the full prognostic state.
- **Read** it at startup when `r_restart` is true (cold start = no restart) and resume **byte-exactly**.
- The `.clock` file already drives cold-vs-restart detection (`mod_clock`); the missing half is the `.clock`
  **write** + the state file write/read.

## 2. State inventory — oracle → FESOM3 (grounded)

FESOM2 `ini_ocean_io` / `ini_ice_io` (`/home/a/a270088/port2/fesom2/src/io_restart.F90:109,236`) define the
canonical restart variable set. Mapping to FESOM3 (✓ = field exists in FESOM3 today):

**Ocean (dynamics + ALE):**
| oracle | FESOM3 | exists | notes |
|---|---|---|---|
| ssh = `eta_n` | `dyn%eta_n` | ✓ | mod_dyn.F90 |
| ssh_rhs_old | `dyn%ssh_rhs_old` | ✓ | |
| hbar / hbar_old | `mesh%hbar` / `mesh%hbar_old` | ✓ | ALE surface elevation (+ lagged) |
| hnode | `mesh%hnode` | ✓ | ALE layer thickness (prognostic) |
| u/v | `dyn%uv(1:2,:,:)` | ✓ | element velocity (nl-1, elem2D) |
| urhs_AB/vrhs_AB (+AB3) | `dyn%uv_rhsAB(AB_order-1,2,:,:)` | ✓ | Adams-Bashforth memory (default AB2) |
| w / w_expl / w_impl | `dyn%w` / `dyn%w_e` / `dyn%w_i` | ✓ | vertical velocity (likely recompute — see F-F) |
| tke (optional) | `dyn%work%tke` | ✓ | **NOT serialized yet** (only when `mix_scheme==5`) |
| iwe / uke / uke_rhs (optional) | — | ✗ | out of FESOM3 v1 scope (IDEMIX/backscatter) |

**Tracers** (`tracers%data(j)`, j=1..num_tracers; T,S + any passive):
| oracle | FESOM3 | exists |
|---|---|---|
| values | `tracers%data(j)%values` | ✓ |
| valuesAB | `tracers%data(j)%valuesAB` | ✓ |
| valuesold (M1/M2) | `tracers%data(j)%valuesold` | ✓ |

**Ice** (`ice%data(1:3)` = a_ice/m_ice/m_snow):
| oracle | FESOM3 | exists | notes |
|---|---|---|---|
| area/hice/hsnow | `ice%data(1:3)%values` | ✓ | + the FCT work arrays values_old/_rhs/etc. (recompute?) |
| uice/vice | `ice%uice` / `ice%vice` | ✓ | |
| (sigma11/12/22) | `ice%work%sigma11/12/22` | ✓ | **EVP stress — NOT in FESOM2 restart; see fork F-C** |
| ice_albedo/ice_temp (optional) | — | ✗ | OIFS/icepack coupling, out of scope |

**Clock:** `mod_clock` globals (year/day/time old+new). `r_restart` read works; the **write** is deferred.

## 3. Already in place — do NOT rebuild

FESOM3 already has unformatted-IO serialization primitives (built for the byte-gate dump infra, `src/infra/mod_dump.F90`):
- `write_t_dyn`/`read_t_dyn` (`src/types/mod_dyn.F90:176,204`) + `write_t_dyn_work`/`read_t_dyn_work` (`:156,166`)
  — already cover eta_n/uv/uv_rhsAB/w/w_e/w_i/ssh_rhs_old (+ work). **tke is NOT in `write_t_dyn_work` (M8-deferred).**
- `write_t_mesh`/`read_t_mesh` (`src/types/mod_mesh.F90:111,177`) — hnode/hnode_new/hbar/hbar_old.
- `write_t_tracer_data`/`read_t_tracer_data` (`src/types/mod_tracer.F90:73,89`) — values/valuesAB/valuesold.
- **No `t_ice` serialization yet** (`src/types/mod_ice.F90`: data(1:3) `:25`, uice/vice `:97,98`, sigma `:37`).
- `r_restart` (`src/infra/mod_clock.F90:36`); `clock_finish`/`clock_newyear` placeholders deferred (`:12-14`).
- Lifecycle already reads `FESOM3_RESTART_IN` → `RestartInPath` (`src/drivers/fesom_lifecycle_native_mr.F90:559-568`).

So the milestone is mostly: a **format/backend decision**, **t_ice serialization**, **tke**, the **clock write**,
**lifecycle wiring** (write-at-end + read-at-start), and the **reproducibility gate**.

## 4. Design forks — BRAINSTORM THESE WITH THE USER

- **F-A — file format / backend (the big one).** FESOM2 ships THREE: `raw` (per-rank Fortran core dump,
  partition-DEPENDENT `np<N>/` path), `bin` (binary), `nc` (netCDF, partition-INDEPENDENT, **gather-based**).
  FESOM3 options: **(A) reuse the M9 Zarr writer** (`mod_io_zarr`/`mod_io_decomp`) — canonical-order, distributed
  chunk-writers, NO rank-0 gather, xarray-inspectable, consistent with the M9 decision; **(B) port FESOM2 netCDF
  restart** — byte-matches FESOM2 restart files but reintroduces the gather the user rejected for output;
  **(C) raw per-rank Fortran dump** — REUSES the existing `write_t_*` primitives, simplest + byte-exact, but
  partition-dependent (FESOM2 itself blesses this for `raw`). Note the M9 output fork landed on (A); restart may
  reasonably differ because a restart is transient/internal, not an analysis product.
- **F-B — partition portability.** Must `np=2`-write → `np=8`-read work? Only the canonical formats (A/netCDF)
  allow it; raw/bin do not (and FESOM2's raw doesn't either). Decides whether F-A can be (C).
- **F-C — EVP `sigma` (gate-shaping).** `ice%work%sigma11/12/22` is elastic memory that **carries across steps**
  (FESOM3 zeroes it only at cold-start, `mod_ice_setup.F90:147`; the EVP solve updates it cumulatively,
  `mod_ice_dyn.F90:187`). **FESOM2's `ini_ice_io` does NOT serialize it** → FESOM2 restart accepts a small ice-stress
  discontinuity. Fork: **match FESOM2 (omit sigma)** — then the gate can NOT be `max|Δ|=0` on ice straight after a
  restart; **or serialize sigma** (beyond FESOM2) for a TRUE byte-exact resume. This choice defines what the gate
  can claim.
- **F-D — cadence.** End-of-run only? Periodic (`restart_length`/unit like FESOM2)? Yearly-split filenames
  (`fesom.<YYYY>.oce.restart`)? Overwrite vs keep-N.
- **F-E — clock write.** Port `clock_finish`/`clock_newyear` from FESOM2 `gen_clock` into `mod_clock` (the `.clock`
  two-line write). Small + mechanical; already stubbed.
- **F-F — prognostic vs recompute.** Which "state" is truly needed vs re-derivable at startup: `w`/`w_e`/`w_i`
  (recomputed by `vert_vel_ale`?), `hnode_new` (work buffer), the ice FCT work arrays, pressure/density (recomputed).
  FESOM2's list is the safe superset; trimming is an optimization to verify against the gate, not guess.

## 5. The gate (restart-reproducibility)

Primary (self-consistency, no FESOM2 oracle needed — FESOM3 state is already byte-exact vs FESOM2 thru M9):
**run N steps straight-through** vs **run K steps → write restart → (fresh process) read restart → run N−K steps**;
compare the full final state → `max|Δ|=0`. Reuse the `mod_dump` gid-keyed dump + `dump_diff.py --glob` machinery.
Secondary (only if F-A/F-B pick a canonical format): **cross-partition** — write at `np=2`, read at `np=8`, same
straight-through compare. The F-C decision sets whether ice is in the `max|Δ|=0` claim or carries a documented caveat.

## 6. Pointers

- **Oracle:** `/home/a/a270088/port2/fesom2/src/io_restart.F90` — public API `:30`
  (`read_initial_conditions`/`write_initial_conditions`/`finalize_restart`); `ini_ocean_io:109`, `ini_ice_io:236`;
  backend select `read_initial_conditions:306` (raw `:360` / bin `:366` / netCDF `:384`); path builders `:58-103`
  (note the `np<N>` partition-dependent dirs for raw/bin). Clock: `gen_clock.F90` (`clock_finish`/`clock_newyear`).
- **FESOM3:** types `src/types/mod_{dyn,mesh,tracer,ice}.F90` (anchors in §3); `src/infra/mod_clock.F90`,
  `src/infra/mod_dump.F90`; lifecycle `src/drivers/fesom_lifecycle_native_mr.F90` (RestartInPath `:559-568`).
- **Memory:** [[project-fesom3-implementation-state]] (build state, M0–M9), [[project-m9-brainstorm-first]]
  (the surface-and-own-the-forks pattern + the no-gather/partition-independent stance), [[stay-close-to-fortran]]
  (read the oracle, gate `max|Δ|=0`), [[reference-fesom2-sequence-facts]].
- **HANDOFF:** `docs/HANDOFF.md` "Next task" points here; the M9 plan (DONE) is `docs/plans/completed/2026-06-28-m9-zarr-output.md`.
