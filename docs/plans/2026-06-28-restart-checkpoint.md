# Restart / Checkpoint — partition-independent Zarr checkpoints (byte-exact resume)

> ✅/⬜ Tick the `- [ ]` boxes (and add a ✅ to a Task header) when a sub-task is **byte-gated DONE**, per
> the project discipline (memory `feedback-tick-plan-checkboxes`). Read the actual FESOM2 `.F90` oracle, not
> summaries, before porting (memory `stay-close-to-fortran`). Gate everything at `max|Δ|=0`.

## Overview

Add restart/checkpoint to FESOM3: **write** a checkpoint holding the full prognostic state (periodically and at
end-of-run), **read** it at startup when `r_restart` is true, and **resume byte-exactly** — `max|Δ|=0` across the
whole state, **sea ice included**. Checkpoints are **partition-independent** (the MPI rank count may change between
segments) and **`ushow`/xarray-viewable** exactly like M9 output.

> **Claim scope (per plan-review).** Bit-identical resume (`max|Δ|=0` on the *continued trajectory*) is a **same-np**
> guarantee (Task 6.1). Resuming at a **different** np restores the state exactly — partition-independent storage,
> `C2 ≡ C8` (Task 6.2) — but then continues on a physically-valid, *not* bit-identical trajectory, because FESOM's
> global reductions (SSH-CG `allreduce_sum`) reduce in comm-size-dependent order. This is inherent to FESOM (FESOM2
> too), not a defect; it is why the cross-np gate is a restore round-trip, not an evolution compare.

This is the milestone after M9 (Zarr output, tagged `m9`). It reuses the M9 Zarr stack (`mod_io_zarr` +
`mod_io_decomp`) and the byte-gate tooling (`mod_dump` + `dump_diff.py`). The brainstorm-prep is
`docs/plans/2026-06-28-restart-handoff.md`; the design was decided in the 2026-06-28 brainstorm.

### Why this design (forks decided in brainstorm — do NOT revisit)

- **F-A / F-B — canonical, partition-independent.** A new `src/io/mod_io_restart.F90` on the M9 primitives
  (`mod_io_zarr` chunk I/O + `mod_io_decomp` canonical redistribute), **not** `mod_io_means` (its growing time-dim /
  per-period stores / cadence / accumulation are the wrong tool). **np=2-write must read at np=8** — the operational
  requirement that ruled out a raw per-rank dump. Consistent with M9's no-gather stance.
- **F-C — include EVP `sigma`.** FESOM2's `ini_ice_io` does **not** serialize `sigma11/12/22` (verified,
  `io_restart.F90:250-254`); it accepts a small ice-stress discontinuity on restart. We **serialize sigma** so the
  gate is a clean `max|Δ|=0` on ice too.
- **F-D — periodic + end-of-run, year-stamped folders.** One **global** cadence (namelist `restart_length` /
  `restart_length_unit`, default `1 'y'`) via `is_due` semantics, plus always at end-of-run.
- **F-F — mirror FESOM2's verified field list (safe superset).** FESOM2 **does** serialize `w`/`w_expl`/`w_impl`
  (`io_restart.F90:228-230`), so we do too. Trimming "recomputables" is a post-gate optimization, not a v1 guess.
- **Registration, not binary dumps.** Because we chose the Zarr path, `mod_io_restart` registers **live array
  pointers** (eta_n, uv, tracers, ice arrays, sigma) exactly like M9's `register_output_var`, and reads back **into**
  those arrays. The handoff's `write_t_dyn`/`write_t_ice` unformatted primitives are the *raw-dump* mechanism (option
  C) — **not used here**; **no new `t_ice` serialization is required**.

## Context (from discovery)

### Files / components involved (FESOM3 — the port)
- **Create:** `src/io/mod_io_restart.F90` (the checkpoint writer/reader); a shared coord helper module
  `src/io/mod_io_coords.F90` (extracted from `mod_io_means`); restart gate scripts under `tools/`.
- **Modify:** `src/io/mod_io_decomp.F90` (add `decomp_gather`, the inverse redistribution); `src/io/mod_io_means.F90`
  (call the shared coord helper); `src/infra/mod_clock.F90` (`clock_finish`/`clock_newyear`); `src/params/mod_config.F90`
  (`RestartOutPath`); `src/drivers/fesom_lifecycle_native_mr.F90` (register / read / write hooks); `test/CMakeLists.txt`
  + root `CMakeLists.txt` (new ctests).
- **Reused as-is:** `mod_io_zarr` (`zarr_create_store`, `zarr_define_array`, `zarr_write_chunk`, `zarr_read_chunk`,
  `zarr_consolidate`); `mod_io_decomp` (`decomp_init_entity`, `decomp_redistribute`, `recv_target` inverse map);
  `src/infra/mod_halo.F90` (`exchange_nod`, `exchange_elem`); `mod_dump`.

### Build mechanics (verified)
- `./configure.sh --compiler intel --precision dp --clean --build` (anchor) + `--compiler gnu`.
- ctest: `cd build_intel_dp && ctest --output-on-failure`.
- Multi-rank: requires `env.sh` (KNEM `single_copy_mechanism=none`; memory `project-levante-mpi-knem-gotcha`).
- Python: `/work/ab0995/a270088/mambaforge/bin/python3`.

### Reference oracle (FESOM2 v2.7.3, `/home/a/a270088/port2/fesom2/src/`)
- `io_restart.F90`: `ini_ocean_io:109-231` (ocean field list, incl. `w`/`w_expl`/`w_impl` `:228-230`);
  `ini_ice_io:236-269` (**no sigma**); `write_initial_conditions:417`; `write_netcdf_restarts:584` (record append
  `rec_count()+1` along an unlimited dim `:629-632` — the time-dim design we are **replacing** with immutable folders);
  read path `:895-910` (reads the **last** record + the `int(ctime)`-matches-clock safety check); `is_due:937`.
- `gen_modules_clock.F90`: `clock_finish:170-200` (the `.clock` write), `clock_newyear:204-213`.

### Tooling
- `mod_dump` (`src/infra/mod_dump.F90`): gid-keyed per-rank dump — `FESOM_DUMP_ALL` (every owned node/elem),
  `FESOM_DUMP_FILE`/`FESOM_DUMP_FILE_ELEM`, `FESOM_DUMP_MINSTEP` (step window).
- `tools/dump_diff.py --glob` (the `max|Δ|=0` comparator), `tools/zarr_diff.py`, `tools/run_clocktest.sh`,
  `tools/run_output_gate.sh`, `tools/run_meshdiag_gate.sh`, `tools/run_zarrsmoke.sh`.
- `ushow` = `/home/a/a270088/toolbox/bin/ushow` (C binary; GUI, no headless). Contract in Technical Details.

## Development Approach
- **Testing == byte-gates** (this project's established discipline, not unit-test-framework TDD). Every task ends
  with a gate that must be GREEN before the next task starts.
- Each gate is the *appropriate* check: a **ctest round-trip** for pure routines (`decomp_gather`), a **byte-gate**
  (`dump_diff.py --glob` → `max|Δ|=0`) for integration, an **xarray/ushow smoke** for store readability, and a
  **no-regression** sweep (`ctest` + M9 output gates + production MR lifecycle gate stay `max|Δ|=0`).
- Small focused changes; read the `.F90` oracle before porting; keep this plan in sync (➕ new tasks, ⚠️ blockers).

## Testing Strategy
- **Primary (self-consistency, no FESOM2 oracle needed — F3 is already byte-exact vs F2 through M9):**
  straight-through N steps vs split (K → checkpoint → fresh process → read → N−K), `FESOM_DUMP_ALL` both, compare
  with `dump_diff.py --glob` → `max|Δ|=0` across the whole state **including ice + sigma**.
- **Secondary (cross-partition, proves F-B):** **restore round-trip** — write `C2` at np=2; a fresh np=8 process
  reads `C2` and immediately re-writes `C8` (zero steps); `C2 ≡ C8` canonically (`max|Δ|=0`). Proves partition-
  independent storage/restore, *not* cross-np evolution byte-identity (physically impossible; see Overview claim scope).
- **Smoke:** `zarr_diff.py`/xarray opens each per-field store (dims, finite `lon`/`lat`, `_ARRAY_DIMENSIONS`); manual
  `ushow` display (+ the glob-animate-across-checkpoints check).
- **No-regression:** `ctest` GREEN both compilers; `run_output_gate.sh` + `run_meshdiag_gate.sh` + `run_zarrsmoke.sh`
  byte-identical (the shared coord helper must not change M9 stores); production MR lifecycle byte-gate stays `max|Δ|=0`.

## Progress Tracking
- Mark `[x]` immediately when a sub-task is byte-gated DONE; add `✅` to the Task header.
- `➕` prefix for newly discovered tasks; `⚠️` for blockers. Update scope here if implementation deviates.

## What Goes Where
- **Implementation Steps** (`[ ]`): code, gates, scripts achievable in this repo.
- **Post-Completion** (no checkboxes): manual `ushow` display, performance tuning at scale.

## Implementation Steps

### Stage 0 — Clock write (F-E) + `RestartOutPath`
*Smallest, independent piece; unblocks the `.clock` write that drives restart detection.*

#### Task 0.1 ✅: Port `clock_finish`/`clock_newyear`, add `RestartOutPath`, gate via `run_clocktest`

**Files:**
- Modify: `src/infra/mod_clock.F90`
- Modify: `src/params/mod_config.F90` (no change needed — `RestartOutPath` already present at `:30`/`:33`)
- Modify: `tools/run_clocktest.sh`
- Modify: `src/drivers/fesom_clocktest.F90` (➕ the round-trip assertions must live in the Fortran driver that
      `run_clocktest.sh` runs — `clock_finish`/`clock_init` are module procedures, not shell-callable)

- [x] `RestartOutPath` **already exists** (`mod_config.F90:30`, in the `/paths/` namelist `:33`) — just add it to
      `mod_clock`'s `use mod_config` list (`mod_clock.F90:17` currently imports only `RestartInPath`); the lifecycle
      sets it from `FESOM3_RESTART=<dir>` in Stage 5
- [x] port `clock_finish` into `mod_clock` **verbatim** from `gen_modules_clock.F90:170-200`: write
      `RestartOutPath//runid//'.clock'`, two lines (`timeold dayold yearold` / `dum_timenew dum_daynew dum_yearnew`)
      with the year-rollover normalization (`daynew==ndpyr .and. timenew==86400` → `0.0 / 1 / yearold+1`)
- [x] port `clock_newyear` (`:204-213`, in-memory rollover used for folder naming)
- [x] drop the "OMITTED / deferred" note at `mod_clock.F90:14-15`
- [x] extend `tools/run_clocktest.sh`: `clock_finish` write → `clock_init` read round-trips **exactly** (both lines;
      explicitly exercise the year-rollover branch `daynew==ndpyr .and. timenew==86400`); equal lines ⇒ cold
      (`r_restart=.false.`), differing ⇒ restart. (List-directed `fmt=*` real I/O round-trips here only because
      `timenew` is a clean multiple of `dt` — note the dependence.)
- [x] **GATE:** `run_clocktest.sh` GREEN (write→read round-trip exact; `r_restart` detection both ways) — PASS
      np=1 **and** np=2 (Intel dp, worktree build); cases A/B/C all `max|Δ|=0`, year-rollover branch fired (C)

### Stage 1 — `decomp_gather` (the inverse redistribution)
*The one genuinely new MPI routine; it is the literal transpose of the proven `decomp_redistribute`.*

#### Task 1.1 ✅: `decomp_gather_2d_r`/`3d_r` + transpose round-trip ctest

**Files:**
- Modify: `src/io/mod_io_decomp.F90`
- Create: `test/test_io_decomp_gather.F90`
- Modify: `test/CMakeLists.txt`

- [x] add a `decomp_gather` interface (`_2d_r`, `_3d_r`) = `decomp_redistribute` transposed: gather
      `rbuf(k) = buf_writer(recv_target(k))`; one `MPI_Alltoallv` with **send/recv swapped** (`recvcounts/rdispls`
      as send, `sendcounts/sdispls` as recv); unpack `field_owned(i) = sbuf(send_pos(i)+1)`
- [x] reuse the **same** `t_io_decomp` plan — no new plan, no extra communication; real WP only (no int variant)
- [x] unit test: deterministic `field_owned` (`real(myList(i),WP)*1.5-0.25`, +per-level offset for 3D; NO RNG, per
      reproducibility) → `decomp_redistribute` → `decomp_gather` reproduces it **identically** (np=1 self-copy + np=2);
      last-chunk padding ignored
- [x] register ctest `test_io_decomp_gather` (np1, np2)
- [x] **GATE:** ctest GREEN — `gather ∘ redistribute == identity` on owned entities — PASS np=1 AND np=2 (Intel dp,
      worktree build); both report `max|Δ|=0` over 2D+3D across 6 (N,C,n_writers) combos incl. partial last chunks,
      writer subsets, and real cross-rank exchange at np=2

### Stage 2 — Shared coord/attr helper (DRY)
*Factor the ushow-viewability embedding so restart stores are byte-identical in shape to M9 output stores.*

#### Task 2.1 ✅: Extract the coord + `_ARRAY_DIMENSIONS` + UGRID-attr embedding from `mod_io_means`

**Files:**
- Create: `src/io/mod_io_coords.F90`
- Modify: `src/io/mod_io_means.F90`

- [x] extract the per-store `lon`/`lat` + `_ARRAY_DIMENSIONS` + UGRID-attr embedding (incl. the elem-centroid
      coordinate path) into a shared routine `(store, entity, mesh, partit, decomp)` → writes coords + attrs —
      4 VERBATIM routines in `mod_io_coords`: `io_coords_compute(entity,mesh,nO,lon,lat,rlon,rlat,nlev)` (node coords
      AND elem-centroid r2g path), `io_coords_init_lonlat(a_lon,a_lat,N,C)`, `io_coords_define_lonlat(store,a_lon,a_lat,hdim)`
      (the `_ARRAY_DIMENSIONS`+CF attrs), `io_coords_put(D,store,arr,owned)` (redistribute+write the coord chunks).
      `nO` is the caller's partit-derived local count; `partit` info enters via `nO`/`decomp`. (auto-built by the
      `src/io/*.F90` GLOB; Fortran scanner orders `mod_io_coords` before `mod_io_means`.)
- [x] `mod_io_means` calls the shared helper (pure refactor, no behavior change) — `means_init` calls `io_coords_compute`
      ×2 (node+elem), `def_field_store` calls `io_coords_init_lonlat`+`io_coords_define_lonlat`, `open_field_store`
      calls `io_coords_put` ×2; the old `put_static` + inline coord/attr blocks removed; unused `r2g` import dropped
- [x] **GATE:** `run_output_gate.sh` + `run_meshdiag_gate.sh` + `run_zarrsmoke.sh` stay byte-identical
      (`zarr_diff.py` `max|Δ|=0`); `ctest` unchanged GREEN — ALL GREEN (Intel dp, worktree build): zarrsmoke
      max|Δ|=0; meshdiag np=1 lon/lat max|Δ|=0 (18 vars, 0 failures); output np 1/2/8 + Task-2.6 knob sweep all
      max|Δ|=0 + partition-independent. **Plus a direct byte-diff: pre-refactor `main` build vs post-refactor
      worktree build of `fesom_outputsmoke` → BYTE-IDENTICAL store trees at np=1 AND np=2** (raw chunk bytes +
      `.zarray`/`.zattrs` + lon/lat/time coords). ctest 20/20 PASS (incl. `test_io_means_np1`)

### Stage 3 — `mod_io_restart` WRITE path
*Per-field M9-shape snapshot stores inside an immutable, atomically-finalized checkpoint folder.*

#### Task 3.1: C filesystem shims — `rename` / `unlink` / `fsync` (`bind(C)`)

**Files:**
- Modify: `src/io/mod_io_zarr.F90` (or a small new `src/io/mod_io_posix.F90`)

- [ ] add `bind(C)` shims for `rename(2)`, `unlink(2)`/`rmdir(2)`, `fsync(2)` mirroring the existing `c_mkdir`
      (`mod_io_zarr.F90:38`) — do **NOT** use `execute_command_line('mv'/'rm')`, which forks and **segfaults after
      MPI_Init** (the M9 lesson that forced `c_mkdir`)
- [ ] a recursive-delete helper for a multi-file Zarr tree (used by tmp-cleanup + keep-N prune) on top of the shims
- [ ] **GATE:** a tiny driver test creates a dir tree, `rename`s it, recursively `unlink`s it — no segfault, exit 0

#### Task 3.2: Writer mechanism — one field end-to-end + folder + `checkpoint.json`

**Files:**
- Create: `src/io/mod_io_restart.F90`

- [ ] a field descriptor (name, units, entity NODE/ELEM, **level kind `nl` vs `nl-1`**, precision, array pointer) + registry
- [ ] `restart_write_field`: `decomp_redistribute` (compute→canonical) → `zarr_write_chunk` per owned chunk into a
      single-variable single-entity **snapshot** store (`<f8`, codec from namelist), coords via the Stage 2 helper
- [ ] checkpoint folder `fesom.<YYYY>.<DDD>.<SSSSS>/` + `checkpoint.json`
      (`format_version, year, day, time_sec, globalstep, fesom_git, npes_wrote`, rank-0)
- [ ] **GATE (smoke):** write one field; `zarr_diff.py`/xarray opens the store (dims, finite `lon`/`lat`,
      `_ARRAY_DIMENSIONS`); `ushow <store> -m fesom.mesh.diag.zarr` noted (manual)

#### Task 3.3: Atomicity + `restart.latest` + keep-N prune (crash-injection gate)

**Files:**
- Modify: `src/io/mod_io_restart.F90`

- [ ] sequence: all writers → `fesom.<tag>.tmp/` → `MPI_Barrier` → rank-0 `fsync`+`rename` to final → rank-0
      atomic-update `restart.latest` (write `.tmp`+`rename`) → rank-0 keep-N prune (`restart_keep`; warn, never abort)
- [ ] **GATE (crash-safety):** inject a stray `*.tmp/` and a finalized-but-not-pointed checkpoint → the reader still
      follows the previous valid `restart.latest`; assert the pointer flips atomically (never a partial `restart.latest`)

#### Task 3.4: Register the full field set (oce + ice incl. sigma)

**Files:**
- Modify: `src/io/mod_io_restart.F90`

- [ ] oce node: `eta_n`, `hbar`, `ssh_rhs_old`, `hnode` (`nl-1`); `w`/`w_e`/`w_i` (**full levels `nl`** — the M9
      `on_full_levels` distinction; wrong level count breaks the shape/gate)
- [ ] tracers `temp`/`salt`(+passive): `values`, `valuesAB`, **`valuesold` M1 mandatory** (`valuesold(1,:,:)`, AB2
      history — oracle serializes it unconditionally `io_restart.F90:224`); **M2 only `if AB_order==3`** (`:225-226`;
      this driver hardcodes `AB_order=2`, so M2 is dead but the conditional must be there)
- [ ] oce element: `uv`→`u`/`v`, `uv_rhsAB`→`urhs_AB`/`vrhs_AB` (+`urhs_AB3`/`vrhs_AB3` if `AB_order==3`)
- [ ] optional `tke` (`dyn%work%tke`, when `mix_scheme==5`) — guarded
- [ ] ice node: `area`/`hice`/`hsnow` (`ice%data(1:3)%values`), `uice`, `vice`; ice element: `sigma11`/`sigma12`/`sigma22`
      (`ice%work%`) — **F-C**
- [ ] **GATE (smoke):** a checkpoint contains all expected stores; xarray shapes match `nod2D`/`elem2D` × correct
      level count; element stores ≈ 2× node count (as in M9)

### Stage 4 — `mod_io_restart` READ path
*The transpose of the write path; restores owned values then halo-exchanges to bit-exact full arrays.*

#### Task 4.1: Per-field reader (`zarr_read_chunk` + `decomp_gather` + halo exchange) + safety check

**Files:**
- Modify: `src/io/mod_io_restart.F90`

- [ ] `restart_read`: resolve the folder via `restart.latest` under `RestartInPath`; read `checkpoint.json`;
      **abort** only on missing/corrupt `restart.latest`/store; on `checkpoint.json` time ≠ clock **WARN and continue**
      (a legitimate `dt`-change restart trips it — oracle `io_restart.F90:904-914` warns, does **not** abort)
- [ ] `restart_read_field`: each writer `zarr_read_chunk` its owned chunks → `buf_writer` → `decomp_gather` →
      `field_owned` → write back **into the live array** → halo exchange using the **same variant the field's in-step
      consumer uses**: `exchange_nod` (nodes); `exchange_elem` for `uv` (`com_elem2D`, `mod_halo.F90:110`),
      `exchange_elem_full` for element fields needing `eDim+eXDim` — a wrong variant leaves stale halo cells
- [ ] node + element + 3D coverage; ice incl. sigma
- [ ] **GATE:** in-process write→read round-trip reproduces every registered array `max|Δ|=0` over the **full local
      extent (incl. eXDim)** — so a wrong exchange variant is caught — including halos after exchange

### Stage 5 — Lifecycle wiring + cadence
*Three hooks PLUS three restart-conditioning fixes the cold-start path currently hardwires. The conditioning is as
load-bearing as the I/O — each of the three, if wrong, makes the Stage 6 gate diverge or the read hook dead (all three
found by plan-review, verified against source).*

#### Task 5.1: Restart-mode detection + skip cold-start `.clock` overwrite + READ hook

**Files:**
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90`
- Modify: `src/io/mod_io_restart.F90`

- [ ] detect restart mode **early** — before the cold-start `.clock` writer (`:574-580`): restart iff `FESOM3_RESTART`
      set AND `restart.latest` exists under `RestartInPath`
- [ ] **skip the unconditional `.clock` overwrite (`:574-580`) when in restart mode** so `clock_init` (`:582`) reads
      the previous segment's `clock_finish`-written `.clock` and sets `r_restart=.true.` — else the overwrite writes
      two equal lines and `r_restart` is **always false** (verified `:574-582` + `mod_clock.F90:126-131`), making the
      read hook dead
- [ ] **READ hook** (after `clock_init`/`clock_nsteps` `:582-583`, after cold-state init, before the loop `:768`):
      `if (r_restart) call restart_read(...)` — resolve folder via `restart.latest`, `decomp_gather`+halo-exchange
      into the live arrays
- [ ] **GATE:** with a checkpoint present, `r_restart` becomes true and `restart_read` runs; with none, cold start unchanged

#### Task 5.2: First-resumed-step AB guard (`lfirst .and. .not. r_restart`)

**Files:**
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90`

- [ ] change the step call (`:793`) `step_oce(n, dt, (n == 1), ...)` → `step_oce(n, dt, (n==1) .and. .not. r_restart, ...)`
      so the first resumed step blends with `ff=ab2` (AB2) over the **restored** `uv_rhsAB`, not forward-Euler
      (`ff=1.0`, `oce_dyn_velrhs.F90:136`) which discards it — the routine's header documents exactly this guard
      (`oce_dyn_velrhs.F90:62-64`)
- [ ] confirm no tracer analog needed (FESOM3 tracer solve has no `lfirst`; FESOM2's only `.not.r_restart` tracer
      guard is the cold-start `valuesold=values` init, which the post-cold-init READ-hook ordering overwrites)
- [ ] **GATE:** folded into 6.1 — the velocity dump on the **first resumed step** is where a wrong guard surfaces

#### Task 5.3: WRITE hooks + cadence + `clock_finish` per write

**Files:**
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90`
- Modify: `src/io/mod_io_restart.F90`

- [ ] **REGISTER** (~`:725`): env-gate `FESOM3_RESTART=<dir>` → `RestartOutPath`; parse `&nml_restart`
      (`restart_length`, `restart_length_unit`, `restart_keep`, `compressor`, `n_writers`) + `FESOM3_*` overrides;
      register the field set
- [ ] **WRITE** in-loop (after `step_oce`, ~`:802`): `if (restart_due(n)) call restart_write(...)`; end-of-run
      (after loop ~`:848`): final `restart_write` (skip if the last step already wrote)
- [ ] **`clock_finish` + `restart.latest` update after EVERY `restart_write`** (periodic AND final), mirroring
      `io_restart.F90:573-578` — NOT only at end-of-run, else a periodic mid-run checkpoint leaves `.clock` at the
      cold-start time → `r_restart` stays false next launch → the periodic checkpoint is silently ignored
- [ ] `restart_due(n)` = **global** `is_due(restart_length_unit, restart_length, n)` — port FESOM2 `is_due` +
      `annual_event`/`step_event` **fresh** (`event_due` is private + coupled to `t_means_clock`, not reusable)
- [ ] **GATE:** a 2-segment run (seg 1: K steps → checkpoint + `.clock`; seg 2: resumes, reads, runs N−K) completes;
      `.clock` chaining + `r_restart` fire; a PERIODIC mid-run checkpoint also yields a resumable `.clock`

### Stage 6 — Reproducibility gates

#### Task 6.1: PRIMARY self-consistency gate (straight-through vs split, np1/np2)

**Files:**
- Create: `tools/run_restart_gate_core2.sh`

- [ ] straight-through: run N steps, `FESOM_DUMP_ALL` final state (node + element incl. ice + sigma)
- [ ] split: run K → `restart_write` → **fresh process** → `restart_read` (`.clock`-chained) → run N−K → `FESOM_DUMP_ALL`
- [ ] `dump_diff.py --glob` → `max|Δ|=0` (whole state incl. ice); at np=1 and np=2
- [ ] **boundary variant:** a second split whose checkpoint lands on a **year/month boundary** (or runs long enough
      that the periodic cadence fires), so forcing-from-clock resume + the `roll_monthly_clim` read-ahead
      (`:1016-1020`, fires at `timenew==86400` / `n==1`) is actually exercised — a mid-run integer K never hits it
- [ ] **GATE:** `run_restart_gate_core2.sh` GREEN (`max|Δ|=0`) on Intel **and** GNU, **both** the mid-run and boundary splits

#### Task 6.2: SECONDARY cross-partition gate — RESTORE round-trip (np=2 write → np=8 restore)

**Files:**
- Create: `tools/run_restart_gate_multirank.sh`

- [ ] write checkpoint `C2` at np=2; a **fresh np=8 process** reads `C2` and **immediately re-writes** `C8`
      (**zero steps**); assert `C2 ≡ C8` canonically (`zarr_diff.py`, `max|Δ|=0`)
- [ ] this proves F-B = partition-independent **STORAGE + RESTORE**, NOT cross-np **evolution** byte-identity — which
      is physically impossible in FESOM (SSH-CG `allreduce_sum` reduces in comm-size-dependent order, `oce_ssh_solve.F90`
      / `mod_halo.F90:25-30`; FESOM2 has the same property; every project byte-gate is same-np). The bit-identical
      **resume** claim is **same-np** (Task 6.1); cross-np gives a physically-valid, not bit-identical, continuation
- [ ] `env.sh` KNEM `single_copy_mechanism=none`
- [ ] **GATE:** `run_restart_gate_multirank.sh` GREEN (`C2 ≡ C8`, `max|Δ|=0`)

#### Task 6.3: No-regression sweep + ctest

- [ ] `ctest` GREEN both compilers (new `test_io_decomp_gather` included)
- [ ] `run_output_gate.sh` + `run_meshdiag_gate.sh` + `run_zarrsmoke.sh` GREEN (shared coord helper didn't regress M9)
- [ ] production MR lifecycle byte-gate (`run_lifecycle_fullynative_gate_multirank.sh`, forced np=2) stays `max|Δ|=0`
- [ ] **GATE:** full sweep GREEN

### Final

#### Task F1: Verify acceptance criteria
- [ ] every Overview requirement met (write periodic+end, read on `r_restart`, byte-exact resume incl. ice,
      partition-independent, ushow-viewable)
- [ ] both gates `max|Δ|=0`; smoke stores open in xarray/ushow

#### Task F2: Docs + memory + tag
- [ ] update `docs/HANDOFF.md` ("Where we are" + "Next task" → restart DONE / next milestone)
- [ ] update memory (`project-fesom3-implementation-state`); tag the milestone
- [ ] move this plan to `docs/plans/completed/`

## Technical Details

### On-disk layout
```
<RestartPath>/
  fesom.<YYYY>.<DDD>.<SSSSS>/      # immutable checkpoint; name = year.day-of-year.sec-of-day (lexically=chronological)
      eta_n.zarr/  hbar.zarr/  ssh_rhs_old.zarr/  hnode.zarr/        # node
      temp.zarr/ temp_AB.zarr/ salt.zarr/ ...                        # node tracers (+ valuesold M1/M2 optional)
      u.zarr/ v.zarr/ urhs_AB.zarr/ vrhs_AB.zarr/ (urhs_AB3 ...)     # element
      w.zarr/ w_e.zarr/ w_i.zarr/                                    # node
      area.zarr/ hice.zarr/ hsnow.zarr/ uice.zarr/ vice.zarr/        # ice node
      sigma11.zarr/ sigma12.zarr/ sigma22.zarr/                      # ice element (F-C)
      checkpoint.json   # {format_version, year, day, time_sec, globalstep, fesom_git, npes_wrote}
  fesom.clock                      # existing FESOM clock (cold-vs-restart detector) — role unchanged
  restart.latest                   # one line: name of the newest finalized checkpoint folder
```
Each `*.zarr` is a single-variable, single-entity, **snapshot** store (no time dimension) — an M9 field store in
every respect, so `ushow`/xarray render it verbatim. History = separate per-checkpoint folders (not records in a
file): immutable, atomic, trivially prunable, branchable. `restart.latest` is the resume pointer (the clock provides
the time + the `checkpoint.json` safety check); naming is cosmetic. Folder names use `int(timenew)` for the
sec-of-day so they are exact integers.

**Time convention (one, end-to-end — per plan-review).** `checkpoint.json.time_sec`, the folder tag, and the
read-time safety compare ALL use the **new**-clock time (`timenew + (daynew−1)*86400`, the value `clock_init` reads
back from `.clock` line 2). Do **not** mix in the oracle's old-time `ctime` (`timeold + (dayold−1)*86400`,
`io_restart.F90:816`) — pick new and use it everywhere so the compare is self-consistent.

**v1 scope (per plan-review).** `globalstep` is stored in `checkpoint.json` for provenance only — nothing in FESOM3
keys off an absolute cumulative step (forcing keys off the clock, `nsteps` is recomputed each segment by
`clock_nsteps`), so it is not read back. **Concurrent field-output + mid-year restart is OUT OF SCOPE for v1:** the
M9 output record counter is in-memory (`mod_io_means.F90:125`, reset to 0 in a fresh process and used to size the
year store), so a resumed process would mis-number/overwrite the current year's output store. Restart and output are
validated independently; document the limitation.

### `decomp_gather` (the inverse of `decomp_redistribute`)
`decomp_redistribute` (write): `sbuf(send_pos(i)+1)=field(i)` → `Alltoallv(send→recv)` →
`buf_writer(recv_target(k))=rbuf(k)`. `decomp_gather` (read) runs it backwards on the **same** `t_io_decomp`:
`rbuf(k)=buf_writer(recv_target(k))` → `Alltoallv(recv→send)` (counts/displs swapped) → `field(i)=sbuf(send_pos(i)+1)`.
3-D loops levels (reusing the 2-D plan), as in the existing code. Store holds only **canonical owned** values; halos
are filled post-read by `exchange_nod`/`exchange_elem` as exact copies → byte-exact incl. halos. Partition-independence
falls out: a different np builds a different plan over the **same** canonical chunks.

### Field set → FESOM3 source mapping
| store | FESOM3 array | entity | dims |
|---|---|---|---|
| `eta_n` / `hbar` / `ssh_rhs_old` | `dyn%eta_n` / `mesh%hbar` / `dyn%ssh_rhs_old` | node | 2-D |
| `hnode` | `mesh%hnode` | node | 3-D |
| `u` / `v` | `dyn%uv(1,:,:)` / `dyn%uv(2,:,:)` | element | 3-D |
| `urhs_AB` / `vrhs_AB` (+`*_AB3`) | `dyn%uv_rhsAB(1,1,:,:)` / `(1,2,:,:)` (`(2,·,·)` if `AB_order==3`) | element | 3-D |
| `temp` / `salt` (+passive) | `tracers%data(j)%values` | node | 3-D |
| `temp_AB` / … | `tracers%data(j)%valuesAB` | node | 3-D |
| `temp_M1` / `_M2` (optional) | `tracers%data(j)%valuesold(1\|2,:,:)` | node | 3-D |
| `w` / `w_expl` / `w_impl` | `dyn%w` / `dyn%w_e` / `dyn%w_i` | node | 3-D |
| `tke` (optional, `mix_scheme==5`) | `dyn%work%tke` | node | 3-D |
| `area` / `hice` / `hsnow` | `ice%data(1:3)%values` | node | 2-D |
| `uice` / `vice` | `ice%uice` / `ice%vice` | node | 2-D |
| `sigma11` / `sigma12` / `sigma22` | `ice%work%sigma11/12/22` | element | 2-D |

### `namelist.io` `&nml_restart` schema
`restart_length` (int), `restart_length_unit` (`y\|m\|d\|h\|s\|off`, default `y`/`1`), `restart_keep` (int; `0`=keep
all), `compressor` (`none\|lz4`), `n_writers` (int). Staged into the rundir like the M9 `namelist.io`; `FESOM3_*` env
overrides.

### ushow contract (the consumer)
`.zgroup`; 1-D `lon`/`lat` of size `n_points` with `units`; data vars carry `_ARRAY_DIMENSIONS`; codec **null/lz4**
(NOT zlib); dtype `<f8`/`<f4`. Topology via `-m fesom.mesh.diag.zarr` (reuse M9's; connectivity is optional so stores
point-render standalone; emit `fesom.mesh.diag.zarr` once if output is off). **No CLI var-selector** and one `-m`
coordinate set → each store must be single-entity (why per-field stores, not grouped `oce.zarr`/`ice.zarr`).
Glob-animate: `ushow "fesom.*/temp.zarr" -m fesom.mesh.diag.zarr` (folders sort chronologically; snapshots = frames).

### Gate commands
- Build: `./configure.sh --compiler intel --precision dp --clean --build` (+ `--compiler gnu`).
- ctest: `cd build_intel_dp && ctest --output-on-failure`.
- Restart gates: `tools/run_restart_gate_core2.sh` (split np1/np2), `tools/run_restart_gate_multirank.sh` (np2→np8).
- Diff: `/work/ab0995/a270088/mambaforge/bin/python3 tools/dump_diff.py --glob …` (`max|Δ|=0`); `tools/zarr_diff.py`.
- Multi-rank: `env.sh` (KNEM `single_copy_mechanism=none`).

### Error handling
- Missing/corrupt checkpoint or `restart.latest` → abort with a clear message (as `clock_init` does for a missing `.clock`).
- `checkpoint.json` time ≠ `.clock` → **warn and continue** (oracle `io_restart.F90:904-914` warns — a legitimate
  `dt`-change restart trips it; never abort on this).
- A stray `.tmp/` from a crashed write is ignored (the reader only ever follows `restart.latest`).
- `restart_keep` prune failure warns, never aborts.
- **np mismatch is not an error** — it is the design goal (F-B).

## Post-Completion
*Informational — external/manual, no checkboxes.*

**Manual verification**
- `ushow` display smoke on a checkpoint store (`ushow <checkpoint>/temp.zarr -m fesom.mesh.diag.zarr`) + the
  glob-animate-across-checkpoints check (the xarray/`zarr_diff.py` round-trip is the automated proxy).

**Performance**
- Measure checkpoint write time + filesystem load at `dist_512`/`dist_864`; tune `n_writers`/chunk shape; compare
  `none` vs `lz4`. Watch inode pressure from per-field stores × kept checkpoints (same profile M9 accepted).
