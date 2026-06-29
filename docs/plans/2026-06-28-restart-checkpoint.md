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
>
> **Claim scope refinement (Task 6.1 gate result, 2026-06-29).** Same-np resume is **exactly `max|Δ|=0`** at **np=1
> for the entire state**, and at **np>1 for the whole state in the quiet (mid-run) case**. At **np>1 with active
> forcing** there is a residual **≤1 ULP** (temp `2.22e-16`, w `4e-22`; propagated from element velocity into temp/w):
> FESOM2 assigns an element to a rank if **any** of its nodes is owned (`gen_comm.F90:265`), so on CORE2/dist_2 **562
> boundary elements are redundantly owned by both ranks** (Σ myDim_elem2D 245221 > nElem2D 244659). `exchange_elem(UV)`
> uses `com_elem2D` (eDim halo) only — FESOM2's UV has **no eXDim slot** (`oce_setup_step.F90:655`) — so these
> node-only-adjacent shared elements are **never synced** and diverge ~1 ULP from order-dependent RHS sums, in **both
> FESOM2 and FESOM3**. The partition-independent checkpoint must dedup them to one value, so resume perturbs the other
> rank's copy by ~1 ULP. **This is FESOM2-faithful** — FESOM2's own netCDF restart `gather_elem3D` deduplicates
> identically. Forcing a live-model sync (eXDim halo + `exchange_elem_full(UV)`) would make FESOM3 diverge from FESOM2
> and **break the M0–M9 byte-identity gates**, so it is deliberately NOT done. Bit-exact same-np resume at np>1 would
> require FESOM2-`raw`-style per-rank element storage (rejected — costs cross-np element restore). The Task 6.1 gate
> therefore asserts **strict `max|Δ|=0` at np=1** and admits a **per-store relative ε floor (1e-12, ≫ the ~2e-16
> residual, ≪ any real regression at ~1e-4)** at np>1. See memory `restart-np-element-1ulp-inherent`.

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

#### Task 3.1 ✅: C filesystem shims — `rename` / `unlink` / `fsync` (`bind(C)`)

**Files:**
- Create: `src/io/mod_io_posix.F90` (new dedicated POSIX-fs module; auto-built by the `src/io/*.F90` GLOB)
- Modify: `src/io/mod_io_zarr.F90` (export `zarr_mkdir` — the shared mkdir -p the gate/`mod_io_restart` reuse)
- Create: `test/test_io_posix.F90`
- Modify: `test/CMakeLists.txt`

- [x] add `bind(C)` shims for `rename(2)`, `unlink(2)`/`rmdir(2)`, `fsync(2)` mirroring the existing `c_mkdir`
      (`mod_io_zarr.F90:38`) — do **NOT** use `execute_command_line('mv'/'rm')`, which forks and **segfaults after
      MPI_Init** (the M9 lesson that forced `c_mkdir`). New `mod_io_posix.F90` binds libc `rename`/`unlink`/`rmdir`/
      `open`/`fsync`/`close` **directly via `iso_c_binding`** (SAME mechanism as `c_mkdir` — the repo compiles **zero**
      `.c` files); public thin wrappers `posix_rename`/`posix_unlink`/`posix_rmdir`/`posix_fsync_dir` each take a
      `character(len=*)`, append `c_null_char`, return the C int status. `posix_fsync_dir` = `open(O_RDONLY)`+`fsync`+
      `close` (directory durability before rename)
- [x] a recursive-delete helper for a multi-file Zarr tree (used by tmp-cleanup + keep-N prune) on top of the shims —
      `posix_rmtree(path)` binds libc `nftw(3)` with `FTW_DEPTH|FTW_PHYS` (post-order) and a `bind(C)` Fortran callback
      passed via `c_funloc` that `remove()`s each entry bottom-up (no `.c` file, no `execute_command_line`); returns 0
      only if `nftw` AND every `remove()` succeeded (module-`save` failure counter, serial rank-0 use)
- [x] **GATE:** a tiny driver test creates a dir tree, `rename`s it, recursively `unlink`s it — no segfault, exit 0 —
      **PASS** (Intel dp, worktree `build_intel_dp`): `test/test_io_posix.F90` builds `scratch/top/a/b/c.bin`+sibling,
      `posix_fsync_dir`→`posix_rename`(top→renamed)→`posix_rmtree`(renamed)+unlink/rmdir round-trip, all shims return 0,
      old path gone after rename, tree gone after rmtree (re-rmtree→ENOENT confirms), no scratch dir left behind.
      `ctest -R io_posix` GREEN (`test_io_posix_np1 ... Passed`, ctest exit 0); direct driver prints
      `ALL PASS (rename/unlink/rmdir/fsync/rmtree, no fork, exit 0)`, **exit code 0** (no segfault)

#### Task 3.2 ✅: Writer mechanism — one field end-to-end + folder + `checkpoint.json`

**Files:**
- Create: `src/io/mod_io_restart.F90`
- Create: `src/drivers/fesom_restartsmoke.F90` (➕ the Stage-3 GATE driver, mirrors `fesom_outputsmoke`)
- Create: `tools/run_restartsmoke.sh` (➕ the gate runner; `F3`/`BUILD`/`RUN` overridable for worktree runs)
- Modify: `tools/zarr_diff.py` (➕ `--restart` mode: snapshot-store + `checkpoint.json` verify)

- [x] a field descriptor (name, units, entity NODE/ELEM, **level kind `nl` vs `nl-1`**, precision, array pointer) + registry
      — `t_restart_field` (`p2d`/`p3d` live POINTERS, `entity`, `ndim`, `on_full_levels`, `dtype` default `<f8`,
      `nlev`/`hdim`/`vdim`) + `t_restart` registry (`f(RESTART_MAXF=64)`, node `Dn`+elem `De` decomps, cached
      node/elem coords + `depth_nz`/`nz1`); `restart_init` (builds decomps + coords once, like `means_init`) +
      `restart_register_field(R, name, units, entity, p2d=/p3d=, on_full_levels, precision)` (pointer dummies →
      strided model sections associate without a copy, stay live for Task 4.1 read-back)
- [x] `restart_write_field`: `decomp_redistribute` (compute→canonical) → `zarr_write_chunk` per owned chunk into a
      single-variable single-entity **snapshot** store (`<f8`, codec from namelist), coords via the Stage 2 helper
      — general node+elem / 2-D `(entity)` + 3-D `(nlev, entity)` (NO time dim, NO mean accumulation); rank-0
      define → barrier → every writer writes its chunks (M9 store-create ordering); 3-D embeds the `nz`/`nz1`
      vertical coord; `io_coords_*` for the lon/lat + `_ARRAY_DIMENSIONS` + UGRID embed
- [x] checkpoint folder `fesom.<YYYY>.<DDD>.<SSSSS>/` + `checkpoint.json`
      (`format_version, year, day, time_sec, globalstep, fesom_git, npes_wrote`, rank-0) — `restart_write` zero-pads
      the folder (`int(time_sec)`=sec-of-day for `SSSSS`), writes each registered field's store, then rank-0 writes
      `checkpoint.json` via plain Fortran string writes (`fesom_git` from `$FESOM3_GIT` else `unknown`). Atomic
      tmp/rename + `restart.latest` deferred to Task 3.3 (plain `zarr_mkdir` + direct write here)
- [x] **GATE (smoke):** write one field; `zarr_diff.py`/xarray opens the store (dims, finite `lon`/`lat`,
      `_ARRAY_DIMENSIONS`); `ushow <store> -m fesom.mesh.diag.zarr` noted (manual) — **PASS np=1 AND np=2** (Intel dp,
      worktree `build_intel_dp`): `fesom_restartsmoke` writes `fesom.2000.001.03600/eta_n.zarr` (`value(g)==g`) +
      `checkpoint.json`; xarray opens `eta_n` dims `(nod2,)`=3140 (NO time dim), embedded lon/lat shape `(3140,)`
      finite, `_ARRAY_DIMENSIONS=['nod2']` (read raw off disk), `value(g)==g` `max|Δ|=0`; `checkpoint.json` parses
      with all 7 keys + correct clock (`npes_wrote=2` at np=2). **np1 vs np2 store tree BYTE-IDENTICAL** (data+coord
      chunks + `.zarray`/`.zattrs`) — canonical write is partition-independent. ushow command printed by the runner
      (manual display; the xarray round-trip is the automated proxy)

#### Task 3.3 ✅: Atomicity + `restart.latest` + keep-N prune (crash-injection gate)

**Files:**
- Modify: `src/io/mod_io_restart.F90`
- Modify: `src/io/mod_io_posix.F90` (➕ `posix_listdir` — opendir/readdir/closedir `bind(C)`, the no-fork
      directory enumeration the keep-N prune needs to find `fesom.*` checkpoint folders; same module as Task 3.1)
- Create: `src/drivers/fesom_restartcrash.F90` (➕ the Task-3.3 GATE driver; mirrors `fesom_restartsmoke`)
- Create: `tools/run_restartcrash.sh` (➕ the gate runner; `F3`/`BUILD`/`RUN` overridable for worktree runs)

- [x] sequence: all writers → `fesom.<tag>.tmp/` → `MPI_Barrier` → rank-0 `fsync`+`rename` to final → rank-0
      atomic-update `restart.latest` (write `.tmp`+`rename`) → rank-0 keep-N prune (`restart_keep`; warn, never abort)
      — `restart_write` now: rank-0 `posix_rmtree` any stale same-tag `.tmp` + `zarr_mkdir(fesom.<tag>.tmp/)` →
      barrier → every writer stages its stores into the tmp + rank-0 stages `checkpoint.json` → barrier → rank-0
      `posix_fsync_dir(tmp)` + `posix_rename(tmp→fesom.<tag>/)` (ATOMIC publish) + `posix_fsync_dir(parent)` →
      `update_restart_latest` (write `restart.latest.tmp` with the bare folder NAME + `posix_rename` onto
      `restart.latest`) → `restart_prune`. Knob `restart_keep` (arg + `FESOM3_RESTART_KEEP`; `0`=keep all) added to
      `restart_init`. Prune: `posix_listdir` → strict `is_checkpoint_name` filter (`fesom.`+4+`.`+3+`.`+5 digits, 20
      chars — excludes `fesom.clock`/`*.zarr`/`restart.latest`/`*.tmp`) → lexical(==chrono) sort → `posix_rmtree` the
      oldest beyond `restart_keep`, but the **current pointer target is always protected** (so a stray later-named
      folder can never cause the just-committed checkpoint to be pruned); every prune failure WARNs, never aborts.
      Resolver `restart_resolve_latest(restart_dir, folder_out, ok)` (public, for Task 4.1) reads `restart.latest`,
      returns `<dir>/<name>`, `ok=.true.` iff the target carries a `checkpoint.json`; it **only follows the pointer**,
      never scans, so stray `*.tmp/` + unpointed `fesom.*/` are ignored by construction.
- [x] **GATE (crash-safety):** inject a stray `*.tmp/` and a finalized-but-not-pointed checkpoint → the reader still
      follows the previous valid `restart.latest`; assert the pointer flips atomically (never a partial `restart.latest`)
      — **PASS np=1 AND np=2** (Intel dp, worktree `build_intel_dp`): `fesom_restartcrash` writes C1
      (`fesom.2000.001.03600`) atomically → `restart.latest`=="fesom.2000.001.03600" + resolve==C1; injects a stray
      `fesom.2000.001.05000.tmp/` + an unpointed LATER `fesom.2000.001.09000/` → resolve STILL returns C1 (follows the
      pointer, ignores both); writes C2 (`fesom.2000.001.07200`, `restart_keep=1`) → `restart.latest`=="…07200" +
      resolve==C2, keep-N=1 pruned C1 (its `checkpoint.json` gone), C2 remains (protect-target despite 09000>07200),
      the later folder survives, the stray `.tmp` is left untouched. All 9 assertions PASS at np 1 AND 2, exit 0; dir
      listing after = `{05000.tmp, 07200, 09000, restart.latest}` (C1 gone). Runner: `tools/run_restartcrash.sh`.
      No-regression: `run_restartsmoke.sh` still GREEN (same final folder, max|Δ|=0); `ctest -R io_` 7/7 PASS.

#### Task 3.4 ✅: Register the full field set (oce + ice incl. sigma)

**Files:**
- Modify: `src/io/mod_io_restart.F90`
- Create: `src/drivers/fesom_restartstate.F90` (➕ the Task-3.4 GATE driver; allocates correctly-shaped REAL
      `t_dyn`/`t_tracer`/`t_ice` + `mesh%hbar`/`hnode`, fills owned slots with a partition-independent formula,
      registers via `restart_register_state`, writes one full checkpoint)
- Modify: `tools/zarr_diff.py` (➕ `--restart-state` mode: verify the full oce+ice store set — shapes/levels/coords/
      values + elem≈2× node)
- Create: `tools/run_restartstate.sh` (➕ the gate runner; AB2+AB3 × np1/np2 + partition-independence compare;
      `F3`/`BUILD`/`RUN` overridable for worktree runs)

- [x] oce node: `eta_n`, `hbar`, `ssh_rhs_old`, `hnode` (`nl-1`); `w`/`w_e`/`w_i` (**full levels `nl`** — the M9
      `on_full_levels` distinction; wrong level count breaks the shape/gate). `hbar`/`hnode` are `real(MP)` (see
      below); registered via `restart_register_field_mp` (lossless MP→WP staging). Store names per the plan
      field-set table (`eta_n`, not the oracle's `ssh`; `w`/`w_expl`/`w_impl`).
- [x] tracers `temp`/`salt`(+passive): `values`, `valuesAB`, **`valuesold` M1 mandatory** (`valuesold(1,:,:)`, AB2
      history — oracle serializes it unconditionally `io_restart.F90:224`); **M2 only `if tracers%data(j)%AB_order==3`**
      (`:225-226`). Names from a `tracer_name(id,j)` helper == the oracle `ini_ocean_io` CASE (1=temp, 2=salt, passive
      by ID, default `tra<j>`). M2 gated by the AB-order driver run (`FESOM3_AB_ORDER=3` ⇒ `temp_M2`/`salt_M2` appear).
- [x] oce element: `uv(1,:,:)`→`u`, `uv(2,:,:)`→`v`, `uv_rhsAB(1,1|1,2,:,:)`→`urhs_AB`/`vrhs_AB` (+`urhs_AB3`/`vrhs_AB3`
      from `uv_rhsAB(2,·,·)` iff `dyn%AB_order==3`) — strided pointer sections associate without a copy
- [x] optional `tke` (`dyn%work%tke`, FULL levels) — guarded `if (mix_scheme==5 .and. allocated(dyn%work%tke))`;
      `mix_scheme` is an optional arg to `restart_register_state` (kept decoupled from `mod_param_phys`; Stage 5 passes
      `mix_scheme_nmb`). Gate exercises it ON (`mix_scheme=5`).
- [x] ice node: `area`/`hice`/`hsnow` (`ice%data(1:3)%values`), `uice`, `vice`; ice element: `sigma11`/`sigma12`/`sigma22`
      (`ice%work%sigma11/12/22` — verified `mod_ice.F90:37`) — **F-C**
- [x] **GATE (smoke):** a checkpoint contains all expected stores; xarray shapes match `nod2D`/`elem2D` × correct
      level count; element stores ≈ 2× node count — **PASS np=1 AND np=2** (Intel dp, worktree `build_intel_dp`):
      `fesom_restartstate` writes `fesom.2000.001.03600/` with all **26** stores (AB2+tke); `zarr_diff.py
      --restart-state` confirms every store's entity×level-kind shape + `_ARRAY_DIMENSIONS` `(vdim,hdim)`, embedded
      finite lon/lat sized to the entity, monotonic positive-down `nz`/`nz1`, finite data, and value==formula
      (`g` 2-D / `g+0.5L` 3-D) **max|Δ|=0 on ALL 26** incl. the MP `hbar`/`hnode` (⇒ MP→WP staging lossless);
      node=3140, elem=5839, elem/node=1.860 (~2×). **AB_order=3 run** ⇒ **30** stores (+`urhs_AB3`/`vrhs_AB3`/
      `temp_M2`/`salt_M2`, all max|Δ|=0) — proves the conditionals. **Partition-independence:** np1 vs np2 checkpoint
      `--output-cmp` byte-value-identical (max|Δ|=0, every store incl. coords). No-regression: `run_restartsmoke.sh`
      + `run_restartcrash.sh` still GREEN np 1/2 after the `restart_register_field` refactor.

> **MP precision handling (for the Task 4.1 read-back).** `mesh%hbar`/`mesh%hnode` are `real(MP)` (`MP=max(WP,4)`;
> `==WP` at the dp anchor but the code is MP-correct for single/half builds). `t_restart_field` gains `mp_src`
> (logical) + `pmp2d(:)`/`pmp3d(:,:)` (`real(MP)` live pointers); `restart_register_field_mp` sets them and forces
> `dtype='<f8'`. `restart_write_field` copies MP→WP into a local staging buffer (`stg2`/`stg3`) before
> `decomp_redistribute` — **never binds an MP array to a WP dummy**. **Task 4.1 read-back:** for `mp_src` fields,
> read into a WP buffer then `f%pmp2d/pmp3d = real(buf, MP)` (the live MP pointer is held in the descriptor).

### Stage 4 — `mod_io_restart` READ path
*The transpose of the write path; restores owned values then halo-exchanges to bit-exact full arrays.*

#### Task 4.1 ✅: Per-field reader (`zarr_read_chunk` + `decomp_gather` + halo exchange) + safety check

**Files:**
- Modify: `src/io/mod_io_restart.F90`
- Create: `src/drivers/fesom_restartroundtrip.F90` (➕ the Task-4.1 GATE driver; mirrors `fesom_restartstate`'s
      full-state setup + adds the in-process write→corrupt→read round-trip, comparing live vs reference over the
      FULL local extent generically over the registry)
- Create: `tools/run_restartroundtrip.sh` (➕ the gate runner; `F3`/`BUILD`/`RUN` overridable for worktree runs)

- [x] `restart_read`: resolve the folder via `restart_resolve_latest` (`restart.latest` under the restart dir); read
      `checkpoint.json` (minimal line parse of `year`/`day`/`time_sec`); **aborts** (`zarr_check`→`error stop`) only on
      missing/unreadable `restart.latest`/incomplete checkpoint, or a missing/corrupt per-field store (inquires
      `<store>/.zgroup`); on `checkpoint.json` time ≠ clock **WARNs and continues** (oracle `io_restart.F90:904-914`).
      Signature `restart_read(R, restart_dir, mesh, partit, [clock_year, clock_day, clock_time_sec])`.
- [x] `restart_read_field`: writers `zarr_read_chunk` their owned chunks `[D%w_first_chunk..D%w_last_chunk]`
      (3-D loops the vertical chunks `[vc,c]` into a full-(cv,C) temp, keeping the valid rows) → `buf_writer` →
      `decomp_gather` (the transpose of the write's `decomp_redistribute`) → owned values → written **into the live
      array** (`f%p2d`/`f%p3d`; MP fields `f%pmp2d/pmp3d = real(buf,MP)`) → halo exchange via the field's recorded
      variant. Signature `restart_read_field(R, f, store_path, partit)`. Read array handle built IDENTICALLY to the
      writer (same dims/chunks/dtype/codec); 'none' **and** 'lz4' codecs round-trip.
- [x] **halo-exchange variant per field class** (recorded in `t_restart_field%halo`, set at registration; matched to
      the production consumer): `RESTART_HALO_NODE` = `exchange_nod` for **all node fields** (node halo is
      owner-consistent at end-of-step — `oce_ale.F90`); `RESTART_HALO_ELEM` = `exchange_elem` (eDim `com_elem2D`) for
      `u`/`v` (production `update_vel` `oce_ale.F90:140`; the full halo CORRUPTS the trajectory — documented there);
      `RESTART_HALO_NONE` = no exchange for `uv_rhsAB`/EVP `sigma11/12/22` (their consumers `compute_vel_rhs`/
      `stress2rhs` read owned-only, never halo-exchanged). MP node fields exchange via a WP staging copy. np=1 = no-op
      (array wholly owned). Node + element + 3-D coverage; ice incl. sigma.
- [x] **GATE:** in-process write→corrupt→read round-trip reproduces every registered array `max|Δ|=0` over the **full
      local extent** — **PASS np=1 AND np=2** (Intel dp, worktree `build_intel_dp`): `fesom_restartroundtrip` fills
      synthetic OWNED values, writes a checkpoint, halo-exchanges the originals (reference), CORRUPTS every live array
      (restart-reconstructed region → wild sentinel `-9.99e30`; un-owned tail → 0 = fresh-allocate), `restart_read`s,
      compares full local extent. AB2 = **26** fields, AB3 = **30** fields (+`urhs_AB3`/`vrhs_AB3`/`temp_M2`/`salt_M2`),
      every field `max|Δ|=0` incl. the MP `hbar`/`hnode` at np 1 AND 2; lz4 codec round-trip also `max|Δ|=0`.
      **Teeth-check:** with the read-path exchange disabled, the 21 halo-bearing fields (19 node + `u`/`v`) report
      `max|Δ|=9.99e30` in the halo/eXDim region → ROUNDTRIP FAIL (exit 1), the 5 owned-only fields stay 0 — so a wrong
      variant / gather bug is genuinely caught. No-regression: `run_restartsmoke.sh` + `run_restartstate.sh` +
      `run_restartcrash.sh` still GREEN np 1/2 after the `restart_register_field` `halo`-arg addition.

### Stage 5 — Lifecycle wiring + cadence
*Three hooks PLUS three restart-conditioning fixes the cold-start path currently hardwires. The conditioning is as
load-bearing as the I/O — each of the three, if wrong, makes the Stage 6 gate diverge or the read hook dead (all three
found by plan-review, verified against source).*

#### Task 5.1 ✅: Restart-mode detection + skip cold-start `.clock` overwrite + READ hook

**Files:**
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90`
- Modify: `src/io/mod_io_restart.F90` (no change needed — `restart_read`/`restart_resolve_latest` already exist from Stage 4)

- [x] detect restart mode **early** — before the cold-start `.clock` writer: `do_restart` = `FESOM3_RESTART=<dir>`
      set; `restart_mode` = `do_restart .and. restart_resolve_latest(RestartInPath)` returns `ok` (restart.latest ->
      a folder carrying checkpoint.json). Sets `RestartOutPath = <dir>/` (trailing-slash stripped from `<dir>` first
      so `restart_write` joins `<dir>/<folder>` cleanly while `clock_finish` keeps the slash)
- [x] **skip the unconditional `.clock` overwrite when in restart mode** — guarded `if (.not. restart_mode .and.
      mype==0)` so `clock_init` reads the previous segment's `clock_finish`-written `.clock` (two DIFFERING lines =>
      `r_restart=.true.`). GATE-CONFIRMED: seg-2 prints `RESTART MODE — ... cold .clock overwrite skipped`, the
      `.clock` shows differing lines `3600.. / 5400..`, `clock_init` sets `r_restart=.true.`
- [x] **READ hook** — after the M9 output register (`~:725`, after the whole cold dyn/tracers/ice + forcing init,
      before the loop): `if (r_restart) call restart_read(rst, RestartInPath, mesh, partit, clock_year=yearnew,
      clock_day=daynew, clock_time_sec=real(timenew,real64))` — resolves the folder via `restart.latest`,
      `decomp_gather`+halo-exchange into the live arrays (overwrites the cold state). `clock_time_sec` = sec-of-day
      (== `timenew`), matching the folder-tag/checkpoint.json convention so the time-vs-clock safety compare is clean
- [x] **GATE (FUNCTIONAL):** `tools/run_restart_lifecycle.sh` np 1 AND 2 (Intel dp, worktree build): with a checkpoint
      present seg-2 sets `r_restart=.true.` (RESTART RUN banner) and `restart_read` restores the state (`state restored
      ... at clock time=5400`); with none seg-1 is a cold INITIALISATION run (no read) — both GREEN, exit 0

#### Task 5.2 ✅: First-resumed-step AB guard (`lfirst .and. .not. r_restart`)

**Files:**
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90`

- [x] changed the step call `step_oce(n, dt, (n == 1), ...)` → `step_oce(n, dt, (n==1) .and. .not. r_restart, ...)`
      so the first resumed step blends with `ff=ab2` (AB2) over the **restored** `uv_rhsAB`, not forward-Euler
      (`ff=1.0`, `oce_dyn_velrhs.F90:135-136`) which discards it — the routine's header documents exactly this guard
      (`oce_dyn_velrhs.F90:62-64`). Cold runs keep the Euler first step (`r_restart=.false.`)
- [x] confirmed no tracer analog needed — VERIFIED in source: `lfirst` is consumed ONLY by `compute_vel_rhs`
      (`mod_step_oce.F90:206`); the tracer solve `solve_tracers_ale` (`:276/:280`) takes no `lfirst`. FESOM2's only
      `.not.r_restart` tracer guard is the cold-start `valuesold=values` init (lifecycle `:338-339`), which the
      post-cold-init READ hook overwrites (restored `<tr>_M1` => `valuesold(1,:,:)`; `valuesold(2)` is unused at the
      default tracer `AB_order==2`, restored only as `<tr>_M2` when `AB_order==3`)
- [ ] **GATE:** folded into **Task 6.1** (Stage 6 byte-gate, NOT run here) — the velocity dump on the **first resumed
      step** is where a wrong guard surfaces. The FUNCTIONAL np1/np2 2-segment gate DID exercise the guarded
      first-resumed step cleanly (lfirst=.false. => AB2 over the restored `uv_rhsAB`; resumed run completed, exit 0),
      but the `max|Δ|=0` velocity check is genuinely Task 6.1

#### Task 5.3 ✅: WRITE hooks + cadence + `clock_finish` per write

**Files:**
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90`
- Modify: `src/io/mod_io_restart.F90` (no change needed — `restart_init`/`restart_register_state`/`restart_write` already exist from Stages 3–4)

- [x] **REGISTER** (after the M9 output register `~:725`): env-gate `FESOM3_RESTART=<dir>` → `RestartOutPath`;
      `restart_init(rst, mesh, partit)` (knobs `restart_keep`/`compressor`/`n_writers`/`chunk_*` enter via its own
      `FESOM3_*` env reads) + `restart_register_state(rst, dyn, tracers, ice, mesh, mix_scheme=mix_scheme_nmb)`.
      `restart_length`/`restart_length_unit` parsed via `FESOM3_RESTART_LENGTH`/`FESOM3_RESTART_UNIT` (env fallback per
      the plan's "env fallbacks are fine if a namelist is heavy" allowance; default `1 'y'`). GATE-CONFIRMED:
      `registered 25 prognostic fields` (full oce+ice incl. EVP sigma; tke correctly ABSENT since the production
      lifecycle runs PP mixing `mix_scheme_nmb=2`, not 5 — the `ms==5` conditional fires as designed)
- [x] **WRITE** in-loop (after the M9 output-eval block, post-`step_oce`): `if (do_restart) then if (restart_due(n))
      call restart_write(rst, restart_dir, cp_year, cp_day, cp_tsec, globalstep=n)`; end-of-run (after the loop):
      a final `restart_write` guarded by `.not. wrote_final` (normally skipped — `restart_due` already treats the last
      step as due). Folder tag / `checkpoint.json` use `clock_finish`'s own year-rollover normalization
      (`daynew==ndpyr .and. timenew==86400 → 0/1/yearold+1`) so the stamps match the chained `.clock`; `time_sec` =
      sec-of-day = `timenew`
- [x] **`clock_finish` after EVERY `restart_write`** — `if (mype==0) call clock_finish()` right after each
      collective `restart_write` (periodic AND final), mirroring `io_restart.F90:573-578` (the `restart.latest` flip is
      inside `restart_write`). GATE-CONFIRMED: seg-1 wrote a PERIODIC mid-run checkpoint at step 2 + the last-step
      checkpoint at step 3, leaving `.clock` line-2 = `5400 1 1948` == the newest-checkpoint tag (resumable)
- [x] `restart_due(n)` = **global** `is_due(restart_length_unit, restart_length, n)` — ported FESOM2 `is_due`
      (`io_restart.F90:937`) + `annual_event`/`monthly_event`/`daily_event`/`hourly_event`/`step_event`
      (`gen_events.F90`) **fresh** into a contained function (the M9 `event_due` is private + coupled to
      `t_means_clock`); reads the host `restart_length`/`unit` + live `mod_clock` state; the LAST step is always due
- [x] **GATE (FUNCTIONAL):** `tools/run_restart_lifecycle.sh` np 1 AND 2 (Intel dp, worktree build): the 2-segment
      run (seg-1 cold 3 steps → 2 checkpoints + chained `.clock` + `restart.latest`; seg-2 fresh process resumes,
      detects `r_restart`, `restart_read`s, runs 3 more, writes its own checkpoints) COMPLETES exit 0; `.clock`
      chaining VERIFIED (`seg-1 .clock=5400 == checkpoint tag=5400 == seg-2 clock_init=5400`); `r_restart` fired;
      the PERIODIC mid-run checkpoint (step 2) yields a resumable `.clock` (clock_finish ran on it too). All 9
      assertions PASS at np 1 AND 2 (BYTE-exact split-vs-straight-through `max|Δ|=0` is Task 6.1, Stage 6)

### Stage 6 — Reproducibility gates

#### Task 6.1 ✅: PRIMARY self-consistency gate (straight-through vs split, np1/np2)

**Files:**
- Create: `tools/run_restart_gate_core2.sh`; relative-ε floor added to `tools/zarr_diff.py` (`output_cmp --rel-floor`)
  + `tools/dump_diff.py` (`compare(..., rel_floor) / --rel-floor=`).

- [x] straight-through: run N steps, `FESOM_DUMP_ALL` final state (node + element incl. ice + sigma)
- [x] split: run K → `restart_write` → **fresh process** → `restart_read` (`.clock`-chained) → run N−K → `FESOM_DUMP_ALL`
- [x] **two complementary comparisons** so the FULL state is covered: **#1** end-of-run CHECKPOINT compare
      (`zarr_diff --output-cmp`, ALL stores incl. ice/sigma/velocity/AB) + **#2** live `FESOM_DUMP_ALL` node dyn/tracer
      (`dump_diff --glob --ignore-step`); at np=1 and np=2.
- [x] **boundary variant:** a second split whose checkpoint lands ON the Jan→Feb month boundary (`82800 31 1948`,
      K=2 → 86400), so forcing-from-clock resume + the `roll_monthly_clim` read-ahead (Feb "slice 2") is exercised —
      asserted via the log; a mid-run integer K never hits it. PASS at np 1 AND 2.
- [x] **Two carried-state fixes found by this gate** (commit `08a4060`): `ice%thermo%t_skin` (Newton ice-surface
      solver seed) + `dyn%d_eta` (SSH-CG initial guess, converges to soltol not machine ε) — both genuine
      read-before-write cross-step prognostic state missing from Task 3.4. Localized by first-resumed-step substep
      dumps (substep 8 SSH_RHS → 9 SSH_SOLVE → clean). See memory `restart-tskin-carried-ice-state`.
- [x] **GATE:** `run_restart_gate_core2.sh` GREEN — np=1 mid+bnd **exactly `max|Δ|=0`** (whole state incl. ice+sigma);
      np=2 mid **exactly `max|Δ|=0`**; np=2 bnd **≤1 ULP** (temp `2.22e-16`, w `4e-22`) admitted by the np>1 relative-ε
      floor (1e-12) — the FESOM2-inherent redundant-element roundoff documented in the Claim-scope refinement above.
      (Intel dp, worktree build. GNU sweep folded into Task 6.3.)

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
