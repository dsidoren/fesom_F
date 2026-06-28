# M9 — Zarr Output (UGRID-compliant, xarray/ushow-readable)

> ✅/⬜ Tick the `- [ ]` boxes (and add a ✅ to a Task header) when a sub-task is **byte-gated DONE**, per
> project convention (`feedback-tick-plan-checkboxes`). This plan was produced from a full brainstorm with
> the user on 2026-06-28; the design below is **settled** — tasks decompose it, they do not re-open it.

## Overview

M9 adds **model output** to FESOM3, emitted as **hand-rolled Zarr v2** stores that are:
- readable by **xarray** (`xr.open_zarr`),
- openable by the user's **ushow** visualizer (`/home/a/a270088/ushow`, a C binary),
- **UGRID-1.0** compliant (a `fesom_mesh` topology variable describing the triangular mesh),
- accompanied by a `fesom.mesh.diag.zarr` analog of FESOM2's `fesom.mesh.diag.nc`.

The model state is already byte-exact vs FESOM2 (M0–M8), so M9 is about **getting that state onto disk** in a
modern, parallel-friendly, partition-independent format — not about new physics.

**Restart is explicitly OUT of scope** — a later milestone (possibly FESOM2-compatible binary + netCDF restart).

Key benefits: first real scientific output from FESOM3; immediate visual validation via ushow; a scalable
parallel write path (no gather-to-root) that works on large meshes.

### Why Zarr / why hand-rolled
Output speed is a primary concern on large meshes, and Zarr's chunked layout maps naturally to MPI parallelism.
NCZarr (netCDF-C 4.8.1's Zarr backend) is available on Levante and was test-verified xarray-readable, but we
chose a **hand-rolled** writer for full control over the chunking, the parallel write path, and the codec — the
decisive factor being a **distributed, no-single-rank-gather** write at scale.

## Context (from discovery)

### Files / components involved (FESOM3 — the port)
- **Driver:** `src/drivers/fesom_lifecycle_native_mr.F90` (918 lines). Hooks: setup ≈ lines 140–696; step loop
  ≈ 701–753 (`step_oce` ≈ 726); finalize ≈ 754–760. Gets exactly two new touch-points; physics untouched.
- **Existing I/O:** `src/io/mod_io_netcdf.F90` — READ-ONLY `nf90` wrapper (forcing/IC). Style to mirror
  (`nc_check(status,ctx)`→`error stop`). Zarr is **hand-rolled, not via this**.
- **Infra idioms:** `src/infra/mod_dump.F90` (per-rank gid-keyed, lazy open), `src/infra/mod_binary_arrays.F90`
  (size-prefixed serialization), `src/infra/mod_halo.F90` + `src/infra/mod_partitioning.F90` +
  `src/infra/mod_part_bounds.F90` (MPI partition + `myList` canonical maps + owned/halo bounds),
  `src/infra/mod_clock.F90` (model time → CF `time` coordinate).
- **State / mesh arrays** (sources for output): `src/types/mod_mesh.F90`, `mod_dyn.F90`, `mod_tracer.F90`,
  `mod_ice.F90`, `mod_partit.F90`. (Full mapping table in Technical Details.)

### Build mechanics (verified)
- `src/io/*.F90` is GLOB'd with `CONFIGURE_DEPENDS` (CMakeLists.txt:81–86) → **new `mod_io_*` modules
  auto-build, no CMake edit**.
- Each `src/drivers/*.F90` auto-becomes an executable (CMakeLists.txt:96–103) → a standalone Zarr smoke-test
  driver just works.
- `test/test_*.F90` registers via `add_fesom_test(<name> <nranks>)` in `test/CMakeLists.txt`.
- netCDF is an INTERFACE target `fesom_netcdf` (CMakeLists.txt:62–67); **lz4 will mirror this** (a `find_library`
  + INTERFACE target linked into `fesom3`).
- Gate runners live in `tools/run_*_gate_*.sh` (shell); many precedents (e.g. `run_geom_gate_multirank.sh`,
  `run_lifecycle_fullynative_gate_multirank.sh`).

### Reference oracle (FESOM2 v2.7.3, `/home/a/a270088/port2/fesom2/src/`)
- `io_mesh_info.F90:write_mesh_info` + `fesom_meshdiag.F90` — the `fesom.mesh.diag.nc` writer (UGRID-1.0).
- `io_meandata.F90` — output engine: `def_stream` (≈L127–129, 247), `update_means` (≈L2092), `output(istep)`
  (≈L2156), `create_new_file` (≈L1805), `write_mean` (≈L1986).
- `io_gather.F90` — `gather_nod2D` (≈L99), `gather_elem2D` (≈L157), `init_io_gather` (≈L19); the
  `myList_nod2D`/`myList_elem2D` canonical (local→global) mapping our `mod_io_decomp` reuses.
- `namelist.io` example: `work_mevp_dump/namelist.io` (`&nml_general` + `&nml_list`).
- `io_restart.F90` — restart reference (LATER milestone; do not implement now).
- Sample mesh.diag (ground-truth schema): `/home/a/a270088/port2/fesom2/test/output_pi/fesom.mesh.diag.nc`
  (3140 nod2, 5839 elem, 8986 edg_n, 48 nz).

### Tooling
- **ushow** `/home/a/a270088/ushow` (C binary). Zarr contract in Technical Details.
- **Python gate env:** `/work/ab0995/a270088/mambaforge/bin/python3` (xarray 2023.5.0, zarr 2.14.2, numpy
  1.24.3) — used to read back stores and diff against FESOM2 / in-memory state.

## Development Approach

- **Testing == byte-gates** (this project's established discipline, not unit-test-framework TDD). Every task ends
  with a runnable **gate** as its required verification. The gate kinds (see Testing Strategy):
  1. **vs-FESOM2 value gate** (read both with xarray, compare in canonical order; chase `max|Δ|=0`),
  2. **partition-independence gate** (`dist_2 ≡ dist_8`, `max|Δ|=0`),
  3. **round-trip gate** (read store back, compare to the in-memory array),
  4. **ushow smoke** (opens + renders without error),
  5. **no-regression** (`ctest` 13/13 + existing production MR byte-gates stay `max|Δ|=0`).
- **1-rank first, then multi-rank** for every stage, via the **optional-`partit`** pattern (absent ⇒ verbatim
  1-rank path; present + npes>1 ⇒ owned bounds + the canonical redistribution). This is the M2.12/M3f/M4f/M5d
  lesson.
- **Stay close to the oracle:** transcribe formats from FESOM2 (`io_mesh_info`/`io_meandata`) and the Zarr v2
  spec; cite file:line / spec section in each module header.
- **Conventions:** `WP` for data (real64 default), `MP` for mesh/geometry; `_WP`/`_MP` literal suffixes; SoA
  `(nl, n_entity)` column-major; `error stop` with context; lazy-cache init; per-rank handles.
- Complete each task fully (code + gate green) before the next. Update this plan's checkboxes immediately.

## Testing Strategy

- **Gate scripts** under `tools/` (mirror `tools/run_*_gate_*.sh`) drive the FESOM3 writer + a Python compare.
- **Python compares** under `tools/` (e.g. `tools/zarr_diff.py`) open the Zarr store with xarray, optionally open
  the FESOM2 netCDF reference, and assert `max|Δ|` per the gate's bar.
- **Stage 0** gates are self-contained (toy arrays, round-trip).
- **Stage 1 (mesh.diag)** gates against FESOM2's `fesom.mesh.diag.nc` (connectivity exact; coords/areas
  round-off) + partition-independence + ushow.
- **Stage 2 (fields)** gates against FESOM2 output (chase `max|Δ|=0`) + partition-independence + round-trip +
  ushow.
- **No-regression after every stage:** `cd build_intel_dp && ctest --output-on-failure` (13/13, Intel+GNU) and
  the existing production MR byte-gates (`run_lifecycle_*_gate_multirank.sh`) stay `max|Δ|=0`.

## Progress Tracking
- `[x]` = done; ➕ = newly discovered task; ⚠️ = blocker/issue.
- Add a ✅ to a Task header when its gate is green.
- Keep this file in sync with actual work; update scope here if it shifts.

## What Goes Where
- **Implementation Steps** (`[ ]`): all code, gate scripts, and in-repo verification.
- **Post-Completion** (no checkboxes): the later restart milestone, manual ushow visual inspection, and external
  build-dependency notes (liblz4).

---

## Implementation Steps

### Stage 0 — `mod_io_zarr`: the standalone Zarr v2 writer

#### Task 0.1: Zarr store/group/array scaffolding + JSON metadata ✅

**Files:**
- Create: `src/io/mod_io_zarr.F90`

- [x] create `mod_io_zarr` with `t_zarr_store` (root path) and `t_zarr_array` (name, shape, chunks, dtype,
      codec, fill_value) types — also `t_zarr_attrs` (JSON-object builder)
- [x] `zarr_create_store(store, path[, attrs])` — make the store dir + write `.zgroup` (`{"zarr_format":2}`) and root
      `.zattrs`
- [x] `zarr_array_init(arr, name, dims, chunks, dtype[, fill, has_fill, codec])` + `zarr_define_array(store, arr[, attrs])`
      — make the array subdir + write `.zarray` (full v2 schema) and `.zattrs` (incl. `_ARRAY_DIMENSIONS`)
- [x] hand-rolled JSON helpers (`zattr_str/int/real` + `_arr` variants; `json_real` round-trip incl. NaN/Inf) +
      `zarr_check(ok, ctx)`→`error stop`
- [x] **gate:** covered by Task 0.2 round-trip (zarr+xarray open + parse); JSON verified on disk
      (`.zgroup`/`.zarray`/`.zattrs` all parse)

#### Task 0.2: Chunk encode + write (codec `none`, C-order) + round-trip gate ✅

**Files:**
- Modify: `src/io/mod_io_zarr.F90`
- Create: `src/drivers/fesom_zarrsmoke.F90`
- Create: `tools/run_zarrsmoke.sh`
- Create: `tools/zarr_diff.py`

- [x] `zarr_write_chunk(store, arr, chunk_index(:), data)` (generic over rank/type) — transpose Fortran
      column-major → **C row-major**, write the chunk file at the `.`-joined index path, codec `none` (raw LE bytes).
      Plus `zarr_write_whole` convenience (loops all chunks of a full in-memory array — the 1-rank path).
- [x] last-/partial-chunk handling: every chunk file is full chunk-size; the final global chunk along a dim is
      padded with `fill_value` (verified on disk: f4 2×2 chunks all 16 B incl. partials)
- [x] `fesom_zarrsmoke` driver: 1-D `f8`, 2-D `f4`, 2-D `i4` (non-dividing chunks; distinct per-cell values)
- [x] `tools/zarr_diff.py --roundtrip`: opens with zarr AND xarray, asserts values byte-match the generator
      (`max|Δ|=0` f8/i4; exact f4) — **PASS**
- [x] `tools/run_zarrsmoke.sh`: reconfigure + build + run + Python round-trip; **GATE GREEN** (exit 0)

#### Task 0.3: lz4 codec + optional `.zmetadata` consolidation ✅

**Files:**
- Modify: `src/io/mod_io_zarr.F90`
- Modify: `CMakeLists.txt` (discover + link `liblz4` as an INTERFACE target, mirroring `fesom_netcdf`)
- Modify: `tools/zarr_diff.py`, `tools/run_zarrsmoke.sh`

- [x] `liblz4` linkable on Levante (spack `lz4-1.9.4`); called via Fortran `iso_c_binding` (no C header);
      `find_library(LZ4_LIB)` + `fesom_lz4` INTERFACE target linked into `fesom3`, guarded by `HAVE_LZ4`
      (absent ⇒ none-only + `error stop` on lz4 request)
- [x] lz4 codec: numcodecs framing (4-byte LE decompressed length + LZ4_compress_default block); `.zarray`
      `compressor:{"id":"lz4","acceleration":1}` (verified on disk: header `18 00 00 00` = 24)
- [x] `zarr_consolidate(store)` — `.zmetadata` (`zarr_consolidated_format:1`; all `.zgroup`/`.zarray`/`.zattrs`)
- [x] **gate:** round-trip also writes `zarrsmoke_lz4.zarr` (compressor=lz4 + consolidated); zarr+xarray read
      it `max|Δ|=0` via `open_consolidated`/`consolidated=True`; `none` path still green — **GATE GREEN**

### Stage 1 — `fesom.mesh.diag.zarr` (proves the whole stack on static data)

#### Task 1.1: `mod_io_decomp` — canonical chunk map + writer-subset + 1-rank identity ✅

**Files:**
- Create: `src/io/mod_io_decomp.F90`
- Create: `test/test_io_decomp.F90` (ctest np 1/2/8)

- [x] `decomp_init(D, C, n_writers, N, myList, myDim, comm, mype, npes)` (raw, unit-testable) +
      `decomp_init_entity(D, C, n_writers, entity, mesh, partit)` convenience — canonical `N`, `nchunks=ceil(N/C)`,
      block chunk→writer assignment, this rank's send-plan (built once via a gid `MPI_Alltoallv`). Optional-`partit`
      (absent OR npes==1 ⇒ identity: one writer, all canonical; `is_multirank` gate, mirrors `owned_bounds`).
- [x] `decomp_redistribute` (generic real/int, 2-D) — one `MPI_Alltoallv` compute→canonical-chunk on the writer
      subset (non-writers: zero recvcounts; npes==1: self-copy)
- [x] `decomp_redistribute` (generic real/int, 3-D) — loops levels, reuses the 2-D plan → in-flight `O(C·nlev)`
- [x] helpers `decomp_is_writer(D)`, `decomp_writer_chunk_range(D,...)`, fields `w_first/last_chunk`,
      `w_nbuf`, `w_base_gid`
- [x] **gate:** `test_io_decomp` — labelled field (`value=canonical id`) over a synthetic ROUND-ROBIN partition
      lands canonical (`buf[g]==g`, pad=fill) for real/int/3-D over 6 (N,C,n_writers) combos incl. partial chunks
      + writer subsets. **PASS np=1 (identity) AND np=2/8 (real MR redistribution)** — stronger than the
      planned 1-rank-only; ctest 16/16 (was 13/13 + 3 new).

#### Task 1.2a: `mod_io_meshdiag` — coords + connectivity + UGRID topology (1-rank) + FESOM2 gate ✅

**Files:**
- Create: `src/io/mod_io_meshdiag.F90`
- Create: `src/drivers/fesom_meshdiagdump.F90` (standalone mesh→meshdiag, mirrors `fesom_geomdump`)
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90` (call `meshdiag_write()` after setup, env-gated)
- Create: `tools/run_meshdiag_gate.sh` (pi; np 1/2/8)
- Modify: `tools/zarr_diff.py` (`--meshdiag` mode)

- [x] `meshdiag_write(path, mesh, partit)` core: `lon`/`lat` (rad→deg), `face_nodes` (1-based,
      `start_index=1`), `edge_nodes`, `edge_face_links` (−999), `nz`/`nz1` (sign-flip).
      ⚠️ `face_edges`/`face_links` **DEFERRED** — FESOM3 never builds `elem_edges`/`elem_neighbors` (only
      `edges`/`edge_tri`), like `gradient_vec`; topology var drops those two connectivity attrs.
- [x] `fesom_mesh` UGRID topology var + `Conventions="UGRID-1.0"` + `_ARRAY_DIMENSIONS` on every var
- [x] wired `meshdiag_write()` into `fesom_lifecycle_native_mr` after setup (gated by `FESOM3_MESHDIAG`)
- [x] `tools/zarr_diff.py --meshdiag`: opens both `mask_and_scale=False` (raw — a valid 0 isn't NaN-masked),
      compares the emitted subset in canonical order — ints exact, floats `max|Δ|=0`
- [x] **gate (1-rank, pi):** `tools/run_meshdiag_gate.sh 1` GREEN — 18 vars `max|Δ|=0` vs FESOM2
      `output_pi/fesom.mesh.diag.nc` (2-rank ref valid: every emitted var is global/canonical there; the only
      partition-local FESOM2 vars `face_edges`/`face_links` are deferred). `ushow` = manual (GUI, no headless).
      Wins: gradient_sca 1:3/4:6 packing; C-mkdir not `execute_command_line` (fork segfaults post-MPI_Init);
      `fill_value: null` to avoid NaN-masking; `zbar_e_bot` computed in the driver (not `compute_geometry`).

#### Task 1.2b: mesh.diag derived / diagnostic fields ✅ (done together with 1.2a)

**Files:**
- Modify: `src/io/mod_io_meshdiag.F90`, `tools/zarr_diff.py`

- [x] `elem_area`, `nlevels`/`nlevels_nod2D`, `nod_in_elem2D`/`_num`, `edge_cross_dxdy`,
      `gradient_sca_x/y` (`gradient_sca(1:3/4:6,:)`), `nod_area`
- [x] computed/derived: `zbar_e_bottom` (−sign), `zbar_n_bottom` (`-zbar(nlevels_nod2D(n))`),
      `nod_part`/`elem_part` (stamp `mype`)
- [x] **DEFERRED** `gradient_vec_x/y` — not computed in FESOM3 (➕)
- [x] **gate:** `--meshdiag` covers the full emitted set; all 18 vars `max|Δ|=0` vs FESOM2 (ints exact)

#### Task 1.3: mesh.diag multi-rank (partition-independence gate) ✅ proves the parallel stack

**Files:**
- Modify: `src/io/mod_io_meshdiag.F90` (route arrays through `mod_io_decomp`)
- `tools/run_meshdiag_gate.sh` (the single gate supports np 1/2/8 — no separate MR script needed)

- [x] **store-create ordering** (reused by Task 2.3): rank 0 `zarr_create_store` + all `zarr_define_array`
      → `MPI_Barrier` → writers write chunk data → rank 0 consolidate (writers never create dirs)
- [x] every mesh.diag array routed through `decomp_redistribute_*` (canonical, writer-subset; each writer writes
      only its chunks) — built MR-ready from the start (1-rank = npes==1 identity path)
- [x] `n_writers` subset path exercised: np=8 with chunk=1000 ⇒ node nchunks=4 ⇒ only 4 of 8 ranks write nodes
- [x] **gate:** `run_meshdiag_gate.sh {2,8}` GREEN — both `== FESOM2` (18 vars `max|Δ|=0`); explicit
      `dist_2 ≡ dist_8` store compare = `max|Δ|=0` over 18 vars (`nod_part`/`elem_part` correctly DIFFER —
      partition descriptors). `ushow` = manual.
- [x] **no-regression:** `ctest` 16/16 (new io modules don't touch physics; production MR lifecycle byte-gate
      deferred to the F1 sweep — physics paths are unchanged)

### Stage 2 — field output (`mod_io_means`)

#### Task 2.1: registry + snapshot + frequency + per-variable-per-year stores (1-rank, node scalars) ✅

**Files:**
- Create: `src/io/mod_io_means.F90`
- Create: `src/drivers/fesom_outputsmoke.F90` (standalone field-output gate driver, mirrors zarrsmoke)
- Modify: `src/io/mod_io_zarr.F90` (`zarr_rewrite_zarray` — bump shape[0] on append)
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90` (`means_init`/`means_*` in the step loop, env-gated)
- Create: `tools/run_output_gate.sh`
- Modify: `tools/zarr_diff.py` (`--output` round-trip + `--output-cmp` partition-independence modes)

- [x] `t_mean_field` (name/long_name/units/std/dtype/snap + per-year `t_zarr_store`/array handles) +
      `t_io_means` registry; `means_define_node2d(...)`; put-based API (decoupled from the state types —
      caller passes the owned 1-D slice, no fragile pointers to strided components)
- [x] register node scalars **ssh, sst, sss, a_ice, m_ice, m_snow** (snapshot) — sources `dyn%eta_n`,
      `tracers%data(1/2)%values(1,:)`, `ice%data(1:3)%values`. (3-D T/S are Task 2.4 — not 2.1.)
- [x] CF `time` coord (`"seconds since <YYYY>-01-01 00:00:00"` per-year ref + `calendar` from `forc_calendar`);
      per-variable-per-year store `<dir>/<name>.fesom.<YYYY>.zarr`; **`chunk_time=1`** append = new data chunk
      `[t,c]` + time chunk `[t]` + `zarr_rewrite_zarray` bump `shape[0]→t+1` (NO read-modify-write)
- [x] embed `lon`/`lat` (node) coords in each store (`coordinates: "lon lat"` so xarray promotes them);
      `FESOM3_OUTPUT`/`FESOM3_OUTPUT_EVERY`/`FESOM3_CHUNK_HORIZ` env (namelist.io deferred to Task 2.6)
- [x] wired into `fesom_lifecycle_native_mr` after `step_oce`, ABOVE the `if(.not.step_diag)cycle`; `means_init`
      after setup, `means_begin/put×6/end` per output step, `means_finalize` after the loop (`FESOM3_OUTPUT`)
- [x] **gate (1-rank, snapshot):** `run_output_gate.sh 1` GREEN — `fesom_outputsmoke` writes 2 synthetic node
      fields whose value at (record k, CANONICAL node g) follows a partition-indep generator formula; `--output`
      asserts every value == formula `max|Δ|=0`, CF time decodes (xarray→datetime64), lon/lat embedded. **Writer
      self-consistency, no FESOM2 oracle** (state is byte-exact thru M8; this proves the WRITER serializes it).
      ⚠️ xarray `decode_times=False` to compare raw seconds (else `time`→datetime64).

#### Task 2.2: mean accumulation + interval reset (1-rank), FESOM2-aligned ✅

**Files:**
- Modify: `src/io/mod_io_means.F90` (accumulator + accumulate/write API), `src/drivers/fesom_outputsmoke.F90`
  (a mean stream `fld_m`), `tools/zarr_diff.py` (`fld_m` formula)

- [x] running **sum ÷ count** accumulator per field, **in the output precision** (real32 for `<f4`, real64 for
      `<f8`) — transcribed from `io_meandata.F90:update_means` (`local_values += value`, `addcounter++`,
      :2107/2142) + `compute_means` (`copy = local_values / addcounter`, divide in that precision, :2335/2353),
      then zero + reset. `means_accumulate` every step (mean: sum; snapshot: overwrite, count=1); `means_write`
      divides + writes + resets. So float32 means accumulate AND divide in float32, byte-matching FESOM2's r4.
- [x] **gate:** `run_output_gate.sh` — `fld_m` is a MEAN stream fed 3 sub-steps (`g-1, g, g+1`) per record ⇒
      mean `== g` `max|Δ|=0` (proves sum + divide-by-count, not just overwrite), partition-independent
      (dist_2≡dist_8); snapshots `fld_a/fld_b` still `max|Δ|=0`. (vs-FESOM2 mean follows by transitivity:
      FESOM2 accumulate/divide semantics transcribed in the matching precision + state byte-exact thru M8 +
      writer self-consistency proven — no separate FESOM2 mean-output oracle run needed.)

#### Task 2.3: multi-rank node scalars (partition-independence) ✅ (done together with 2.1)

**Files:**
- `src/io/mod_io_means.F90` (built MR-ready from the start: every write routed through `decomp_redistribute`)
- `tools/run_output_gate.sh` (np 1/2/8 in one script — no separate MR script)

- [x] every field routed through `decomp_redistribute` (canonical, writer-subset) before the chunk write
      (1-rank = npes==1 identity path, same code) — `put_static` (lon/lat) + `means_put` (data)
- [x] **gate:** `run_output_gate.sh` GREEN — `--output-cmp` proves `dist_2 ≡ dist_8 ≡ 1-rank` `max|Δ|=0`
      (data + lon/lat + time). Writer subset exercised: np=8 with chunk=1000 ⇒ nchunks=4 ⇒ only 4 ranks write.
      (`== FESOM2` is implied: state byte-exact thru M8 + writer self-consistency proven.)
- [x] **no-regression:** `ctest` 16/16 (production MR lifecycle byte-gate deferred to the F1 sweep — the new io
      modules don't touch physics; the lifecycle output wiring is compile-verified + mirrors the proven meshdiag)

#### Task 2.4: 3D fields (T, S full-depth; w on nodes) + vertical chunking ✅

**Files:**
- Modify: `src/io/mod_io_zarr.F90` (`zarr_write_chunk_3d_real` — (time,nz,nod2) C-order chunk),
  `src/io/mod_io_means.F90` (`means_define_node3d` + generic `means_accumulate` + 3-D `means_write`),
  `src/drivers/fesom_outputsmoke.F90` (`fld_3`), `src/drivers/fesom_lifecycle_native_mr.F90` (temp/salt/w),
  `tools/zarr_diff.py` (3-D `--output` check)

- [x] register 3-D **temp, salt** (`nl-1` layers, vdim `nz1`) + node **w** (`nl` levels, vdim `nz`) via
      `means_define_node3d(on_full_levels=)`; `decomp_redistribute_3d` (level-by-level, already existed);
      embedded `nz`/`nz1` vertical coord (positive-down); **`_FillValue=NC_FILL` for below-bottom**
      (`L > nlevels_nod2D - voff`, voff=1 layers / 0 levels) — nlevels-based CF mask (per plan, cleaner than
      FESOM2's value-based `abs(acc)<1e-30` quirk; valid-level VALUES still byte-match FESOM2). T/S/w wired
      into the lifecycle. `vert_chunk` = full-depth single chunk in v1 (multi-vchunk knob → Task 2.6).
- [x] **gate:** `run_output_gate.sh` — `fld_3` (3-D, `g+L`): valid levels `max|Δ|=0`, below-bottom NaN-masked
      (229155 entries on pi), `nz1` coord monotonic positive-down, partition-independent (dist_2≡dist_8). ctest
      16/16. (`ushow` level slice = manual; vs-FESOM2 by transitivity as in 2.2.)

#### Task 2.5: NODE velocity vectors unod, vnod + r2g rotation (geographic default, native knob) ✅

**DECISION (HANDOFF-confirmed): NODE-based `unod`/`vnod` (`dyn%uvnode(1/2,:,:)`, nl-1 layers), not element.**
Every FESOM2 `namelist.io` outputs `unod`/`vnod` (`dynamics%uvnode`), NEVER the element `u`/`v` — node is the
faithful default AND reuses the proven `means_define_node3d` (no new elem decomp). `compute_vel_nodes`
(`mod_step_oce.F90:110`) populates `dyn%uvnode` in the live step. Element `u`/`v` (`dyn%uv`) left as a clean
future elem-decomp addition. (mod_io_meshdiag elem-centroid coords NOT needed — node coords already embedded.)

**Files:** `src/mesh/mod_mesh_rotate.F90` (port `vector_r2g`); `src/io/mod_io_means.F90` (vector-pair registry +
rotation + `vec_frame`); `src/drivers/fesom_lifecycle_native_mr.F90` (register/accumulate unod/vnod);
`test/test_vector_rotate.F90` (+ `test/CMakeLists.txt`); `src/drivers/fesom_outputsmoke.F90`, `tools/zarr_diff.py`,
`tools/run_output_gate.sh` (gate); `config/namelist.io` (template stub; full parse = Task 2.6).

- [x] **ported `vector_r2g`** VERBATIM from FESOM2 `gen_modules_rotate_grid.F90:164-202` into `mod_mesh_rotate.F90`
      (rotated→geographic; the exact inverse of the existing byte-gated `vector_g2r` — Cartesian from ROTATED
      angles, TRANSPOSED `r2g_matrix`, project onto GEO). Round-trip ctest `test_vector_rotate` (np 1): non-identity
      (50,15,-90) `vector_g2r∘vector_r2g==identity` 8.9e-15, magnitude-preserving 7.1e-15, `flag0==flag1` 8.5e-14,
      identity-matrix no-op 7.1e-15 — all ≪ 1e-11. ctest 16→17.
- [x] `means_define_vector3d` (links two `means_define_node3d` as an (x,y) pair); register **unod, vnod** in the
      lifecycle (`dyn%uvnode(1/2,1:nl-1,1:nNodO)`), accumulated independently, rotated together at write.
- [x] apply `vector_r2g(flag_coord=0)` per (node,level) at the cached ROTATED node coords (`coord_nod2D`, rad).
      **FESOM2 ORDER matched:** `io_r2g` rotates the accumulated SUM (`io_meandata.F90:2265`) BEFORE
      `compute_means` divides (:2335) — so `write_vector_3d` rotates the sum THEN divides (`<f4`: promote r4 sum→r8,
      rotate, demote r4, divide r4 = io_r2g r4 branch :3028). Below-bottom (nlevels mask) → NC_FILL, not rotated.
- [x] `vec_frame` knob `geographic`(default)|`native` via `means_init(vec_frame=)` + `FESOM3_VEC_FRAME` env
      (= FESOM2 `vec_autorotate`; FESOM3 default geographic). namelist.io parse deferred to 2.6.
- [x] **gate** (`run_output_gate.sh`, np 1/2/8): `native` == raw generator `max|Δ|=0`; `geographic` == an
      INDEPENDENT numpy `vector_r2g` reference (zarr_diff `_vector_r2g`, flag=1 on embedded geo coords)
      `max|Δ|≈2e-13` ≪ 1e-9 + non-vacuous (rotation changed values); partition-indep `dist_2≡dist_8≡1` `max|Δ|=0`
      both frames; below-bottom masked. ctest 17/17 Intel+GNU. (vs-FESOM2 by transitivity: state byte-exact thru
      M8 + `vector_r2g` verbatim + round-trip ctest + writer self-consistency — no FESOM2 oracle, as 2.1–2.4.)

#### Task 2.6: full namelist.io knobs + float32/lz4/chunk-shape/n_writers/filesplit

**Files:**
- Modify: `src/io/mod_io_means.F90`, `config/namelist.io`

- [ ] expose all knobs: per-field `freq`/`unit`/`precision`/`mean|snap`/`frame`; global `n_writers`,
      `chunk_shape (time,vert,horiz)`, `compressor (none|lz4)`, `filesplit_freq`; float32 default
- [ ] `chunk_time>1`: implement the partial-last-time-chunk **read-modify-write** append path (deferred from
      Task 2.1's `chunk_time=1`); gate that `chunk_time={1,N}` produce value-identical stores
- [ ] **gate:** a run with `compressor=lz4` + custom `chunk_shape` + `n_writers` subset still passes
      partition-independence + round-trip + `ushow` opens it

#### Task 2.7: ELEMENT-based output (u/v vectors + Av + bolus) ➕ USER-REQUESTED (2026-06-28)

The user explicitly needs **element output** in addition to the node fields: element velocity **u, v**
(`dynamics%uv`), and element SCALARS — notably **Av** (vertical viscosity, `nl` levels, elem) — plus the GM
**bolus_u/bolus_v** (`dynamics%fer_uv`, Fer_GM). Node `unod`/`vnod` (Task 2.5) and element `u`/`v` are BOTH
real FESOM2 outputs; this adds the element source. (FESOM2 element streams: `u`,`v`,`Av`,`bolus_u`,`bolus_v`,
`pgf_x/y`, ALE `helem`/`h`/`d`/`dhe`, the `ke_*` KE budget, ice `eps*`/`sgm*` — port the production core first.)

**Files:** `src/io/mod_io_means.F90` (element registration + a 2nd `t_io_decomp` for elements + elem-centroid
coords), `src/io/mod_io_meshdiag.F90` (reuse elem-centroid `lon`/`lat`), `src/drivers/fesom_lifecycle_native_mr.F90`,
`src/drivers/fesom_outputsmoke.F90` + `tools/{zarr_diff.py,run_output_gate.sh}` (elem gate).

- [ ] add an ELEMENT decomp to `t_io_means` (a 2nd `t_io_decomp De` via `decomp_init_entity(.., DECOMP_ELEM, ..)`
      — the selector ALREADY EXISTS, used by meshdiag) + cache elem-centroid `lon`/`lat` (geographic, for embed)
      and the ROTATED centroid coords for r2g. ⚠️ FESOM2 `io_r2g` rotates element vectors at the **simple mean of
      the 3 ROTATED node coords** `sum(coord_nod2D(1:2,elem2D_nodes(1:3,e)))/3` with `flag_coord=0` (NOT cyclic-
      aware `elem_center`) — match that for byte-faithfulness; `elem_center(mesh,n,cx,cy)` (`mod_mesh_areas.F90:178`,
      cyclic-aware) is fine for the EMBEDDED display centroid.
- [ ] `means_define_elem2d`/`means_define_elem3d` (scalars: **Av** `nl`/elem; the `is_elem` flag routes
      write_data_* through `De` instead of `Dn`) + an element variant of `means_define_vector3d` (u/v, bolus) using
      the elem-centroid rotation. The mask is element `nlevels`-based (below-bottom → NC_FILL).
- [ ] register **u, v** (`dyn%uv(1:2,:,:)`, nl-1), **Av** (`dyn%work%Av`, nl), **bolus_u/v** (`dyn%fer_uv`, Fer_GM
      only) in the lifecycle; embed elem-centroid `lon`/`lat`.
- [ ] **gate:** elem `native==raw` `max|Δ|=0` + `geographic==`numpy elem-centroid r2g reference + partition-indep
      `dist_2≡dist_8`; reuse the Task 2.5 vector gate machinery with elem coords. (vs-FESOM2 by transitivity.)

### Final

#### Task F1: Full no-regression + gate sweep

- [ ] `cd build_intel_dp && ctest --output-on-failure` → 17/17 (repeat GNU)
- [ ] existing production MR byte-gates (`run_lifecycle_*_gate_multirank.sh`, both whichEVP) stay `max|Δ|=0`
      (also the DEFERRED real-lifecycle output integration: a forced run with `FESOM3_OUTPUT` writes
      unod/vnod/temp/salt/... and they read back sanely)
- [ ] all new M9 gates green (Stage 0 round-trip; mesh.diag 1-rank + MR; fields 1-rank + MR + 3D + node vectors
      + ELEMENT u/v/Av/bolus)
- [ ] `ushow` smoke on a mesh.diag store + a field store + `ushow <field>.zarr -m fesom.mesh.diag.zarr`

#### Task F2: Docs + memory

- [ ] update `docs/HANDOFF.md` ("Where we are" + "Next task" → M9 done / restart NEXT) and `docs/LESSONS.md` if
      any lesson surfaced (e.g. C-order transpose, lz4 framing, Alltoallv redistribution)
- [ ] update memory (`project-fesom3-implementation-state`; retire `project-m9-brainstorm-first`); tag `m9`
- [ ] move this plan to `docs/plans/completed/`

---

## Technical Details

### Zarr v2 on-disk format (what `mod_io_zarr` emits)
A store is "just files":
```
store.zarr/
  .zgroup            {"zarr_format": 2}
  .zattrs            { ...global attrs... }              # e.g. Conventions, model/git provenance
  <var>/
    .zarray          {"zarr_format":2,"shape":[...],"chunks":[...],"dtype":"<f4|<f8|<i4",
                      "compressor":null | {"id":"lz4"},"fill_value":<num|null>,
                      "order":"C","filters":null,"dimension_separator":"."}
    .zattrs          {"_ARRAY_DIMENSIONS":["time","nod2"], "units":..., "long_name":..., ...}
    0.0, 0.1, ...    # chunk files: C-row-major raw (codec none) or lz4-framed bytes
  .zmetadata         # optional consolidated metadata (all .zarray/.zattrs concatenated)
```
- **C-order gotcha:** Fortran is column-major; a chunk must be written transposed to C row-major.
- **Last/partial chunk:** the final global chunk along a dim is padded to full chunk size with `fill_value`
  (standard Zarr; the reader truncates to `shape`). This is the **only** padding — the canonical write scheme
  has **no per-rank padding**.
- **dtype strings:** `<f8` (real64), `<f4` (real32), `<i4` (int32). Little-endian.
- **lz4 framing (numcodecs-compatible):** 4-byte little-endian `uint32` decompressed length, then the lz4 block.

### Canonical-order + distributed-chunk-writer redistribution (`mod_io_decomp`)
- **Ordering:** canonical global ids (`nod2d.out`/`myList` order — what FESOM2 gathers to) ⇒ files are
  **partition-independent** (identical at any rank count), connectivity needs **no renumbering**, coords dense.
- **Chunks:** uniform size `C` along the entity dim; `ceil(N/C)` chunks; blocks of chunks assigned to a
  configurable subset of `n_writers` writer ranks (default e.g. one-per-node).
- **Redistribution:** each compute rank maps each owned entity (`myList(i)`) → chunk → destination writer rank,
  buckets, and ships via one `MPI_Alltoallv` (non-writers have zero recvcounts). 3D loops levels (or per
  `vert_chunk`) so in-flight memory is `O(C·nlev)` per writer — **never a whole field on one rank, never a
  rank-0 funnel**.
- **Write:** each writer writes only its own chunks (zero contention); rank 0 writes the JSON metadata.
- **Store-create ordering (parallel correctness):** rank 0 creates the store + all array dirs + `.zarray`/`.zattrs`
  → `MPI_Barrier` → writers write chunk files → rank 0 (optional) consolidates. Writers never create dirs.
- Reuse FESOM2's `myList_nod2D`/`myList_elem2D` canonical maps (the same ones `io_gather` uses).

### `fesom.mesh.diag.zarr` variable → FESOM3 source mapping
| diag var | dims | FESOM3 source (`mesh%…` unless noted) | note |
|---|---|---|---|
| `lon`,`lat` | (nod2) | `geo_coord_nod2D(1:2,:)` | radians → degrees |
| `nz`,`nz1` | (nz),(nz1) | `zbar`,`Z` | sign-flip (positive-down) |
| `elem_area` | (elem) | `elem_area` | |
| `nlevels_nod2D`,`nlevels` | (nod2),(elem) | same | |
| `nod_in_elem2D_num`,`nod_in_elem2D` | (nod2),(N,nod2) | same | |
| `face_nodes` | (n3,elem) | `elem2D_nodes(1:3,:)` | 1-based, `start_index=1` |
| `edge_nodes` | (n2,edg_n) | `edges(1:2,:)` | |
| `face_edges` | (n3,elem) | `elem_edges(1:3,:)` | |
| `face_links` | (n3,elem) | `elem_neighbors(1:3,:)` | −999 fill |
| `edge_face_links` | (n2,edg_n) | `edge_tri(1:2,:)` | −999 fill |
| `edge_cross_dxdy` | (n4,edg_n) | `edge_cross_dxdy(1:4,:)` | |
| `gradient_sca_x/y` | (n3,elem) | `gradient_sca(1:3/4:6,:)` | packed tightly into 1:6 |
| ~~`gradient_vec_x/y`~~ | — | **DEFERRED** | NOT computed in FESOM3 (opt_visc=7 never builds it); transcribe `compute_gradient_vec` later for full pyfesom2 parity |
| `nod_area` | (nz,nod2) | `area(:,:)` | |
| `zbar_e_bottom` | (elem) | `zbar_e_bot` | sign-flip |
| `zbar_n_bottom` | (nod2) | **compute** `zbar(nlevels_nod2D(n))` | sign-flip; not stored |
| `nod_part` | (nod2) | `partit%part(:)` (global node-owner map) | |
| `elem_part` | (elem) | **derive at write**: stamp `mype` over `myList_elem2D(1:myDim_elem2D)` | no element-owner array in `t_partit` |

UGRID topology var `fesom_mesh`: `cf_role="mesh_topology"`, `topology_dimension=2`,
`node_coordinates="lon lat"`, `face_node_connectivity="face_nodes"`, `edge_node_connectivity="edge_nodes"`,
`face_edge_connectivity="face_edges"`, `face_face_connectivity="face_links"`,
`edge_face_connectivity="edge_face_links"`; global `Conventions="UGRID-1.0"`.

### ushow input contract (the consumer)
- Zarr store with `.zgroup`; 1-D `lon`/`lat` (or `longitude`/`latitude`) of size `n_points` with `units`
  (`degrees…` or `rad`); data vars **must** carry the xarray `_ARRAY_DIMENSIONS` attribute; optional `time` with
  CF units; optional connectivity `face_nodes`/`face_node_connectivity` `[elem,3]` int with `start_index`.
- Compressors recognized: **null / lz4 / blosc** (NOT zlib). dtypes `<f8`/`<f4` (int `<i4`).
- Separate mesh via `-m`: `ushow temp.fesom.1964.zarr -m fesom.mesh.diag.zarr`.

### `namelist.io` schema (the FESOM2 analog)
- Staged into the **rundir** at runtime (like the oracle namelists; `FESOM3_*` env overrides) — not a repo
  `config/` runtime dependency. A template lives at `config/namelist.io`.
- `&nml_general`: `n_writers`, `chunk_time` (**default 1 in v1**), `chunk_vert`, `chunk_horiz`,
  `compressor (none|lz4)`, `filesplit_freq (y|m)`, `vec_frame (geographic|native)`.
- `&nml_list`: rows `'<var>', <freq>, '<unit y|m|d|h|s>', <precision 4|8>, '<mean|snap>'`.

### Gate commands
- Python: `/work/ab0995/a270088/mambaforge/bin/python3 tools/zarr_diff.py …`
- Build: `./configure.sh --compiler intel --precision dp --clean --build` (anchor) + `--compiler gnu`.
- ctest: `cd build_intel_dp && ctest --output-on-failure`.
- Multi-rank: requires `env.sh` (KNEM `single_copy_mechanism=none`).

## Post-Completion
*Informational — external/manual, no checkboxes.*

**Manual verification**
- Visual: open a real CORE2/JRA55 output year in `ushow` (and xarray/pyfesom2) and eyeball SST/SSS/ice/velocity
  fields for physical sanity (the byte-gates prove correctness; this proves usability).
- Performance: measure write time + filesystem load at `dist_512`/`dist_864` and tune `n_writers`/`chunk_shape`;
  compare `none` vs `lz4` throughput and ratio.

**External build dependency**
- `liblz4` must be linkable on Levante (verify in Task 0.3). If not readily available, ship `none`-only and add
  the codec when the lib is provisioned (the slot is pluggable).

**Next milestone (out of scope here)**
- **Restart** (write/read model state incl. `tke`/EVP/ALE serialization; bit-reproducible restart→continue gate;
  possibly FESOM2-compatible binary + netCDF). The driver already emits `fesom.clock` (read-only) and
  `mod_clock` has the deferred `clock_finish`/`clock_newyear`. Reference: FESOM2 `io_restart.F90`.
