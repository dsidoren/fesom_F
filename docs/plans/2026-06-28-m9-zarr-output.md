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

#### Task 0.1: Zarr store/group/array scaffolding + JSON metadata

**Files:**
- Create: `src/io/mod_io_zarr.F90`

- [ ] create `mod_io_zarr` with `t_zarr_store` (root path) and `t_zarr_array` (name, shape, chunks, dtype,
      codec, fill_value) types
- [ ] `zarr_create_store(store, path)` — make the store dir + write `.zgroup` (`{"zarr_format":2}`) and root
      `.zattrs`
- [ ] `zarr_define_array(store, arr, name, shape, chunks, dtype, attrs, fill_value, codec)` — make the array
      subdir + write `.zarray` (full v2 schema, see Technical Details) and `.zattrs` (incl. `_ARRAY_DIMENSIONS`)
- [ ] hand-rolled JSON helpers (ints/reals/strings/1-D arrays of each) + `zarr_check(ok, ctx)`→`error stop`
- [ ] **gate:** the Stage-0 smoke driver (Task 0.2) asserts `.zgroup`/`.zarray`/`.zattrs` exist and parse via
      Python `json` + `zarr` (deferred to 0.2 where data is written) — for now, a tiny inline self-check that the
      JSON strings round-trip through `json.loads` in `tools/zarr_diff.py`

#### Task 0.2: Chunk encode + write (codec `none`, C-order) + round-trip gate

**Files:**
- Modify: `src/io/mod_io_zarr.F90`
- Create: `src/drivers/fesom_zarrsmoke.F90`
- Create: `tools/run_zarrsmoke.sh`
- Create: `tools/zarr_diff.py`

- [ ] `zarr_write_chunk(store, arr, chunk_index(:), data)` — transpose Fortran column-major → **C row-major**,
      write the chunk file at the `dimension_separator`-joined index path, codec `none` (raw little-endian bytes)
- [ ] last-/partial-chunk handling: pad the final global chunk to full chunk size with `fill_value` (the *only*
      padding; standard Zarr — see Technical Details)
- [ ] `fesom_zarrsmoke` driver: write a 1-D `f8`, a 2-D `f4`, and a 2-D `i4` array (known values) to a store
- [ ] `tools/zarr_diff.py`: `--roundtrip` mode opens the store with xarray/zarr and asserts values byte-match the
      generator (`max|Δ|=0` for f8/i4; exact f4)
- [ ] `tools/run_zarrsmoke.sh`: build + run the driver + run the Python round-trip; **gate green** = exit 0

#### Task 0.3: lz4 codec + optional `.zmetadata` consolidation

**Files:**
- Modify: `src/io/mod_io_zarr.F90`
- Modify: `CMakeLists.txt` (discover + link `liblz4` as an INTERFACE target, mirroring `fesom_netcdf`)
- Modify: `tools/zarr_diff.py`, `tools/run_zarrsmoke.sh`

- [ ] confirm `liblz4` (+ headers) linkable on Levante; add `find_library(LZ4 ...)` + `fesom_lz4` INTERFACE
      target linked into `fesom3` (guarded: if absent, compile `none`-only and `error stop` on lz4 request)
- [ ] lz4 codec: numcodecs framing (4-byte little-endian decompressed length + lz4 block), `.zarray`
      `compressor:{"id":"lz4"}`
- [ ] `zarr_consolidate(store)` — optional `.zmetadata` (concatenate all `.zarray`/`.zattrs`) for fast opens
- [ ] **gate:** extend the round-trip to also write with `compressor=lz4` and `--consolidated`; xarray reads both
      and values match; `none` path still green

### Stage 1 — `fesom.mesh.diag.zarr` (proves the whole stack on static data)

#### Task 1.1: `mod_io_decomp` — canonical chunk map + writer-subset + 1-rank identity

**Files:**
- Create: `src/io/mod_io_decomp.F90`

- [ ] `decomp_init(C, n_writers, partit, mesh)` — compute, for node and elem entities: canonical global size
      `N`, chunk count `ceil(N/C)`, block assignment of chunks → `n_writers` writer ranks, and this rank's
      send-plan (for each owned entity, its canonical id `myList(i)` → chunk → destination writer rank). Lazy
      cache. Optional-`partit` (absent ⇒ 1-rank identity: one writer holds all in canonical order).
- [ ] `decomp_redistribute_2d(field_local, buf_writer, entity)` — one `MPI_Alltoallv` from compute layout →
      canonical-chunk layout on the writer subset (zero recvcounts for non-writers)
- [ ] `decomp_redistribute_3d(field_local, buf_writer, entity)` — loop levels (or per `vert_chunk`), reuse the 2-D
      redistribute per level → bounded in-flight memory `O(C·nlev)` per writer
- [ ] helper `decomp_is_writer(rank)` + the writer's chunk-index range + canonical placement offsets
- [ ] **gate:** 1-rank — a labelled array (`value(i)=myList(i)`) round-trips through `decomp_*` and lands in
      canonical order (`buf[g]==g`); deferred multi-rank assertion to Task 1.3

#### Task 1.2a: `mod_io_meshdiag` — coords + connectivity + UGRID topology (1-rank) + FESOM2 gate

**Files:**
- Create: `src/io/mod_io_meshdiag.F90`
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90` (call `meshdiag_write()` after setup, env/namelist-gated)
- Create: `tools/run_meshdiag_gate_core2.sh`
- Modify: `tools/zarr_diff.py` (mesh.diag-vs-FESOM2 mode)

- [ ] `meshdiag_write(path, partit, mesh)` core: `lon`/`lat` (`geo_coord_nod2D` rad→deg), `face_nodes`
      (`elem2D_nodes(1:3,:)`, 1-based, `start_index=1`), `edge_nodes`, `face_edges`, `face_links` (−999),
      `edge_face_links`, `nz`/`nz1` (sign-flip)
- [ ] the `fesom_mesh` UGRID topology variable (cf_role/topology_dimension/node_coordinates/
      face_node_connectivity/…) + `Conventions="UGRID-1.0"` global attr + `_ARRAY_DIMENSIONS` on every var
- [ ] wire `meshdiag_write()` into the driver after setup (gated by `FESOM3_MESHDIAG`/namelist)
- [ ] `tools/zarr_diff.py --meshdiag`: open `fesom.mesh.diag.zarr` + FESOM2 `fesom.mesh.diag.nc`; assert
      connectivity/integer fields **exact**, coords `max|Δ|` ≤ round-off — over the **emitted subset only**
      (cavity/partial-cell vars `ulevels*`/`zbar_*_surface` are OFF in this config and intentionally omitted;
      `gradient_vec_x/y` deferred — see 1.2b)
- [ ] **gate (1-rank, CORE2):** `tools/run_meshdiag_gate_core2.sh` green + `ushow <store>` opens (smoke)

#### Task 1.2b: mesh.diag derived / diagnostic fields

**Files:**
- Modify: `src/io/mod_io_meshdiag.F90`, `tools/zarr_diff.py`

- [ ] add `elem_area`, `nlevels`/`nlevels_nod2D`, `nod_in_elem2D`/`_num`, `edge_cross_dxdy`,
      `gradient_sca_x/y` (`gradient_sca(1:3/4:6,:)`), `nod_area`
- [ ] computed/derived: `zbar_e_bottom` (sign-flip), `zbar_n_bottom` (= `zbar(nlevels_nod2D(n))`, sign-flip),
      `nod_part` (`partit%part`), `elem_part` (stamp `mype` over `myList_elem2D`)
- [ ] **DEFER** `gradient_vec_x/y` — not computed in FESOM3 (opt_visc=7 never builds it); transcribe
      `compute_gradient_vec` from FESOM2 `oce_mesh.F90` only if full pyfesom2 parity is later required (➕)
- [ ] **gate:** extend `--meshdiag` to the full emitted set; `max|Δ|` ≤ round-off vs FESOM2 (ints exact)

#### Task 1.3: mesh.diag multi-rank (partition-independence gate) ✅ proves the parallel stack

**Files:**
- Modify: `src/io/mod_io_meshdiag.F90` (route arrays through `mod_io_decomp`)
- Modify: `src/io/mod_io_decomp.F90` (fixes surfaced at MR)
- Create: `tools/run_meshdiag_gate_multirank.sh`

- [ ] **store-create ordering** (reuse for Task 2.3): rank 0 `zarr_create_store` + all `zarr_define_array`
      (dirs + `.zarray`/`.zattrs`) → `MPI_Barrier` → writers write chunk data → rank 0 (optional) consolidate —
      avoids writers racing a not-yet-created `<var>/` dir
- [ ] route each mesh.diag array through `decomp_redistribute_*` (canonical, writer-subset); each writer writes
      only its chunks
- [ ] verify `n_writers` subset path (e.g. `n_writers=2` with 8 ranks) writes the same store as all-writers
- [ ] **gate:** `dist_2 ≡ dist_8` store value-identical (`max|Δ|=0`, partition-independence) **and** still
      `== FESOM2 mesh.diag`; `ushow` opens the MR-written store
- [ ] **no-regression:** `ctest` 13/13 + one production MR lifecycle byte-gate still `max|Δ|=0`

### Stage 2 — field output (`mod_io_means`)

#### Task 2.1: registry + snapshot + frequency + per-variable-per-year stores (1-rank, node scalars)

**Files:**
- Create: `src/io/mod_io_means.F90`
- Create: `config/namelist.io` (FESOM3 output config)
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90` (`output_init` + `output(istep)` in the step loop)
- Create: `tools/run_output_gate_core2.sh`
- Modify: `tools/zarr_diff.py` (field-vs-FESOM2 + round-trip modes)

- [ ] `t_output_field` (name, long_name, units, location node|elem, ndims, source pointer, freq, mean|snap,
      precision, chunk_shape) + `t_output_stream` registry; `output_register_field(...)`
- [ ] register node scalars **T, S, ssh, sst, sss, a_ice, m_ice, m_snow** (snapshot first) with source-array
      pointers (`tracers%data`, `dyn%eta_n`, `ice%data`)
- [ ] frequency check (`y/m/d/h/step` + interval) from `mod_clock`; CF `time` coord (`"seconds since <start>"`
      **+ `time:calendar`** from `forc_calendar` — `gregorian` (JRA55) / `noleap` (CORE2); L49-sensitive);
      per-variable-per-year store create + append. **Pin `chunk_time=1`** in v1: one record ⇒ one fresh chunk
      file per writer per step, so append = write new chunks + bump `.zarray` `shape[0]` + extend the `time`
      coord (NO read-modify-write of a partial time-chunk; `chunk_time>1` deferred to Task 2.6)
- [ ] embed `lon`/`lat` (node) coords in each data store; parse `namelist.io` (staged into the **rundir** like
      the oracle namelists; `FESOM3_*` env overrides — no repo `config/` runtime dependency)
- [ ] wire `output_init` + `output(istep)` into the driver — **immediately after `step_oce` (~line 730), ABOVE
      the `if (.not. step_diag) cycle` at line 740** (else output is skipped on every production step); gated by
      `FESOM3_OUTPUT`/namelist
- [ ] **gate (1-rank, snapshot):** store opens in xarray (coords+time present); values `==` in-memory state
      (round-trip, `max|Δ|=0`); **vs FESOM2 snapshot output** chase `max|Δ|=0` (float32)

#### Task 2.2: mean accumulation + interval reset (1-rank), FESOM2-aligned

**Files:**
- Modify: `src/io/mod_io_means.F90`

- [ ] running **sum ÷ count** accumulator per field; accumulate each step; divide + reset each output interval
      (transcribe `io_meandata.F90:update_means` semantics: which steps, divide timing)
- [ ] **gate:** monthly-mean field vs FESOM2 monthly mean — chase `max|Δ|=0` (float32); round-trip still green

#### Task 2.3: multi-rank node scalars (partition-independence)

**Files:**
- Modify: `src/io/mod_io_means.F90` (route writes through `mod_io_decomp`)
- Create: `tools/run_output_gate_multirank.sh`

- [ ] route the accumulated field through `decomp_redistribute_2d` (writer subset) before the chunk write
- [ ] **gate:** `dist_2 ≡ dist_8` (`max|Δ|=0`, partition-independence); still `== FESOM2`; round-trip
- [ ] **no-regression:** `ctest` 13/13 + a production MR lifecycle byte-gate `max|Δ|=0`

#### Task 2.4: 3D fields (T, S full-depth; w on nodes) + vertical chunking

**Files:**
- Modify: `src/io/mod_io_means.F90`, `src/io/mod_io_decomp.F90`

- [ ] register 3D **T, S** (`nl-1` layers) and node **w** (`nl`); use `decomp_redistribute_3d` (level-by-level);
      honor `vert_chunk`; `_FillValue` for below-bottom levels (`> nlevels_nod2D`) — standard CF masking, not
      chunk padding
- [ ] **gate:** 3D field vs FESOM2 (chase `max|Δ|=0`); partition-independence; round-trip; `ushow` opens a level
      slice

#### Task 2.5: element vectors u, v + r2g rotation (geographic default, native knob)

**Files:**
- Modify: `src/mesh/mod_mesh_rotate.F90` — **port `vector_r2g`** (rotated→geographic vector transform)
- Modify: `src/io/mod_io_means.F90` (vector-pair registry + rotation), `src/io/mod_io_meshdiag.F90` (elem-centroid
  `lon`/`lat` coords if not already emitted)
- Modify: `config/namelist.io` (vector-frame knob)

- [ ] **port `vector_r2g`** from FESOM2 `gen_modules_rotate_grid.F90` into `mod_mesh_rotate.F90` (cite file:line):
      the existing module has only `vector_g2r` (hard-wired geo→rotated, **not** flag-invertible); the scalar
      `r2g` already there is reused for the elem-centroid geographic coords
- [ ] register element **u, v** (`dyn%uv(1:2,:,:)`) as a **vector pair**; emit/embed elem-centroid `lon`/`lat`
- [ ] apply `vector_r2g` to the (u,v) pair at write time using the elem-centroid coords
- [ ] `vec_frame` namelist knob: `geographic` (default) | `native`
- [ ] **gate:** u, v vs FESOM2 (matching `vec_autorotate`) chase `max|Δ|=0`; `native` variant matches the
      unrotated element values; partition-independence holds

#### Task 2.6: full namelist.io knobs + float32/lz4/chunk-shape/n_writers/filesplit

**Files:**
- Modify: `src/io/mod_io_means.F90`, `config/namelist.io`

- [ ] expose all knobs: per-field `freq`/`unit`/`precision`/`mean|snap`/`frame`; global `n_writers`,
      `chunk_shape (time,vert,horiz)`, `compressor (none|lz4)`, `filesplit_freq`; float32 default
- [ ] `chunk_time>1`: implement the partial-last-time-chunk **read-modify-write** append path (deferred from
      Task 2.1's `chunk_time=1`); gate that `chunk_time={1,N}` produce value-identical stores
- [ ] **gate:** a run with `compressor=lz4` + custom `chunk_shape` + `n_writers` subset still passes
      partition-independence + round-trip + `ushow` opens it

### Final

#### Task F1: Full no-regression + gate sweep

- [ ] `cd build_intel_dp && ctest --output-on-failure` → 13/13 (repeat GNU)
- [ ] existing production MR byte-gates (`run_lifecycle_*_gate_multirank.sh`, both whichEVP) stay `max|Δ|=0`
- [ ] all new M9 gates green (Stage 0 round-trip; mesh.diag 1-rank + MR; fields 1-rank + MR + 3D + vectors)
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
