# FESOM3 — Lessons (mistakes + fixes, never repeated)

Append-only log of non-obvious gotchas. Each entry: what bit us, why, the fix.

## L1 — Anchor compile flags come from FESOM2 v2.7.3, NOT the tracer_dwarf

The byte-gate oracle is **FESOM2 v2.7.3** (`/home/a/a270088/port2/fesom2/src/CMakeLists.txt`).
The tracer_dwarf is a *structural* reference whose flags differ and were never byte-gated.
Differences that would silently break bit-identity if we copied the dwarf:

- dwarf adds `-no-prec-sqrt -ip` to **ifx** too; FESOM2 restricts them to **ifort classic**.
- dwarf adds `-march=core-avx2 -mtune=core-avx2` on Levante Intel; **FESOM2 v2.7.3 COMMENTS
  THAT LINE OUT** (`src/CMakeLists.txt:356`). So the FESOM2 Levante Intel build targets the
  SSE2 baseline (no AVX2, **no FMA contraction**). Adding `-march=core-avx2` enables FMA →
  different bits. Our anchor therefore uses **no `-march`** on Levante.

Faithful anchor (Intel classic, Levante), from FESOM2 v2.7.3 `src/CMakeLists.txt:335-340`:
```
-O3 -r8 -i4 -fp-model precise -no-prec-div -fimf-use-svml -init=zero
-no-wrap-margin -fpe0 -fpp  -no-prec-sqrt -ip
```
GNU (FESOM2 `:396`,`:434`): `-O3 -ffloat-store -finit-local-zero -finline-functions
-fimplicit-none -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -cpp
-march=znver3 -mtune=znver3 -ftree-vectorize -flto`.

Codified in `cmake/fesom_flags.cmake`. `-fpe0` implies flush-to-zero (FTZ) of denormals —
itself a bit-affecting behavior we must match, so keep it. `-init=zero` zeroes locals
(prevents spurious `-fpe0` traps on uninitialised reads). When a real byte-gate runs
(M0.7+), re-verify the actual FESOM2 build's flags from its build dir, not just the CMake.

## L5 — Mesh geometry (M0.7): ocean-only meshes, units, deferred byte-gate

- **pi is OCEAN-only** (no land elements). Σ elem_area = ocean area ≈ 3.40e14 m² ≈
  0.67·(4πr²), NOT the full sphere. Don't assert full-sphere area; assert the ocean
  fraction (~0.5–0.85). pi extent is global (lon 0–360, lat −78..89), nl=48.
- **FESOM2's geometry convention is radians × r_earth**: coords stored in radians,
  `elem_area *= r_earth²`. A Cartesian analytic mesh must therefore store coords as
  *radians-like* (physical metres / r_earth) so the shared pipeline yields physical
  m². Test a linear field with the physical position = coord·r_earth.
- **gradient_sca assumes CW node order** — run CW enforcement (test_tri) BEFORE
  computing geometry, on file AND analytic meshes. (pi files are pre-oriented: 0 swaps.)
- **FESOM2 geometry byte-gate is DEFERRED to M1**: it needs an instrumented-FESOM2
  run dumping reference elem_area/gradient_sca/areas, which wasn't producible this
  session. M0.7 is gated by self-consistency (CW, ocean-area, gradient annihilates
  constants, gradient reproduces a linear field, adjacency). Deferred geometry not yet
  transcribed: gradient_vec (M2 momentum), mesh_resolution smoothing (M4 GM), and the
  multi-rank mesh remap (M2.12). Edge left/right orientation on the analytic mesh is
  not FESOM2's (analytic only; pi reads edge_tri from file).

## L4 — Type-design choices (M0.3) deviating from the plan letter

Documented here per the plan's "update scope when implementation deviates":
- **Field names kept FESOM2-verbatim** (`elem2D_nodes`, `nod_in_elem2D`, `nlevels`,
  `gradient_sca`, ...) so kernels transcribe line-for-line. The plan's `elem_nodes`/
  `elem_nnodes` naming → realized as `elem2D_nodes(MAX_NV,:)` + new `elem2D_nnodes(:)`.
- **node→element adjacency is DENSE** (`nod_in_elem2D(MAX_ADJACENT,:)` + `_num`), not
  CSR as D2 suggested. FESOM2 and the dwarf both use dense; matching it preserves
  kernel iteration order for byte-identity. (CSR was a memory aspiration, not in the refs.)
- **Types are lean/growable**: dropped beyond-v1 fields (cavity, iceberg, OASIS, split-
  explicit `se_*`, backscatter/UKE, energy diagnostics `ke_*`, DVD). Aux 3D WORK fields
  (density/N²/Kv/Av/hpressure/sw_α-β) are added to `t_dyn_work` at M2.1 with verified
  FESOM2 names rather than guessed now. GM `fer_uv/fer_w` declared, serialized at M4.
- **Serialization = derived-type `read/write(unformatted)`** (DTIO) over an `access='stream'`
  unit; size-prefixed via `mod_binary_arrays`. Only WP-real variants exist — at DP/SP
  MP==WP so MP mesh arrays serialize through them; HP compiles IO out. Full restart-file
  orchestration (`mod_io_restart`) deferred to M2.11.

FESOM2's GNU Levante flags include `-flto`, but FESOM2 links one big executable;
FESOM3 links a static `libfesom3.a` + separate driver/test executables. `-flto` puts
LTO objects in the archive, and gfortran 8.5's `ar`/`ranlib` need the LTO plugin
(`gcc-ar`/`gcc-ranlib`) plus `-flto` at link — the default tools fail with
"plugin needed to handle lto object" and unresolved module symbols. GNU is the
portability build, not the bit-identity anchor (that's Intel), and `-flto` can reduce
FP reproducibility, so we omit it. Revisit only if GNU↔FESOM2-GNU byte-identity is wanted.

## L2 — Precision: `-r8` promotion AND explicit `_WP` kinds move together

FESOM2 relies on `-r8`/`-fdefault-real-8`: unsuffixed literals like `1.0` are **double**.
Transcribe literals verbatim (keep unsuffixed where FESOM2 leaves them unsuffixed) AND build
with `-r8`, or the bits differ. The cpp define `USE_SINGLE_PRECISION` (sets `WP=4`) and the
`-r4` flag are flipped together by `configure.sh --precision sp`. `MP = max(WP,4)`.
