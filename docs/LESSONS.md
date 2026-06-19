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
**VERIFIED at the M1 geometry gate:** both build with ifort 2021.5.0 via the same mpif90,
and the FESOM2 build flags (`build/src/CMakeFiles/fesom.dir/flags.make`) are byte-for-byte
the anchor list above. Same compiler + same flags ⇒ scalar SVML transcendentals match.

## L7 — `-no-prec-div`: divide by a LITERAL, not a runtime `real(nv)` (arity trap)

The single hardest bug in the M1 geometry gate. The anchor flag `-no-prec-div` lets ifort
implement `x/y` as `x * recip(y)`. For a **compile-time constant** divisor (`x/3.0_WP`)
ifort folds an exact `1/3` and the result is one value; for a **runtime** divisor
(`x/real(nv,WP)`, nv from `elem2D_nnodes`) it emits a hardware reciprocal approximation
that differs by **1 ULP**. That 1 ULP in a centroid latitude `sum(lat)/3` flowed into
`cos`/`tan` (steep near the rotated pole → up to ~25% of elements differ), then
`elem_area`, `gradient_sca`, `area`, `edge_cross_dxdy` — the whole geometry failed the
gate. **Fix:** divide by the literal `3.0_WP`/`3.0_MP` (what FESOM2 writes), NOT
`real(nv,WP)`. This is the plan's arity caveat made concrete: the `/3 → /n_vert`
generalization is NOT universally bit-safe. v1 is triangles only (nv==3); quad support
must branch on nv and divide by the matching literal. **Non-causes ruled out** (cost a
day if not): the intrinsic `sum()` vs an explicit accumulation loop are bit-identical here;
auto-vectorization of `cos` was NOT the cause (`!DIR$ NOVECTOR` had zero effect — removed);
coords + `elem2D_nodes` were proven `max|Δ|=0` first, which is what localised it to the
divide. Lesson: when a geometry byte-gate shows ~1-ULP transcendental diffs, suspect a
runtime divisor under `-no-prec-div` before anything else.

## L8 — 1-rank FESOM2 oracle (hand-crafted dist_1) + geometry byte-faithfulness

To byte-gate anything that accumulates per node/element (`area`, later FCT scatter), the
oracle must run **1-rank**, because FESOM2 builds `nod_in_elem2D` in LOCAL element order
(oce_mesh.F90 find_neighbors), so a multi-rank run reorders the per-node `Σ elem_area/3`
→ ULP diffs vs FESOM3's global order. METIS can't produce `dist_1`; hand-craft it
(see HANDOFF "Geometry byte-gate"). FESOM2's reader accepts blank lines for zero-size
halo arrays (zero-trip `read(*,*)` skips a record). 1-rank forcing init hangs on a login
node — irrelevant: dump geometry at `mesh_setup` and STOP before forcing.
Geometry transcription gotchas the gate caught (all now `max|Δ|=0`):
- `elem_center`/`edge_center` wrap arithmetic must be VERBATIM (FESOM2 wraps each lon vs
  `amin`; edge_center shifts a(1)/b(1) asymmetrically). An "equivalent" rewrite diverged.
- **Vertex (CW) order — gated only on the 0-swap case so far.** pi has 0/5839 elements
  needing the clockwise swap (it's pre-oriented), so the geometry gate did NOT exercise
  `enforce_cw_orientation`'s reorder path. soufflet has 228/5700, CORE2 ~certainly >0.
  Within-element node order is FP-order-sensitive (centroid `sum(lat)/3`; `gradient_sca`
  column→vertex map; elem_area is swap-invariant). FESOM3 `enforce_cw_orientation` is
  byte-identical to FESOM2 runtime `test_tri` (oce_mesh.F90:1706-1726: same b/c, same
  `trim_cyclic` on comp 1, same `r=b1*c2-b2*c1`, same `r>0`→swap 2,3) operating on the
  same elem2d.out + byte-proven coords ⇒ identical swaps by construction. **But empirically
  confirm on a swap mesh (soufflet, or the M2.11 single-rank CORE2 gate) before trusting it.**
- `elem_area`/`area`/`areasvol` accumulate UNSCALED then `*r_earth²` ONCE at the end
  (a single deferred scaling), matching mesh_areas — not per-element scaling.
- `edge_tri.out` stores **-999** for boundary (no 2nd triangle); FESOM2 `load_edges`
  does `where(edge_tri<0) edge_tri=0`. Replicate, or the boundary marker leaks.

## L6 — Oracle is runnable; how to drive it (proven 2026-06-19)

Running the prebuilt instrumented FESOM2 on pi to produce reference dumps works and is
fast (~0.25 s, 2 steps, 2 ranks). Gotchas, all handled by `tools/run_oracle_pi.sh`:
- `fesom.clock` must be non-empty even for a fresh start: two lines `0 1 1948` (the
  forcing year; columns = seconds, day-of-year, year). Empty → `severe (24): end-of-file`.
- `run_length_unit='s'` means STEPS (not seconds) — set `run_length=2` for a 2-step run.
- `ClimateDataPath` in the shipped namelist is relative (`..//test/input/global/`); from a
  scratch run dir it must be absolute. Forcing paths were already absolute.
- Login-node MPI needs `--mca pml ob1 --mca btl self,vader` (UCX/IB settings fail).
- A fresh run byte-matches a prior run's dump (`max|Δ|=0`) → FESOM2 is deterministic, so
  `max|Δ|=0` byte-gating is well-posed. The shipped pi config runs KPP/GM (substep MIXING
  dumps `Kv`); the **reduced M2 config** (PP, no GM/Redi, linfs, opt_visc=7) is a separate
  namelist to apply when byte-gating M2.

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
