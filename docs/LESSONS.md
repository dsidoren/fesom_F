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

## L9 — Prior gates validate hidden orderings; the kernel-gate recipe (M1.1, PASSED first try)

M1.1 (horizontal tracer advection: upw1 + MUSCL + scatter) hit `max|Δ|=0` vs FESOM2 on
the FIRST gate run — no bug-hunt. Two reasons, both reusable:

- **A passing gate transitively proves orderings the next gate depends on.** M1.1's
  `find_up_downwind_triangles` picks a triangle by `atan2` comparisons, and `fill_up_dn_grad`
  area-weight-AVERAGES `tr_xy` over `nod_in_elem2D(:,node)` — both sensitive to the ORDER of
  `nod_in_elem2D`. That order was never dumped directly, but the **geometry gate's `area`
  match (`area(nz,n)=Σ_k elem_area(nod_in_elem2D(k,n))/3`, an FP sum) already pinned it**:
  if the k-order differed, `area` would have differed. So `edge_up_dn_tri` (integer indices)
  and `edge_up_dn_grad` matched by construction. Lesson: when designing a gate, list the
  order-sensitive inputs and check whether an earlier `max|Δ|=0` already constrains them —
  it often does, and that's why faithful transcription "just works."

- **The kernel-gate oracle recipe** (now proven, reuse for M1.2+): wire an env-gated,
  npes==1 dump shim into FESOM2 at the LATEST setup point that has what you need but is
  still before the 1-rank forcing hang (M1.1 = end of `ocean_setup`, after
  `init_thickness_ale`). The shim PRESCRIBES analytic inputs from the byte-proven coords
  (so FESOM2 and FESOM3 generate identical `ttf`/`vel` independently — no captured-input
  file) and calls the **real FESOM2 kernels**, then STOPS. Gate every intermediate, not
  just the final target, so a divergence localises itself. Pin scalar knobs (dt, num_ord)
  as shared constants in BOTH codes, not from a namelist.

- **`helem` at init is `which_ale`-independent.** Both linfs and zstar `init_thickness_ale`
  reduce to `helem(nz,e)=zbar(nz)-zbar(nz+1)` (full cells) at the initial state because
  `hbar=eta=0` there. So a kernel gate that dumps at init can compute `helem` from `zbar`
  directly without importing ALE thickness evolution (deferred to M2.7). Verify it as a
  dumped field anyway (M1.1 did — `max|Δ|=0`), since partial cells/cavity would break it.

## L10 — Vertical advection (M1.2): QR4C wants ALE 3D depths, not static Z/zbar (PASSED first try)

M1.2 (vertical advection: upw1 + QR4C + scatter) hit `max|Δ|=0` vs FESOM2 on the FIRST gate
run, like M1.1. Reusable specifics:

- **QR4C divides by the *per-node ALE* mid-depths `Z_3d_n`/`zbar_3d_n`, NOT the static 1-D
  `Z`/`zbar`.** FESOM2 builds them in `init_ale` (oce_ale.F90:531-566). At the initial state
  (no cavity, full cells) they reduce to: `zbar_3d_n(nz,n)=zbar(nz)` for all nz; interior
  `Z_3d_n(nz,n)=Z(nz)`; and the surface/bottom interfaces use
  `Z_3d_n(nzmin)=zbar(nzmin)+(zbar(nzmin+1)-zbar_n_srf)/2` and
  `Z_3d_n(nzmax-1)=zbar(nzmax-1)+(zbar_n_bot-zbar(nzmax-1))/2` with `zbar_n_srf=zbar(nzmin)`,
  `zbar_n_bot=zbar(nzmax)`. The FESOM2 shim used the LIVE `mesh%Z_3d_n` (init_ale ran earlier
  in ocean_setup); the FESOM3 driver REPLICATES that formula. Both `zbar_3d_n` and `Z_3d_n`
  were dumped as gated fields (`max|Δ|=0`) so the replication is proven, not assumed. NOTE the
  geom gate did NOT cover `Z`/`zbar` (it dumps level *counts* only); the advhor gate now does
  (via the 3-D arrays). FESOM3 `read_mesh` builds `Z=0.5*(zbar(k)+zbar(k+1))` in ONE step vs
  FESOM2's two steps (`Z=zbar(k)+zbar(k+1); Z=0.5*Z`) — byte-identical because ×0.5 is exact.

- **`/2` is byte-exact under `-no-prec-div`; `/3` is not (cf. L7).** The boundary `Z_3d_n`
  `/2` and the QR4C `/3.0_WP` literal divisors fold to exact/identical constants. And QR4C's
  `qc/qu/qd = (ttf diff)/(Z_3d_n(k)-Z_3d_n(k+1))` are RUNTIME divisors — fine here because
  BOTH codes divide by the *same* runtime operand (byte-identical `Z_3d_n`), so the reciprocal
  approximation matches. L7's trap was a runtime-vs-LITERAL *mismatch* between the two codes,
  not runtime division per se. Transcribe the divisor expression verbatim and operands match.

- **The "D=2 shallow-column double-write" is moot on pi** (its min column is `nlevels_nod2D=5`
  = 4 layers; `sort -n nlvls.out | head -1`). The QR4C 2-layer double-write (2nd-layer and
  bottom-1 both hitting interface `nzmin+1`, which nets the flux there to ~0) only triggers at
  `nlevels_nod2D=3`; a 1-layer column (`=2`) would read `ttf(0)` OOB. pi has neither, so
  faithful statement-order transcription is automatically correct. Confirmed with a Debug
  (`-check all`) run of `fesom_advhordump` (EXIT 0, no OOB/FPE). Re-examine for meshes with
  shallower minima (soufflet, CORE2) — there the double-write/OOB order matters and must match.

- **Debug build is NOT byte-comparable to the Release oracle.** A `-check all` (`-O0`, no
  `-no-prec-div`) FESOM3 dump differs from the `-O3` FESOM2 oracle on EVERY FP-computed field
  (horizontal too, and even `area` at ~3.7e-4) — a uniform flag effect, not a bug. Use Debug
  ONLY to catch OOB/uninit/FPE; byte-gate Release-vs-Release (anchor flags both sides, L1).

## L11 — FCT (M1.3): config knobs, the ttf-gradient quirk, the bignumber bottom-fill (PASSED first try)

M1.3 (FCT/Zalesak limiter `oce_tra_adv_fct` + the `use_lo` scatter) hit `max|Δ|=0` vs FESOM2 on
the FIRST gate run, like M1.1/M1.2. The 15 new fields (fct_LO, clipped adf_h/adf_v, fct_ttf_max/
min, fct_plus/minus, del_ttf_*_fct, hnode/hnode_new, standalone MFCT) all matched. Specifics:

- **The pi FCT config is MFCT/QR4C/FCT with `opth=0.0`, `optv=1.0`** (namelist.tra tracer 1/2:
  `1,'MFCT','QR4C','FCT ',0.,1.` → hor.Ord=0 ⇒ 3rd-order MFCT, vert.Ord=1 ⇒ 4th-order QR4C), NOT
  the `num_ord=0.75` M1.1/M1.2 used for the standalone MUSCL/QR4C dumps. Pin opth/optv as shared
  constants in BOTH codes. The blended num_ord term is already proven (M1.1 MUSCL@0.75); a
  standalone `adv_flux_hor_mfct`@0.75 dump finally gates the MFCT kernel (no-clamp variant) fully.

- **`edge_up_dn_grad` is the gradient of `values` (ttf), NOT `valuesAB` (ttfAB).** FESOM2
  `init_tracers_AB` (oce_tracer_mod.F90:126-127) has the `valuesAB` gradient call COMMENTED OUT
  and uses `values`; then `do_oce_adv_tra` reconstructs `valuesAB` (ttfAB) in MUSCL/MFCT/QR4C
  using that ttf-gradient (note line 131 "WHY NOT AB HERE? DSIDOREN!"). So the gate reuses
  `eudg=grad(ttf)` (already built for M1.1) for the `MFCT(ttfAB)` call — do NOT recompute
  grad(ttfAB). LO fluxes use ttf, HO fluxes use ttfAB, both with eudg=grad(ttf).

- **The `a2` bignumber fill clobbers EVERY element's bottom layer.** `do nz=nu1,nl1-1` (real
  bounds) then `if(nl1<=nl-1) do nz=nl1,nl-1` with `nl1=nlevels(elem)-1`, so layer
  `nlevels(elem)-1` (the element's own bottom layer) is set to ∓bignumber, leaving the deepest
  node layer's admissible increment unconstrained. An off-by-one (`nl1` vs `nl1-1`) here would
  diverge — pi has shallow elements so the fill is exercised and gated.

- **AUX scratch: FESOM2 reuses `edge_up_dn_grad`; FESOM3 uses a LOCAL.** Byte-identical on pi
  because `a2` writes `AUX(1:2, ulevels(elem):nl-1, elem)` for every element and (with
  ulevels==1, no cavity) `a3` only reads written entries — the uninitialised local never reaches
  a result. On a CAVITY mesh `a3` reads `AUX(:,nz,elem)` at `nz<ulevels(elem)` (FESOM2 = stale
  edge_up_dn_grad there); unreproducible with a fresh array → **re-gate FCT on a cavity mesh
  (M2+).** In the SHIM, pass a SEPARATE `aux_scratch` (not edge_up_dn_grad) to the real
  `oce_tra_adv_fct` so edge_up_dn_grad survives intact for its own dump record. `dmax1/dmin1`→
  generic `max/min` (identical IEEE at WP=8). exchange_nod(fct_plus,fct_minus) is a 1-rank no-op.

- **`hnode/hnode_new` at init = `zbar(nz)-zbar(nz+1)`** (node analog of helem, full cells;
  `hnode_new=hnode`, oce_ale.F90:1217), built like helem and gated. Used by the LO vertical
  update, the b2 limiter divisor, and the `use_lo` scatter reconstruction. Runtime divisors
  `areasvol`/`hnode_new` match FESOM2 (same byte-identical operands; cf. L7/L10).

- **A strong FCT gate needs ttfAB sharp+large vs the smooth ttf, or the limiter never clips.**
  The b2 increment is `flux·dt/(areasvol·hnode_new)` with `dt/(areasvol·hnode_new) ~ 1.8e-10` on
  pi, so for smooth realistic-magnitude fields the antidiffusive increment sits far below the
  admissible bound ⇒ `ae≡1` ⇒ b1/b2/b3 are computed but multiplied by 1, masking transcription
  bugs in the limiting *selection*. A high-wavenumber, amplitude-20/10 ttfAB (vs amplitude-1 ttf)
  drove clipping at ~41% of nodes (fct_plus/minus min=0), changing ~14% of the antidiffusive
  fluxes — so the sign-based clip selection was genuinely exercised. Verify clipping fraction
  post-run (count fct_plus<1); byte-identity holds regardless, but gate STRENGTH needs it.

## L12 — Assembled step (M1.4): gate the REAL driver, not an inline copy (PASSED first try)

M1.4 (the assembled advection step: `do_oce_adv_tra` + `init_tracers_AB` + the `adv_tracers_ale`
`del_ttf += advhoriz+advvert` accumulation) hit `max|Δ|=0` vs FESOM2 on the FIRST gate run.
Unlike M1.1–M1.3 (the shim INLINED the orchestration), M1.4 drives FESOM2's OWN `init_tracers_AB`
+ `do_oce_adv_tra`. The driver-gate recipe (reuse for every future assembled-routine gate):

- **Run the real FESOM2 routine in the oracle shim AFTER the inline-kernel records are dumped.**
  The driver overwrites the shared `tracers%work` (adv_flux_*, fct_LO, fct_*, edge_up_dn_grad as
  FCT scratch) and `tracers%data(1)`. Putting the driver section LAST (just before `close(u)`)
  needs ZERO edits to the proven M1.1–M1.3 dump — purely additive, regression-safe (re-run the
  whole gate: the 34 prior fields must still PASS, which they did). The alternative (driver-first)
  works only because the inline section rebuilds work from scratch; last-is-simpler.

- **Prescribe the driver's INPUTS, not its intermediates, or the new numeric never runs.** M1.3
  prescribed `valuesAB`(ttfAB) directly; that never exercises FESOM2's AB interpolation. M1.4
  prescribes `values`(smooth ttf) + `valuesold(1)`(sharp ttfAB) so `init_tracers_AB` COMPUTES
  `valuesAB = -(0.5+ε)·valuesold + (1.5+ε)·values` (ε=0.1) and the gate compares that.

- **The AB-offset `+` did NOT trip the L7 literal-vs-runtime trap — but gate it anyway.** Feared:
  FESOM2's `epsilon` is a runtime module var (`(1.5_WP+epsilon)` computed live), so if FESOM3
  folds `1.5+0.1→1.6` at compile time the bits could differ under fast-math (as `/` does under
  `-no-prec-div`). They DON'T: `1.5d0+0.1d0` rounds to exactly the literal `1.6d0` (and
  `0.5+0.1→0.6d0`), so `valuesAB` matched byte-for-bit. `+`/`*` reassociation is far less fragile
  than `/` here. Defensive choice kept regardless: `ab_epsilon` is a non-parameter module var in
  mod_config (mirrors FESOM2), so the compiler can't fold even in principle. Lesson: a runtime
  scalar constant in a SUM is usually fold-safe, but the only proof is the oracle gate.

- **pi runs `use_wsplit=.true.`; the gate must FORCE `.false.` on both sides.** With w-split on,
  the FCT path calls `adv_tra_vert_impl` (implicit vertical, NOT ported until M2) + recomputes the
  LO vertical on full `w`. Forcing `dynamics%use_wsplit=.false.` in the shim (= the FESOM3 `t_dyn`
  default) gates the matched EXPLICIT path (`w==w_e`, both passed the single prescribed wvel). Do
  not gate against the production `.true.` setting — it would hit an unported kernel / the FESOM3
  `error stop` guard. The implicit w-split path is a separate M2 gate.

- **Non-FCT (`do_zero_flux`) + per-tracer order knobs closed in the same gate.** A second config
  on tracer 1 (MUSCL/QR4C/NON, ph=pv=0.75) exercises the `do_zero_flux=.true.` dispatch (HO scheme
  applied directly to valuesAB, scatter WITHOUT use_lo) and the tra_adv_ph/pv knobs the M1.3 gate
  fixed at 0/1. FCT vs non-FCT `del_ttf` differ by max 26 yet BOTH byte-match FESOM2 — proving the
  branch selection, not just one path. (The HANDOFF had flagged this gap for M1.4.)

- **WP pointers onto MP work arrays bind by kind VALUE.** `do_oce_adv_tra` associates `real(WP)`
  pointers (`ttf`, `fct_LO`, `adv_flux_*`, ...) with the `real(MP)` `tracers%work`/`tracers%data`
  components and passes them to the WP kernels; legal because MP==WP==8 at the DP/SP anchor (the
  pointer/argument match is on the kind integer 8, not the parameter NAME). The `tracers`/`w`/`we`
  dummies need `target`. Revisit only for FP16 (WP=2, MP=4) — a deferred precision decision.

## L13 — Pressure/EOS/N² (M2.1): split EOS + the horizontal smoother both byte-match (PASSED first try)

M2.1 (`pressure_bv`: split Jackett-McDougall EOS → `density_m_rho0`, top-down `hpressure`, N²
`bvfreq` + the horizontal `smooth_nod` sweep) hit `max|Δ|=0` vs FESOM2 on the FIRST gate run, like
all of M1. The 10 fields (3 inputs T/S/density_ref, 3 depths, density_m_rho0, hpressure, bvfreq raw
+ smoothed) all matched. Reusable specifics:

- **`density_ref == density_0` on pi.** `use_density_ref` defaults `.false.`, so `arrays_init` sets
  `density_ref = density_0 = 1030` and `init_ref_density` is NOT called (`oce_setup_step.F90:218`).
  The plan's "subtract the `density_ref(nz,node)` ARRAY, not the scalar" is therefore bit-TRIVIAL on
  pi (the array is a constant 1030), but transcribe the array form anyway — a cavity/`use_density_ref`
  mesh (M2.11) makes it non-constant. N² divides by the SCALAR `density_0`; the PGF anomaly subtracts
  the ARRAY. Two different density_0 uses in one routine — don't conflate.

- **The horizontal `smooth_nod` (mass-matrix sweep) byte-matched by faithful transcription** — the one
  horizontally-coupled step in M2. Same reason as M1 (L9): it weights by `elem_area` (geom-proven) and
  iterates `nod_in_elem2D(:,n)` / `elem2D_nodes(:,el)` (order area-gate-proven). The `1/(3*Σarea)`
  runtime reciprocal matches because both codes divide the same byte-identical operand (L7/L10). 1-rank:
  the per-cycle `exchange_nod(bvfreq)` is a no-op, dropped (lifted at M2.12). With `N2smth_hidx=1` only
  the first sweep runs (the `do q=1,N_smooth-1` loop is empty). The smoother changed 100% of valid
  entries on the prescribed T/S, so it is a NON-vacuous gate (verify: `max|bvfreq_raw - bvfreq|>0`).

- **Dump bvfreq BOTH pre- and post-smoothing (`bvfreq_raw` + `bvfreq`) to localize.** A two-call toggle
  (`pressure_bv` with `N2smth_h=.false.` then `.true.`) separates an EOS/N²-difference bug (shows in
  bvfreq_raw) from a smoother bug (shows in bvfreq-but-not-raw). Same "gate every intermediate" rule as
  M1. density_m_rho0/hpressure are identical across the two calls (smoothing only touches bvfreq).

- **The caller pre-zeros the three outputs; the kernel does NOT.** FESOM2 `pressure_bv` writes only
  valid levels `nz=nzmin..nzmax-1` (hpressure/density) / `..nzmax` (bvfreq) and leaves the below-bottom
  entries untouched (not zeroed at allocation — heap, not `-init=zero` stack). So the dump's below-bottom
  region is indeterminate unless zeroed. BOTH the FESOM2 shim and the FESOM3 driver `=0` the outputs
  before calling → deterministic 0 there → the rectangular dump compares clean. Faithful: the FESOM3
  kernel is `intent(inout)` and also doesn't zero (matches FESOM2); the zeroing is a gate-harness concern.

- **Force the gate knobs in the shim; KPP/GM vs reduced-M2 does NOT change the EOS fields.** EOS density
  is a pure function of T/S/Z (Jackett-McDougall), so the shipped KPP/GM pi namelist gives the same
  density/hpressure/bvfreq as reduced-M2 — no need to assemble the reduced namelist for M2.1. The shim
  pins only what `pressure_bv` reads: `state_equation=1`, `which_ale='linfs'` (routes hpressure through
  this routine; zstar/zlevel compute it in `pressure_force_4_zxxxx`, M2.x), `N2smth_v=.false.`/
  `N2smth_hidx=1`, `ldiag_dMOC=.false.` (skip density_dmoc), `mix_scheme_nmb=-1` (`mixing_kpp=.false.`
  → dbsfc untouched). MLD1/2/3 + dbsfc + dMOC are KPP-only diagnostics that do NOT feed the gated fields
  (density_m_rho0 is computed BEFORE rho_surf/dbsfc1; the MLD logic only READS bvfreq/rhopot) — omitted.

- **Shared-lib rebuild: `make fesom.x` rebuilds `libfesom.so` but does NOT relink the exe — and that's
  fine.** The shim lives in `libfesom.so`; `fesom.x` resolves it at runtime (`ldd` → build/lib64). After
  editing a shim, `make -C build fesom.x` updates the `.so` (timestamp moves) while the exe stays old;
  running it loads the new `.so`. Confirm with `nm -D libfesom.so | grep <shim_symbol>`. (Same as M1.4.)

## L14 — Hydrostatic PGF (M2.2) byte-matched; and the configure.sh `--debug` Release-clobber footgun

M2.2 (`oce_pgf.F90` `pressure_force_4_linfs_fullcell`: the `gradient_sca` contraction of the M2.1
`hpressure` → element `pgf_x`/`pgf_y`) hit `max|Δ|=0` vs FESOM2 on the FIRST gate run, like all of M1
+ M2.1. The PGF gate is now the 2-field tail (pgf_x/pgf_y) of the same `tools/run_pressure_gate.sh`
(12 fields total). Specifics:

- **PGF is the SAME contraction shape as M1.1's `tracer_gradient_elements`, on a pre-gated input.**
  `pgf_{x,y}(nz,elem) = Σ_k gradient_sca({1:3,4:6},elem)·hpressure(nz,elnodes_k)/density_0`. Both
  operands were already `max|Δ|=0`: `gradient_sca` from the geometry gate, `hpressure` from the M2.1
  pass (computed in the smoothed pressure_bv call — hpressure is identical across the raw/smoothed
  calls, only bvfreq changes, L13). So faithful transcription "just works" (the L9 transitive-gate
  pattern again). `density_0` is a RUNTIME divisor but byte-identical on both sides (=1030), so the
  `-no-prec-div` reciprocal matches (L7/L10). Transcribe VERBATIM: the `/density_0` is INSIDE the sum
  (each of the 3 products divided, then summed), not `sum(...)/density_0`. The shim calls the REAL
  FESOM2 `pressure_force_4_linfs_fullcell` directly (the `pressure_force_4_linfs` dispatcher is pure
  namelist branching — no numerics — so it is not gated, cf. M1.4 which DID gate `do_oce_adv_tra`
  because that dispatch had order-knob numerics). Caller pre-zeros pgf_x/pgf_y (the kernel writes only
  `ule..nle`; FESOM2 leaves below-bottom at its allocation-zero), same as the M2.1 output treatment.
  Gate STRENGTH check: pgf is 66% non-zero at ~1e-5 m/s² (realistic PGF accel) on the prescribed T/S.

- **`configure.sh --debug` overwrote the Release anchor build dir → a spurious `area=3.7e-4` gate FAIL.**
  The footgun: `configure.sh` named the build dir `build_<compiler>_<precision>` with NO build-type
  suffix, so `--debug` reconfigured **build_intel_dp itself to CMAKE_BUILD_TYPE=Debug** and rebuilt its
  binaries at `-O0`/no-`-no-prec-div`. Running a parallel `--debug` build alongside the M1 advhor gate
  meant the gate executed the freshly-clobbered DEBUG `fesom_advhordump` → every FP field diverged,
  with `area` at the tell-tale `3.66e-4` (the exact L10 Debug-vs-Release signature). It LOOKED like an
  M2.2 regression in M1; it was a build-dir collision. **Fix:** `configure.sh` now appends `_debug` to
  the build dir for Debug builds (`build_intel_dp_debug`), so Debug never touches the Release anchor.
  **Diagnosis tell:** if a previously-`max|Δ|=0` gate suddenly shows `area≈3.7e-4` + diffs on EVERY FP
  field, you are running a Debug binary — check `CMAKE_BUILD_TYPE` in the build dir's CMakeCache.txt and
  the binary mtime/size (Debug `fesom_advhordump` ~3.1 MB vs Release ~1.75 MB) BEFORE suspecting a code
  regression. (The real M2.2 gates, run against the restored Release anchor, are all `max|Δ|=0`.)

## L15 — Partial vel_rhs (M2.3): the Debug `-check all` catches what Release + a passing gate hide

M2.3 (`compute_vel_rhs` Coriolis + AB2 + PGF + SSH-gradient, the non-advection part; momentum
advection deferred to M2.4) byte-matched FESOM2 `max|Δ|=0` on the FIRST Release gate run (19 fields).
But the Release gate PASSED **with a latent non-conformance** that only the Debug `-check all` run
exposed — the single most reusable lesson here:

- **`elem2D_nodes` is `(MAX_NV=4, elem2D)` in FESOM3 (quad-capable, L4), `(3, elem2D)` in FESOM2.**
  The faithful-looking transcription `elnodes = mesh%elem2D_nodes(:,elem)` (FESOM2 writes exactly
  that) assigns a 4-extent RHS into the 3-extent `elnodes(3)`. In **Release** ifort silently copies
  the LHS extent (3) — reading `elem2D_nodes(1:3)`, the correct triangle nodes — so the gate is
  `max|Δ|=0`. In **Debug** (`-check all`) it traps: `severe (408) Shape mismatch: extent of dim 1 of
  ELNODES is 3 and ... MESH is 4`. **Fix:** slice explicitly `elem2D_nodes(1:3,elem)` (what `oce_pgf`
  already did). Lesson: a green Release byte-gate does NOT prove conformance; ALWAYS run the Debug
  `-check all` pass (L10/L13/L14 said so for OOB/FPE — add *array-shape* to that list). Any FESOM2
  `(:,elem)` on `elem2D_nodes`/`gradient_sca`/other `MAX_NV`-dim arrays must become `(1:3,elem)` in
  FESOM3.

- **`coriolis` byte-matches via `r2g` transitively — a NEW geometry field gated for free (L9 again).**
  `coriolis(elem)=2·omega·sin(lat_geo)` where `lat_geo` is `r2g` applied to the rotated element
  centroid (`elem_center`). The geometry gate proved `coord_nod2D` (built by `g2r` with the rotation
  matrix), and `r2g` reuses the SAME matrix with the inverse operations — so `coriolis` is `max|Δ|=0`
  by construction, no new transcendental risk. `elem_center` uses the same literal `/3.0_WP` (L7) it
  uses for `elem_cos`. Guard the cartesian/analytic path (no `r2g` on cartesian coords → asin>1 NaN
  under `-fpe0`): for `cartesian=.true.` fill from the stored lat directly (ungated, benign).

- **Two-call gate covers BOTH AB ff branches.** `compute_vel_rhs` has a `save :: lfirst`; the first
  call uses `ff=1.0` (Euler start, `lfirst .and. .not. r_restart`), every later call `ff=ab2=1.6`.
  Driving the REAL routine TWICE in the shim (dump `uv_rhs_eul` then `uv_rhs_ab2`) exercises both —
  they differ over 66% of entries, so the ff selection is non-vacuous. The 2nd call naturally consumes
  the 1st call's `UV_rhsAB` (= Coriolis) as its "previous", a self-consistent step (no reset needed).
  FESOM3 takes `lfirst` as an explicit argument (`.true.` then `.false.`) since it has no `save` state;
  `r_restart` is folded into it (v1 has no restart — M2.11).

- **Forcing the gated config in the shim: three knobs that bite.** (1) `dynamics%ldiag_ke` **defaults
  to `.true.`** in FESOM2 `MOD_DYN` — leave it and `compute_vel_rhs` writes the `ke_*_AB` diag arrays
  (extra work, and FESOM3 has no `ke_*`); force `.false.`. (2) `momadv_opt` is `2` on pi (calls
  `momentum_adv_scalar`); force `0` so the partial gate isolates Coriolis+AB2+PGF (the FESOM3 kernel
  simply omits the momadv call). (3) `dt` is pinned to a shared constant (`1800.0_WP`) in BOTH codes,
  **not** pi's namelist `dt=86400/36=2400 s` — the gate only needs both sides equal in the final
  `dt*(…)/elem_area` scaling, and pinning removes the namelist dependency (L9).

- **The `use_pice=0` path lets a minimal fake `ice` satisfy the real routine.** `compute_vel_rhs`
  unconditionally associates `m_ice => ice%data(2)%values` / `m_snow => ice%data(3)%values` at entry,
  but with `which_ale='linfs'` → `use_pice=0` they are never read. `ice_setup` runs AFTER `ocean_setup`
  (where the shim fires), so no live `ice` exists yet — the shim builds a throwaway `type(t_ice)` with
  just `data(2:3)%values` allocated. The pointer targets exist; the values are irrelevant. (FESOM3 v1
  has no ice at all, so its kernel drops the `p_ice` term entirely — byte-identical since p_ice=0.)

- **Same transitive-gate reason M2.3 "just worked":** every operand was already byte-pinned —
  `coriolis` (r2g/geometry), `elem_area`/`gradient_sca` (geometry gate), `pgf_x`/`pgf_y` (M2.2),
  `ab_epsilon`-derived `ab1=-0.6`/`ab2=1.6` (L12 fold-safe), and the runtime `/elem_area` divisor is
  byte-identical on both sides (L7/L10/L14). The only genuinely new arithmetic is the AB blend and the
  SSH-gradient `sum(gradient_sca·(-g·eta))`, both simple `+`/`*`/`sum` — fold-safe.

## L16 — Momentum advection (M2.4): the full UV_rhs gate; OpenMP-off oracle; gate the intermediate, not the perturbation (PASSED first try)

M2.4 (`momentum_adv_scalar` → the FULL `UV_rhs`: Coriolis + AB2 + PGF + SSH-grad + momentum advection)
byte-matched FESOM2 `max|Δ|=0` on the FIRST Release gate run (21 fields), like all of M1 + M2.1–M2.3.
`momentum_adv_scalar` lives in `src/oce/oce_dyn_velrhs.F90` (private, called by `compute_vel_rhs` when
`momadv_opt==2` — mirroring FESOM2's own `oce_ale_vel_rhs.F90` layout, NOT the plan's separate
`oce_dyn_momadv.F90`). Reusable specifics:

- **The FESOM2 oracle is built `ENABLE_OPENMP=OFF`, so momadv's `omp_set_lock(partit%plock)` +
  `!$OMP ORDERED` paths compile OUT.** momadv is the FIRST gated kernel with OpenMP locks/ordered
  reductions (the M1/M2.3 kernels had `!$OMP DO` but no locks). With OpenMP off they are dead code and
  the per-node edge-scatter runs SERIALLY in natural edge order (1..edge2D) — which is exactly FESOM3's
  serial loop. So the order-sensitive accumulation matches by construction (the L8/L9 argument: edges
  order is geometry-gate-proven, `nod_in_elem2D` order is area-gate-proven). Check `ENABLE_OPENMP` in
  the oracle `build/CMakeCache.txt` before trusting a parallel-reduction kernel's serial transcription;
  if it were ON with `__openmp_reproducible`, the `!$OMP ORDERED` would still serialize, but plain
  OpenMP would reorder → ULP diffs (re-gate would need the reproducible build or 1-thread).

- **Gate the operator's own intermediate, not just its (small) contribution to the dominant field.**
  momadv adds only ~`elem_area·avg(uvnode_rhs)` ≈ 8 to a Coriolis `UV_rhsAB` of ≈1e4 (≈0.07%), and ~1%
  of the final `uv_rhs`. A small perturbation to a dominant field is a WEAK gate *if you only watch the
  dominant field*. Fix: dump `uvnode_rhs` (the post-`areasvol_inv` nodal advection) as its OWN record —
  it is 68.9% non-zero with BOTH signs (the edge-scatter `+nod(1)` / `−nod(2)` branches), gated at its
  native ~1e-6 magnitude. So a bug in EITHER the vertical (`w·du/dz`) OR horizontal (`u·du/dx`) pass
  shows in `uvnode_rhs` directly, independent of how little it moves `uv_rhs`. (The EXACT `max|Δ|=0`
  gate would catch any ULP diff regardless — but the separate intermediate LOCALISES it, the L9/L13
  "gate every intermediate" rule.) And because the gate runs the REAL FESOM2 vertical+horizontal passes,
  FESOM3 matching `uvnode_rhs` proves BOTH passes byte-identical *by construction* (a no-op vertical
  pass would leave `uvnode_rhs` horizontal-only ≠ FESOM2 → FAIL).

- **Prescribe the velocity fully-defined on `1:nl-1` to keep the `0·below-bottom` products clean.** The
  horizontal scatter reads `un2(nz)*UV(:,nz,el2)` for `nz` up to `max(nl1,nl2)` — i.e. BELOW el2's
  bottom, where `un2(nz)=0` (zero-filled) but `UV(:,nz,el2)` is whatever sits there. `0·finite=0`
  (clean) but `0·NaN=NaN`. The driver/shim prescribe `UV` for ALL `nz=1..nl-1` (every element, every
  level), so the below-bottom reads are finite zeros-of-the-analytic-formula, never uninitialised heap →
  no NaN under `-fpe0`. (Same defensive-prescription reason M1.2 set `wvel` on the full column.)

- **`w_e` prescribed like M1.2's `wvel`** — `1e-4·sin(2·lon)·cos(lat)·cos(0.3·nz)`, sign varying in
  space AND depth so the `w·du/dz` finite difference `wu(nz)−wu(nz+1)` is non-trivial. Dumped as an
  input record (`w_e`) so the prescription itself is gated, like `wvel`/`uv_in`.

- **The momadv divisors follow the established rules.** `/(3._WP·hnode(nz,n))` is a RUNTIME divisor but
  byte-identical operands on both sides → reciprocal matches (L7/L10). The vertice→element `/3.0_WP` is
  a LITERAL divisor (exact 1/3 fold, byte-safe; triangles only — the L7 arity caveat). `elem2D_nodes`
  accessed PER-COMPONENT (`(1,el)`/`(2,el)`/`(3,el)`, scalar) not as a `(:,el)` slice → no MAX_NV=4-vs-3
  shape mismatch (the L15 trap; Debug `-check all` confirmed clean).

## L17 — Biharmonic viscosity (M2.4, opt_visc=7): no new geometry; large |du| to exercise the flow-aware max (PASSED first try)

M2.4 biharmonic viscosity (`viscosity_filter(7)` → `visc_filt_bidiff`, FESOM2 `oce_dyn.F90:591-744`)
byte-matched FESOM2 `max|Δ|=0` on the FIRST Release gate run (24 fields now), like all of M1 + M2.1-M2.4.
It lives in `src/oce/oce_dyn_visc.F90` (mirroring FESOM2's `oce_dyn.F90` — a SEPARATE file/operator run
AFTER `compute_vel_rhs`, NOT inside it like momadv). Reusable specifics:

- **opt_visc=7 needs NO new geometry — the HANDOFF's `gradient_vec` worry was unfounded for THIS scheme.**
  `visc_filt_bidiff` is a pure edge-based ∇² applied TWICE (a biharmonic): pass 1 scatters the across-edge
  velocity jump `u1=UV(el1)-UV(el2)` × a flow-aware coeff into the element field `U_c`/`V_c`; pass 2
  scatters the across-edge jump of `U_c` into `UV_rhs/elem_area`. It reads only `edge_tri`, `elem_area`,
  `ulevels`/`nlevels`, `edge2D_in` (+ `dynamics%uv`/`uv_rhs`/`work%u_c`/`v_c`) — every one already gated.
  `gradient_vec` is still deferred but is NOT required here (it would be for a vector-Laplacian scheme; the
  edge-difference biharmonic sidesteps it). So the L9 transitive-gate pattern again → `max|Δ|=0` first run.

- **First gated kernel to use `edge2D_in` (interior-edge-only loop; free slip on boundary edges).** FESOM2
  `if(myList_edge2D(ed)>edge2D_in) cycle`. At 1-rank `myList_edge2D` is identity AND the mesh orders interior
  edges first (`1..edge2D_in`, boundary `edge2D_in+1..edge2D` — `fvom_init` writes them so, line 554-555), so
  the lift is `do ed=1,edge2D; if(ed>edge2D_in) cycle`. `edge2D_in` is byte-pinned transitively: same
  `edgenum.out` (FESOM3 reads line 2), same edge ordering proven by the geometry gate (`edge_tri`/
  `edge_cross_dxdy` were `max|Δ|=0`). For interior edges BOTH `edge_tri(1/2,ed)>0`, so no `el2<0` guard is
  needed (unlike momadv, which looped ALL edges). `edge_tri(:,ed)` is the size-2 left/right slot — NOT a
  MAX_NV dim, so no L15 shape trap, and `elem2D_nodes` is never touched. Debug `-check all` EXIT 0.

- **The biharmonic is TWO SEPARATE edge loops with `U_c`/`V_c` between — the serial structure replaces the
  `!$OMP BARRIER`.** Pass 1 must FULLY fill `U_c` before pass 2 reads `U_c(el1)-U_c(el2)`; FESOM2 enforces
  this with an `!$OMP BARRIER` between the loops, FESOM3's two sequential `do ed` loops get it for free.
  `exchange_elem(U_c/V_c)` between them is a 1-rank no-op (lifted at M2.12). `U_c`/`V_c` are zeroed over all
  elements first. The OpenMP-off oracle (L16, re-confirmed `ENABLE_OPENMP=OFF` in the oracle CMakeCache)
  compiles out the `omp_set_lock`/`!$OMP ORDERED` scatter → serial edge order → matches FESOM3's serial loop.

- **Gate STRENGTH: the flow-aware coeff `max(γ0, max(γ1·|du|, γ2·|du|²))` needs LARGE |du|, or only γ0 fires.**
  The pi gammas are `γ0=0.003` (background; OVERRIDES the `t_dyn` default 0.03 — pin the namelist value!),
  `γ1=0.1`, `γ2=0.285`, `γ0_h=γ1_h=0` (pure biharmonic). The winners cross at `|du|=γ0/γ1=0.03` (γ0→γ1) and
  `|du|=γ1/γ2=0.351` (γ1→γ2). At pi-realistic velocities (~0.1 m/s) `|du|≲0.02` so ONLY γ0 (background) is
  ever selected → the γ1/γ2 branches are computed-but-not-selected, masking a transcription bug there (the
  L11 weak-gate trap). Fix: bump the SHARED prescribed `UV` amplitude to 2.0/1.5 m/s (a strong-current
  stress test; the operator gate has no stability constraint) so `|du|` spans all three — measured
  selection **γ0/γ1/γ2 = 20.4/78.6/1.0%** of edge-levels (γ2 on ~2000 levels, a clear margin). A
  driver-side, NON-gated diagnostic counts the actual `max` winner to verify this. NOTE: γ2 is a near-dead
  path in PRODUCTION pi (the namelist even labels γ2 "only used for opt_visc=5/8") — it fires only on
  synthetic |du|>0.351; the gate exercises it deliberately. Bumping UV re-gated M2.3/M2.4 with new values
  (automatic — identical formula on both sides; all prior fields stayed `max|Δ|=0`).

- **Gate the operator's own intermediate (the L16 rule).** Dump `visc_u_c`/`visc_v_c` (the first-stage
  Laplacian, the element field pass 1 builds) as their OWN records, so a pass-1 bug localises there while a
  pass-2 bug shows in `uv_rhs_visc`-but-not-`visc_u_c`. The post-viscosity `uv_rhs_visc` is the gate target;
  viscosity only READS `UV` and ADDS its increment into the incoming `uv_rhs` (= `uv_rhs_ab2`, already gated),
  so running it after the 2nd `compute_vel_rhs` is well-defined.

- **The viLapl (harmonic-addition) term is transcribed but contributes 0 on pi (γ_h=0) — a deferred sub-gate.**
  `viLapl=dt·max(γ0_h, γ1_h·|du|)·len`; with γ0_h=γ1_h=0 it is exactly 0 (no FPE: `0·sqrt(...)`=0). The
  γ_h>0 combined harmonic+biharmonic path is ungated on pi (like the cavity/zstar/use_density_ref deferrals);
  re-gate it when a config with γ_h>0 is exercised. All divisors follow the established rules: `len=
  sqrt(sum(elem_area))`, `-dt·sqrt(...)`, final `/elem_area(el)` are RUNTIME divisors byte-identical on both
  sides (L7/L10/L14). Local intermediates are `real(WP)`; FESOM2 hard-codes them `real(kind=8)` — byte-equal
  at the DP anchor (WP=8), a deferred SP-precision nuance.

## L18 — Implicit vertical viscosity (M2.5): a sequential TDMA byte-matches by pure transitivity (PASSED first try)

M2.5 implicit vertical viscosity (`impl_vert_visc_ale`, FESOM2 `oce_ale.F90:3315-3526`) byte-matched FESOM2
`max|Δ|=0` on the FIRST Release gate run (28 fields now), like all of M1 + M2.1-M2.4. It lives in a new
`src/oce/oce_dyn_ivertvisc.F90`, a SEPARATE operator run AFTER `viscosity_filter` (FESOM2 `oce_ale.F90:3874`).
Reusable specifics:

- **A tridiagonal Thomas solve is a strictly SEQUENTIAL recurrence — there is NO summation/scatter order
  ambiguity to defend (unlike the edge/node scatters of M2.4/M2.4-visc).** Per element the column is
  independent; forward elimination then back substitution are recurrences with one well-defined order. So the
  gate reduces ENTIRELY to "are the operands byte-identical?" — and they are: `UV`, the incoming `uv_rhs` (=
  the M2.4-visc-gated `uv_rhs_visc`), the prescribed `w_i`/`Av`/`stress_surf`, the geometry-proven
  `helem`/`zbar_e_bot`/`ulevels`/`nlevels`, and `C_d`/`density_0`/`dt`. The L9 transitive-gate pattern again →
  `max|Δ|=0` first run. The many runtime divisors (`zinv`, the Z/zbar differences, the Thomas `1/b`, `1/m`)
  are byte-identical operands on both sides, so the `-no-prec-div` reciprocal matches (L7/L10/L14).

- **Inputs not yet produced by the core (`Av` from PP mixing M2.8, `stress_surf` from forcing M2.10) → pass
  as EXPLICIT dummy args, not as type fields.** Honest about provenance, minimal diff, defers the
  type-placement decision; when M2.8/M2.10 land the caller just sources them (`dyn%work%Av`, `forcing%...`) and
  the kernel is unchanged. FESOM2 pulls both from `o_ARRAYS`; the shim sets those module arrays directly.
  Prescribe `Av` STRICTLY POSITIVE (a real viscosity: `[1e-3,1.4e-2]` m²/s here) over all `nz=1..nl` (the TDMA
  reads `Av(nzmax)`), and `stress_surf` sign-varying (~0.1 N/m², both signs → the surface BC is exercised
  fully).

- **New mesh field `zbar_e_bot` + the (missing) `helem` — both full-cell trivial; serialization is
  unallocated-safe.** The pressure driver built `hnode` but NOT `helem`; M2.5 needs both. Full cells (pi,
  `use_partial_cell=.false.`): `helem(nz,e)=zbar(nz)-zbar(nz+1)` (the element analog of `hnode`),
  `zbar_e_bot(e)=zbar(nlevels(e))` (FESOM2 `init_bottom_elem_thickness`). The kernel rebuilds the per-column
  `zbar_n`/`Z_n` bottom-up from these. Added `zbar_e_bot` to `t_mesh` + write/read serialization;
  `write_bin_array` writes `0` for an unallocated array and `read_bin_array` leaves it unallocated, so the
  ctest mesh round-trip (where `zbar_e_bot` is never built) stays green — 13/13.

- **The single-layer-column OOB is a REAL benign-OOB in FESOM2 that pi avoids — CHECK THE MESH (`elvls.out`)
  before assuming `-check all` clean.** For `nlevels(elem)==2` the "last row" reads `Z_n(nzmax-2)=Z_n(0)` and
  `UV(:,nzmax-2)=UV(:,0)` — values that are overwritten by the "first row" block / multiplied by `a(top)=0`,
  so the Release result is byte-correct (FESOM2's Release oracle tolerates it), but `-check all` would TRAP.
  pi has NO such columns (min element `nlevels=5` per `elvls.out` — i.e. ≥4 layers), so the EXACT transcription
  is both byte-identical AND `-check all` clean (EXIT 0, confirmed). A `nlevels==2` shelf column (CORE2) needs
  a guard matched to FESOM2's behaviour — deferred to M2.11. Determining this up front (read `elvls.out`)
  avoided a guess.

- **`w_i` (implicit vertical velocity) prescribed sign-varying to exercise both upwind branches; the TDMA must
  be shown NON-vacuous.** A DISTINCT analytic formula from `w_e` (`sin(lon)·cos(2·lat)·cos(0.25·nz)`) so its
  sign varies in space AND depth → `min(0,wu)`/`max(0,wu)` (+ the `wd` pair) all fire. Driver diagnostic
  (non-gated): `w_i>0/<0 = 75360/75360` (both signs) and `max|d(uv_rhs)|=0.985` (the solve substantially
  transforms the rhs, not a near-identity) — the L11 weak-gate guard.

- **No auto-generated `*_interface` module for `impl_vert_visc_ale` (only the `_vtransp` variant) → declare an
  explicit interface block in the shim.** FESOM2's interface-gen tool skips it because it is declared in an
  explicit interface block inside `oce_ale.F90:269`; the other gated routines (`pressure_bv`,
  `compute_vel_rhs`, `visc_filt_bidiff`) DO get `<name>_interface.mod`. So the shim copies the interface block
  verbatim. Pin `C_d=0.0025` in the shim (the routine reads it from `o_PARAM`; matches the FESOM3
  `mod_param_phys` default).

- **Gate target is `UV_rhs`, NOT `UV`.** The routine OVERWRITES `UV_rhs` with the Thomas solution; `UV` is
  read-only (the velocity update `UV += UV_rhs` happens later in the timestep). The HANDOFF "post-solve UV"
  note was imprecise — the dumped gate target is `uv_rhs_ivv`. Out-of-pi branches dropped: `ldiag_ke`
  `ke_wind`/`ke_drag` (no `ke_*` in FESOM3 `t_dyn`, as `compute_vel_rhs`), and the `toy_ocean`
  dbgyre/neverworld2/soufflet constant-`C_d` friction (pi `toy_ocean=.false.` → quadratic bottom drag
  `-C_d·|UV_bottom|`).
