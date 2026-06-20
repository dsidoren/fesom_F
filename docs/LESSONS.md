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

## L19 — SSH stiffness + CG solve (M2.6): the first iterative solver byte-matches; the incremental-build ULP-drift footgun (PASSED after a clean rebuild)

M2.6 SSH (CSR stiffness `init_stiff_mat_ale` + `compute_ssh_rhs_ale` + preconditioned-CG `solve_ssh_ale`,
FESOM2 `oce_ale.F90:1584/2012/3272` + `solver.F90`) byte-matched FESOM2 `max|Δ|=0` on all 4 new fields
(`ssh_stiff_diag`, `ssh_Aeta`, `ssh_rhs`, `d_eta` → **32 fields** total). Built `src/oce/oce_ssh_rhs.F90`
(stiffness + rhs) + `src/oce/oce_ssh_solve.F90` (precond + CG). The FIRST iterative solver in the port.

- **THE FOOTGUN (cost me one red gate): after ADDING NEW module files, an incremental FESOM3 build is
  byte-UNRELIABLE — clean-rebuild before gating.** The first gate run FAILED with a *uniform ~few-ULP
  relative* delta on EVERY field — including prescribed trig-derived inputs (`temp` 7e-15, `coriolis` 5e-20)
  while integer/geometry fields (`zbar`, `hnode`) stayed exactly 0. That split is the L10 "uniform flag/
  codegen effect" signature, NOT a kernel bug. Root cause: CMake's `GLOB CONFIGURE_DEPENDS` detected the 2
  new `oce_ssh_*.F90`, re-ran `cmake` configure mid-incremental-build, and the partial recompile linked
  objects built against a MIX of stale/fresh `.mod` interfaces → ULP-drifted FP codegen vs the FESOM2
  oracle. `./configure.sh --clean --build` (full recompile, consistent `.mod` set) restored `max|Δ|=0` on
  all 32 fields on the next run. **Rule: when a milestone adds NEW source files (not just edits), do a CLEAN
  Release rebuild before trusting the byte-gate.** (Edits-only incremental builds stayed byte-exact through
  M2.1–M2.5; it is specifically NEW files + the reconfigure that drift.)

- **A CG solve byte-matches by the same L9 transitivity as the M2.5 TDMA — given byte-identical operands the
  iteration is a deterministic recurrence.** Same matrix + rhs + x0 ⇒ identical per-iteration scalars ⇒ the
  convergence test (`sqrt(Σrr²/nod2D) < soltol·sqrt(Σb²/nod2D)`, soltol=1e-5) fires at the SAME iteration
  (here **37**) ⇒ identical `d_eta`. So a byte-matching `d_eta` proves the entire 37-iteration recurrence
  (matvecs + dot-products + preconditioner) is bit-identical. The operands are all pre-gated: `ssh_stiff`
  (geometry-proven `edge_tri`/`edges`/`gradient_sca`/`edge_cross_dxdy`/`zbar_e_bot`/`areasvol` + the edge-
  scatter accumulation order), `ssh_rhs` (same edge order + the M2.5-gated `uv_rhs`/`helem`), `g`/`dt`/`α`/`θ`.

- **REDUCTION FORM is part of the bits.** The FESOM2 oracle is `ENABLE_OPENMP=OFF` and `__openmp_reproducible`
  undefined, so its dot-products compile to the explicit serial `DO row; s=s+…; END DO` (the `!$OMP REDUCTION`
  is an inert comment). Transcribe THAT serial form for `s_old`/`s_aux`/`sprod` — NOT a `sum()` intrinsic
  (the compiler may reduce it in a different order). The matrix-vector products DO use `sum()` over a CSR
  slice, exactly as FESOM2 — both codes `sum()` the same slice ⇒ byte-identical. (Same OpenMP-off reasoning
  as L16/L17's serial edge scatter.)

- **dt for the stiffness is the pi NAMELIST dt (=86400/36=2400 s), NOT the shim's `dt=1800`, and must be
  computed the SAME way (not a `2400.0` literal).** `init_stiff_mat_ale` runs at `ocean_setup` (oce_setup_step
  :140) where dt is still the namelist value, BEFORE the shim overrides it to 1800 for the vel_rhs gates. dt
  enters M2.6 ONLY via the stiffness (`factor=g·dt·α·θ`, mass `areasvol/dt`); compute_ssh_rhs has no dt. So
  the FESOM3 driver builds the stiffness with `dt_ssh = 86400._WP/real(36,WP)` (FESOM2's `gen_model_setup`
  formula) — a literal `2400.0_WP` could differ by an ULP under `-no-prec-div` (the L7 trap) and break the
  whole matrix.

- **`α=θ=1.0` on pi (FESOM2 default, namelist doesn't set them) ⇒ two simplifications:** the linfs water-flux
  term `(1-α)·ssh_rhs_old` VANISHES (so `ssh_rhs_old` is multiplied by 0 — set it to 0, don't bother
  prescribing it), and `factor = g·dt`. Force `α=θ=1.0` in BOTH shim and driver defensively (the stiffness
  already used 1.0; compute_ssh_rhs re-reads them).

- **linfs builds the stiffness ONCE and never updates it** (`oce_ale.F90:3921` calls `update_stiff_mat_ale`
  only for NON-linfs). So the single `init_stiff_mat_ale` assembly with the unperturbed depth IS the
  production matrix on pi/reduced-M2 — no `update_stiff_mat_ale` port needed for M2.

- **1-rank collapses the global-numbering machinery to the identity.** The CG uses `colind_loc`/`rowptr_loc`
  (LOCAL CSR); the FESOM2 global remap (`rpart.out` mapping, the per-PE nza offset, `myList_nod2D`) is
  identity at 1-rank and is dropped — build the local CSR and set `colind=colind_loc`, `rowptr=rowptr_loc`.
  `exchange_nod`(diag/rr/pp/x) and the `MPI_Allreduce`(s_old/s_aux/sprod) are no-ops (the local serial sum
  IS the global sum), dropped. The new multi-rank bit-identity risk — the cross-rank dot-product reduction
  ORDER — is an M2.12 concern.

- **Gate the matrix with a MATVEC, not just the diagonal.** `ssh_stiff_diag` (the row-diagonal = mass +
  self-stiffness) localises the assembly, but `ssh_Aeta = A·eta_n` (full CSR matvec against the prescribed
  `eta_n`) exercises EVERY non-zero — the off-diagonals are the bulk of the operator. Both `max|Δ|=0`.

- **The `Σ ssh_rhs` telescoping check is ~1e-13 RELATIVE, not absolute.** `ssh_rhs` is an edge-divergence
  (`+/-（c1+c2)` scattered into the two edge nodes) so `Σ_nodes` telescopes to 0 in exact arithmetic; the FP
  roundoff scales with the `ssh_rhs` magnitude. With UV bumped to 2.0/1.5 m/s (the M2.4-visc strength test)
  `ssh_rhs ~ 1e8`, so `Σ = -3.7e-5` (≈1e-13 relative) — consistent, non-vacuous. The plan's "~1e-13" was the
  small-UV estimate.

- **The prescribe-and-stop shim STILL works for M2.6** (the worry it might need to "run past forcing" was
  unfounded): `init_stiff_mat_ale` runs at `ocean_setup` line 140, BEFORE the end-of-`ocean_setup` shim, so
  the matrix is already assembled; the shim just prescribes `UV`/post-TDMA `UV_rhs`/`d_eta`=0/`ssh_rhs_old`=0
  and drives the REAL `compute_ssh_rhs_ale` + `solve_ssh_ale`. The reduced-M2 namelist is still NOT needed.
  No auto-gen `*_interface` module for `compute_ssh_rhs_ale`/`solve_ssh_ale` (internal interface blocks) →
  declare explicit interfaces in the shim (as for `impl_vert_visc_ale`, L18).

- **`zbar_e_srf` = `zbar(ulevels(elem))` = `zbar(1)` = 0 on pi (no cavity).** init_stiff's H factor is
  `(zbar_e_bot − zbar_e_srf)`; computed inline from `zbar`+`ulevels` (no new mesh field — faithful to the
  FESOM2 non-cavity default `oce_ale.F90:525`). The CG `nod2D` divisor (`rtol`, exit test) is the global
  `nod2D=3140`. Debug `-check all` clean (RC=0): `n_pos(12,nod2D)` is wide enough (pi max node degree < 11),
  no CSR-slice OOB.

## L20 — ALE velocity/SSH/thickness-W update (M2.7): the post-CG tail byte-matches; in-place overwrite of prescribed inputs; the Debug I/O stack-temp overflow (PASSED first try)

M2.7 (the linfs post-CG tail: `update_vel` + `compute_hbar_ale` + the `eta_n` blend + `vert_vel_ale`
→ `compute_CFLz` + `compute_Wvel_split`) byte-matched FESOM2 `max|Δ|=0` on the FIRST Release gate run
(all 11 new fields → **43 fields**), like all of M1 + M2.1-M2.6. Built `src/oce/oce_ale.F90` (FESOM3
consolidates `update_vel` — FESOM2 keeps it in `oce_dyn.F90` — with the `oce_ale.F90` ALE routines into
one M2.7 module). Reusable specifics:

- **The whole post-CG chain is `max|Δ|=0` by pure L9 transitivity — every operand was already pinned.**
  `update_vel` is a `gradient_sca` contraction of `-g·θ·dt·d_eta` (d_eta from the M2.6 CG, gradient_sca
  geometry-gated) added to the post-TDMA `UV_rhs` (M2.5; the CG solve never touches it) into `UV`. The
  `compute_hbar_ale` + `vert_vel_ale` edge-divergences reuse the SAME edge order + `helem`/`edge_cross_dxdy`/
  `areasvol`/`area` the M2.6 `compute_ssh_rhs_ale` already gated. So nothing genuinely new arithmetic-wise →
  faithful transcription "just works" (the recurring L9 pattern). `hbar` is prescribed identically on both
  sides. `dt=1800`/`θ=1`/`α=1` pinned (as M2.3-M2.6); the OpenMP-off oracle (L16/L17/L19) serial edge order
  matches the serial loops.

- **Three prescribed inputs are OVERWRITTEN in place by the M2.7 kernels — save copies BEFORE the chain or
  the INPUT dump records become wrong.** `update_vel` overwrites `UV` (but `uv_in` was already saved at
  M2.3); the `eta_n` blend overwrites `dynamics%eta_n`; `compute_Wvel_split` (inside `vert_vel_ale`)
  overwrites `dynamics%w_e`/`w_i` — the SAME arrays prescribed as the M2.4/M2.5 inputs. Since the dump fires
  at the END (after all compute), the `eta_n`/`w_e`/`w_i` INPUT records would otherwise dump the M2.7 output.
  Fix: save `eta_n_in`/`w_e_in`/`w_i_in` right before the chain (they are fully consumed by M2.4/M2.5 first),
  dump those for the input records, and dump the post-M2.7 values as NEW records (`eta_n_upd`/`w_split_e`/
  `w_split_i`). This is correct physics: in the real timestep `vert_vel_ale` recomputes `w_e`/`w_i` for the
  NEXT step's momadv/ivertvisc. Same save-the-input pattern as `uv_in` (M2.3) and the M2.1 raw/smoothed split.

- **`compute_CFLz` keeps its TWO-statement form for byte-reproducibility (the L19 "reduction form is part of
  the bits" rule, made explicit by FESOM2).** `CFL_z(nz)=CFL_z(nz)+c1` then `CFL_z(nz+1)=c2` — FESOM2's own
  comment says this exact split (vs. folding both into one accumulate) is "for the sake of reproducibility …
  (rounding error)". Transcribe it verbatim; do NOT combine. `c1`/`c2` are scalars here (a DIFFERENT `c1`
  from the `vert_vel_ale` per-level array).

- **linfs collapses `vert_vel_ale` to W-only.** `which_ale='linfs'` ⇒ the `zlevel`/`zstar` thickness-
  redistribution branches are NOT taken, so `hnode_new` stays = `hnode` (its init value — gate it anyway to
  confirm the linfs path leaves it untouched) and there is no surface Wvel/hnode correction; `compute_hbar_ale`'s
  water-flux term (`.not. linfs`) vanishes. So `vert_vel_ale` is just: zero W → edge-scatter `div(UV·h)` →
  cumsum bottom-up → `/area`. Fer_GM (`fer_UV`/`fer_Wvel`) and ldiag_ke (`ke_*`) branches dropped (no GM / no
  `ke_*` in v1, as `compute_vel_rhs`/`impl_vert_visc_ale`). The full-free-surface thickness evolution gets its
  own later gate.

- **Non-vacuity needs `use_wsplit=.true.` (pi production) + a large UV — else `compute_Wvel_split` is the
  trivial `Wvel_e=Wvel` branch.** The split only deviates where `CFL_z > wsplit_maxcfl` (=1.0 on pi). With the
  M2.4-visc-bumped UV (2.0/1.5 m/s) the W divergence is large enough that `CFL_z>1` on **13253** (nz,node) →
  the `dd`/`Wvel_e`/`Wvel_i` formula is genuinely exercised (a driver diagnostic counts it; the L11/L17 weak-
  gate guard). Force `use_wsplit=.true.`/`wsplit_maxcfl=1.0` in BOTH shim and driver (the namelist sets them,
  but assert — defensive like the M2.6 α=θ=1). `max|w|=0.042` m/s (physical), `max|uv_upd|=2.5` (the SSH-grad
  correction + UV_rhs visibly moved UV) — all non-vacuous.

- **The "`K_v⁻` deformation bound" plan bullet was a MISLABEL — there is no Kv bound in the post-CG ALE path.**
  `Kv` (vertical diffusivity) is produced by PP mixing (M2.8), not the velocity/SSH/thickness update. Grep of
  the `oce_ale.F90` step tail (`update_vel`→`compute_hbar_ale`→`eta_n`→`vert_vel_ale`) confirmed no Kv/
  deformation bound there. Folded the bullet into M2.8.

- **THE DEBUG FOOTGUN (not a kernel bug): the I/O dump writer SEGFAULTs under Debug `-check all` on the default
  8 MB stack — `ulimit -s unlimited` fixes it.** The compute ran CLEAN under `-check all` (every M2.7 diagnostic
  printed → all six kernels completed, no OOB/shape/FPE; the L15 `elem2D_nodes(1:3,·)` slices are correct). The
  crash was at `mod_advhor_dump.F90:68` `write(u) real(a, real64)` on the FIRST big `wr_r3` (`uv_in`, an
  EXISTING record) — ifort `-O0` puts the ~4.4 MB `real(...,MP)` array temporary on the stack, overflowing the
  8192 KB default (the extra M2.7 records pushed cumulative pressure over the edge). `ulimit -s unlimited` →
  EXIT 0, dump written. Release (`-O3`) is unaffected (different temp handling). **Diagnosis tell:** a segfault
  in `wr_r{2,3}`/`real(a,real64)` AFTER all compute diagnostics printed is a stack-temp overflow, NOT a kernel
  fault — raise the stack, don't hunt the kernel. (Byte-gate is Release-vs-Release; Debug is OOB-only, L10.)

## L21 — PP vertical mixing (M2.8): three sequential passes byte-match by transitivity; prescribe `uvnode` (not `UV`) to isolate the add + control the Ri-factor range (PASSED first gate run)

M2.8 (Pacanowski-Philander Richardson-number mixing: `oce_mixing_pp` → `Kv` on nodes, `Av` on elements)
byte-matched FESOM2 `max|Δ|=0` on the FIRST Release gate run (3 new fields `uvnode`/`pp_Kv`/`pp_Av` →
**46 fields**), like all of M1 + M2.1-M2.7. Built `src/oce/oce_ale_mixing_pp.F90` (`oce_mixing_pp` +
`Kv0_background_qiang` + `Kv0_background`, mirroring the oracle filename). Reusable specifics:

- **Three SEQUENTIAL passes, byte-match by pure L9 transitivity.** Pass 1 (nodes): `Kv := factor =
  shear/(shear + 5·max(N²,0) + 1e-14)`, the inverse-Richardson factor `1/(1+5·Ri)` with `shear =
  |d(uvnode)/dz|²`. Pass 2 (elements): `Av = mix_coeff_PP·mean₃(factor²) + A_ver`. Pass 3 (nodes): `Kv =
  mix_coeff_PP·factor³ + K_ver`. The ordering is LOAD-BEARING — `Kv` is the scratch that holds the factor
  (pass 1) AND the final diffusivity (pass 3), so the `Av` elem loop MUST run between them (Av reads factor²,
  Kv overwrites with factor³). No scatter/accumulation (`sum(Kv(nz,elnodes)**2)` is a fixed 3-element array
  order, elnodes geom-gated) → no order ambiguity → every operand pinned (`bvfreq` M2.1, `uvnode`/`Z_3d_n`
  geometry) ⇒ `max|Δ|=0` first run. The recurring L9 pattern (like M2.5's TDMA).

- **Prescribe `dyn%uvnode` DIRECTLY (a new input), do NOT compute it from the prescribed `UV`.** PP reads the
  nodal velocity `dynamics%uvnode` (the area-weighted elem→node average `compute_vel_nodes` fills each step,
  `oce_dyn.F90:177`). Two wins from prescribing it rather than running `compute_vel_nodes`: (1) **isolation** —
  the existing `UV` (and hence all M2.3-M2.7 records) is UNTOUCHED, so M2.8 is a pure additive gate (no
  re-gate churn); (2) **non-vacuity control** — the PP factor depends on the VERTICAL shear `Δuvnode/Δz`, and
  the existing `UV`'s depth term (`-0.005·nz`, tuned so it cancels in the across-EDGE viscosity `du`) gives a
  shear ~3e-8 → factor ~1e-3 (technically non-zero in dp, but a weak gate). A bespoke `uvnode` with a strong
  vertical shear (`1.2·cos(lat)sin(lon)·sin(0.5·nz)` etc., per-layer jump ~0.6 m/s) makes the factor span
  **[5.8e-7, 0.93]** (62.4% of node-levels > 0.1) → `f²`/`f³` exercised across their full range, `max|Kv|`=8e-3
  (800× the `K_ver`=1e-5 background), `max|Av|`=8.7e-3 (87× `A_ver`). This is the M2.5 precedent (prescribe the
  not-yet-sourced input: `Av`/`stress_surf` there): `compute_vel_nodes` is gated with the step at M2.9, where
  `impl_vert_visc_ale` will also source its `Av` from `dyn%work%Av` instead of the prescribed analytic one.

- **`A_ver` is a pi NAMELIST override (1e-4), NOT the `mod_param_phys` default (1e-3) — force it.** The pi
  `namelist.oce` sets `A_ver=1.e-4`; FESOM3's `mod_param_phys` default is `0.001`. PP's `Av=…+A_ver` reads it,
  so the driver+shim must FORCE `A_ver=1e-4` (+ `mix_coeff_PP=0.01`, `K_ver=1e-5`, `Kv0_const=.true.`, all pi
  defaults). Proof `A_ver` was unused by M2.1-M2.7: M2.5 gated `max|Δ|=0` with the oracle at `A_ver=1e-4` and
  FESOM3 at its 1e-3 default — if `impl_vert_visc_ale` read `A_ver` they would have diverged. (`Kv0_const=.false.`
  Qiang lat/depth background + the cavity `nzmin>1` path are transcribed but UNGATED on pi — M2.11.)

- **PP OVERWRITES the shared `Av`; save the M2.5 prescribed `Av` first (the L20 save-the-input pattern).** In
  the FESOM2 shim `Av`/`Kv` are o_ARRAYS globals: PP clobbers the SAME `Av` the M2.5 section prescribed, and the
  dump fires at the END → the `Av` record would dump PP's output. Save `Av_in = Av` before PP, dump `Av_in` for
  the M2.5 record, dump post-PP `Av`/`Kv` as NEW records `pp_Av`/`pp_Kv`. (FESOM3 has no conflict: PP writes
  `dyn%work%Av`, distinct from the driver's local M2.5 `Av` — but the DUMPED values match.) Both sides pre-zero
  `Kv`/`Av` before PP (it writes only `nz∈[nzmin+1,nzmax-1]`; surface/bottom/below stay a deterministic 0).

- **`target` on the `dyn` dummy is required** for the `UVnode=>dyn%uvnode` / `Kv=>dyn%work%Kv` pointer aliases
  (ifort `#6796`) — FESOM2 declares `dynamics` `target` too; mirror it. The only compile error; otherwise the
  transcription + Debug `-check all` were clean first try (the L15 `elem2D_nodes(1:3,elem)` slice avoids the
  MAX_NV=4 trap; FESOM2's `elem2D_nodes(:,elem)` is `(3,·)` on its runtime path so its `:` is already 1:3).

## L22 — Convective adjustment (M2.8b): the FIRST gate to deliberately perturb the shared T/S (an unstable band) for non-vacuity — the cascade auto-re-verifies (PASSED first gate run)

M2.8b (`mo_convect`: where N²<0, floor `Kv` nodes / `Av` elements to `instabmix_kv`=0.1 — the convective
overturning proxy, run after PP, FESOM2 `oce_ale.F90:3729`) byte-matched FESOM2 `max|Δ|=0` on the FIRST gate
run (2 new fields `moc_Kv`/`moc_Av` → **48 fields**). Built `src/oce/oce_mo_conv.F90` (`mo_convect`). Specifics:

- **The convective branch is UNTESTABLE on a stable T/S — perturb the shared input to make `bvfreq<0`, and the
  cascade re-verifies for free.** The standing prescribed T/S (T↓, S↑ with depth) is statically STABLE
  everywhere → `bvfreq≥0` → the `if (bvfreq<0) Kv=max(Kv,0.1)` floor NEVER fires → a vacuous gate (passes
  trivially, a transcription bug in the condition/`max` uncaught). Fix: add a LOCALIZED unstable band to T — a
  warm subsurface lens `+8·max(0,cos lat·cos lon)·min(nz,8)/8` whose downward warming over the top 8 levels
  overcomes the −0.20·nz(T)/+0.03·nz(S) stabilisation in the warm hemisphere → `bvfreq<0` on **5165 node-levels
  / 9936 elem-levels**, the floor genuinely fires (`max|ΔKv|`=0.090, `max|ΔAv|`=0.096, i.e. ~0.01→0.1). This is
  the FIRST gate to deliberately change the shared T/S — it CASCADES (T→density→hpressure→pgf→vel_rhs→…→every
  downstream field, 40+ records), but **all M2.1-M2.8 records re-verify `max|Δ|=0` automatically** (both sides
  use the identical T) — the cascade is harmless, only the new `moc_Kv`/`moc_Av` are added. Keep it BOUNDED
  (≤+8 °C → T≤~28 °C) and S>0 so the EOS `sqrt(s)` stays well-posed; verify no NaN. (Side effect, noted not
  re-committed: the unstable region clamps PP's `5·max(N²,0)=0` → the M2.8 Ri factor now reaches exactly **1.0**,
  vs the [5.8e-7, 0.93] of the M2.8 commit's stable T/S — `pp_Kv`/`pp_Av` re-verify regardless.)

- **`mo_convect` overwrites `Kv`/`Av` in place → save the PP output before calling it (the L20/L21 pattern).**
  `pp_Kv`/`pp_Av` (the M2.8 records) must dump the PRE-adjustment PP output, but the dump fires at the end after
  `mo_convect` has floored `Kv`/`Av`. Save `pp_Kv_save`/`pp_Av_save` right after PP, dump those for `pp_Kv`/
  `pp_Av` and the live post-adjustment `Kv`/`Av` as NEW records `moc_Kv`/`moc_Av`. The difference localises the
  gate: `pp_*` isolates PP, `moc_*` isolates the convective floor.

- **Port only the reachable branch; OMIT (don't stub) the forcing-coupled `use_momix` block.** FESOM2
  `mo_convect` has three enhancements: `use_instabmix` (convective, M2-live), `use_windmix` (near-surface, uses
  only params+`Kv`/`Av` → transcribe guarded-off) and `use_momix` (TB04 Monin-Obukhov). The momix block reads
  forcing/ice fields not ported yet (`water_flux`/`heat_flux`/`stress_node_surf`/ice/`mo`/`mixlength`/`mo_length`/
  `pmlktmo`) — they don't EXIST, so it cannot compile → OMIT it entirely (a comment marks it deferred to M2.10),
  not a dead `if`-guard. Force `use_momix=.false.` in the gate. **Oracle-side gotcha:** the FESOM2 `mo_convect`
  signature is `(ice, partit, mesh)` and it pointer-assigns `u_ice=>ice%uice` / `v_ice=>ice%vice` /
  `a_ice=>ice%data(1)%values` UNCONDITIONALLY (before the `use_momix` guard), so the shim's `ice_dummy` must have
  those three allocated (size `nod2D`) even though `use_momix=.false.` never dereferences them — extend the M2.3
  `ice_dummy` (which only had `data(2)`/`data(3)` for `compute_vel_rhs`'s `m_ice`/`m_snow`). FESOM3's ported
  `mo_convect(dyn, mesh)` drops `ice` entirely (no momix → no ice).

## L23 — Tracer-solve assembly (M2.9a): the FIRST `Kv` consumer (sourced LIVE), `NaN·0` poisoning from Redi-off slope reads, module-procedure vs free-subroutine symbol mangling (PASSED first gate run)

M2.9a (the tracer diffusion solve: `diff_tracers_ale` = `diff_part_hor_redi` horizontal diffusion + ALE
reconstruct + `diff_ver_part_impl_ale` implicit vertical-diffusion TDMA + `bc_surface`) byte-matched FESOM2
`max|Δ|=0` on the FIRST gate run (9 new `tsol_*` records → **57 fields**). Built into `src/oce/oce_ale_tracer.F90`
(mirrors FESOM2's file, NOT the plan's `oce_solve_tracers.F90` — the L16 layout precedent). Specifics:

- **The TDMA is the FIRST consumer of a PP output — source `dyn%work%Kv` LIVE (post-`mo_convect`), don't
  prescribe it.** Every prior M2 gate PRESCRIBED its diffusivity/viscosity (M2.5 `Av`, M2.8 `uvnode`→`Kv`); here
  `diff_ver_part_impl_ale` reads the `Kv` that PP+`mo_convect` already wrote into `o_ARRAYS Kv` / FESOM3
  `dyn%work%Kv` earlier in the SAME gate. Append M2.9a after the M2.8b section so `Kv` is live; both sides consume
  `Kv∈[0, 0.1]` (the `mo_convect` floor) — non-vacuous AND it realises the integration the step (M2.9b) needs.
  Gate `max|dT|`=1.26 °C / `max|dS|`=0.11 confirms the solve actually moves T/S.

- **On pi the tracer-solve diffusion reduces to the implicit vertical TDMA ALONE; transcribe THAT path, defer the
  rest.** `K_hor=0` (horizontal diffusion `diff_part_hor_redi` is a no-op in production), `i_vert_diff=.true.`
  (skip the explicit `diff_ver_part_expl_ale`), T/S use `'FCT'` → `do_wimpl=.false.` (the implicit vertical
  ADVECTION terms are off — FCT carries advection explicitly via the full `w`), `Redi=.false.` → `isredi=0` kills
  every isoneutral `slope_tapered`/`Ki` term, `mix_scheme='PP'` → `use_kpp_nonlclflx=.false.` (no nonlocal flux),
  `smooth_bh_tra=.false.` (no biharmonic). So the FESOM3 TDMA = `Kv` vertical diffusion + the `bc_surface` row;
  the do_wimpl advection block is transcribed-but-unexercised, the rest deferred. **linfs makes the ALE
  reconstruct's `del_ttf += T·(hnode-hnode_new)` term vanish** (`hnode_new==hnode`) → `T* = T + del_ttf/hnode`
  (the zstar non-trivial case waits for M2.11). To exercise `diff_part_hor_redi` non-vacuously anyway, PRESCRIBE
  `Ki>0` (the L21 prescribe-the-unsourced-input pattern: `Ki` needs `mesh_resolution`, deferred to M4) and source
  `tr_xy` from the gated M1.1 `tracer_gradient_elements` (re-run per tracer; `do_oce_adv_tra` normally sets it).

- **Redi-off reads `slope_tapered`/`tr_z` BEFORE multiplying by `isredi=0` → `NaN·0=NaN` poisons the oracle.**
  The shipped namelist has `Fer_GM=.true.`/`Redi=.true.` so `slope_tapered`/`tr_z` are ALLOCATED but not yet
  FILLED at end-of-`ocean_setup`; FESOM2 computes `Fx=Kh·(Tx + SxTz·isredi)` where `SxTz=Σ(Tz·slope_tapered)/2` —
  if `slope_tapered`/`tr_z` are uninitialised garbage/NaN, `SxTz·0 = NaN` and `Fx=NaN`, while FESOM3 (which OMITS
  the slope terms for Redi-off) computes a finite `Fx=Kh·Tx` → the gate would MISMATCH (or NaN-propagate). Fix:
  the shim explicitly `slope_tapered=0`/`tr_z=0` (guarded by `allocated`) before the solve, so `SxTz=0` cleanly
  and both sides agree on `Fx=Kh·Tx` (`Tx+0.0·0.0 == Tx` to the bit). General rule: when transcribing a
  feature-OFF branch that DROPS a `×flag` term FESOM2 still evaluates, zero the dropped term's operands in the
  oracle so its `×0` is a clean finite 0, not `NaN·0`.

- **Oracle-side: distinguish module procedures from free subroutines — the symbol mangling differs.** The shim
  drives REAL FESOM2 routines via either a `use` of an auto-generated `*_interface` module (e.g.
  `diff_tracers_ale_interface`) OR an explicit `interface` block (for free subroutines with no `*_interface`
  module, like `impl_vert_visc_ale`/`mo_convect`). `tracer_gradient_elements` is a MODULE PROCEDURE of `o_tracers`
  → its symbol is mangled (`o_tracers_mp_tracer_gradient_elements_`), so an explicit external `interface` block
  emits a call to the UNMANGLED `tracer_gradient_elements_` → **`undefined symbol` at load time** (links fine,
  fails at runtime `symbol lookup error`). Fix: `use o_tracers, only: tracer_gradient_elements` (let the module
  provide the mangled name). Check `module … contains` membership before choosing interface-block vs `use`.

- **Two transcription gotchas the build caught.** (1) Fortran specification-expression order: a dummy whose
  array bound references ANOTHER dummy (`tr_xy(2, mesh%nl-1, …)`) must have that dummy (`mesh`) DECLARED EARLIER
  in the spec part — Intel `#6158`/`#6415`; reorder the type-decls (the arg LIST order is independent and stays
  matching the call). (2) `real_salt_flux` lives in `g_forcing_arrays`, NOT `o_ARRAYS` with its sibling surface
  fluxes `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt` — `use` the right module per symbol. Keep FESOM2's
  `zinv1=zinv2` backup bookkeeping in the TDMA verbatim (`1/dz` of the layer above) rather than recomputing
  `1/(Z_n(nz-1)-Z_n(nz))` inline — byte-identical, but faithful transcription removes reviewer doubt (L7).

## L24 — `step_oce` ASSEMBLY (M2.9b): the whole ocean step byte-matches the REAL `oce_timestep_ale` via its own built-in dumps; physical UV to avoid the `check_blowup` stop (PASSED first gate run)

M2.9b assembles every M2.1-M2.9a leaf kernel into one driver (`mod_step_oce::step_oce`, mirroring FESOM2
`oce_timestep_ale` + the `compute_vel_nodes` the main loop runs just before it) with **LIVE data flow** — each
kernel reads the PREVIOUS kernel's output, not a prescribed input — and gates the WHOLE step substep-by-substep.
- **The oracle dumps come for free from the REAL driver.** FESOM2's `oce_timestep_ale` already has built-in
  `dump_shim_record_node` calls at every substep, and FESOM3's `mod_dump` is byte-format-identical to that shim
  (M0.6). So the gate just prescribes the clean state + forcing arrays + reduced-M2 config and calls the REAL
  `compute_vel_nodes` + `oce_timestep_ale` ONCE — no new oracle dump code. FESOM3's `step_oce` emits the same
  records inline (same substep IDs/names), `dump_diff.py` compares (13 NODE fields × 5 probes = 65 records).
  This is the **mod_dump / `dump_diff.py`** gate vehicle (the run_oracle path), DISTINCT from the per-kernel
  FADVHDMP / `pressure_diff.py` gate — both coexist; the step gate is its own `run_step_gate.sh`.
- **Use the pi namelist (KPP/GM/Redi) + FORCE the reduced-M2 dispatch in the shim** rather than authoring a
  reduced-M2 namelist. The shipped pi namelist guarantees ALL the GM/Redi/KPP arrays (`sw_alpha`/`sigma_xy`/
  `neutral_slope`/…) are ALLOCATED at `ocean_setup`; the shim then just redirects the DISPATCH (`mix_scheme_nmb=2`,
  `Fer_GM=.false.`, `Redi=.false.`, `opt_visc=7`, …). The dead-in-M2 producers (`sw_alpha_beta`/`compute_sigma_xy`/
  `compute_neutral_slope`) still RUN inside `oce_timestep_ale` (FESOM3's `step_oce` OMITS them), so they dump the
  SW_AB substep (id=2) that FESOM3 never emits → pass `--ignore-substep=2` to `dump_diff.py` (a clean, documented
  deferral) rather than faking the records. The "reduced-M2 namelist running past forcing" the HANDOFF anticipated
  was NOT needed — the prescribe-and-stop shim drives `oce_timestep_ale` standalone (it reads forcing ARRAYS, not
  files; the 1-rank forcing-file hang L8 is bypassed exactly as for M2.1-M2.9a).
- **`max|Δ|=0` first run by L9 transitivity** — every operand was already byte-pinned per-kernel; M2.9b only adds
  the WIRING (uvnode← `compute_vel_nodes`, Av← PP into `impl_vert_visc_ale`, Kv← PP into the tracer TDMA, the full
  `solve_tracers_ale` loop, the linfs `update_thickness_ale` no-op). The data flow being correct is exactly what a
  whole-step substep match proves. `compute_vel_nodes` is validated transitively: PP's Kv (a dumped substep)
  consumes it, so Kv matching ⇒ uvnode matched.
- **`update_thickness_ale` is a no-op for linfs** (FESOM2 `oce_ale.F90:1226`): neither the zlevel nor the zstar
  thickness-redistribution branch is taken (hnode/helem/zbar_3d_n/Z_3d_n fixed), only `exchange_elem(helem)` —
  a 1-rank no-op. The `hnode` THICKNESS dump equals the init value on both sides.
- **`use_wsplit=.false.` is FORCED (M1.4 precedent).** `do_oce_adv_tra`'s `use_wsplit=.true.` path needs the FCT
  implicit vertical-advection correction `adv_tra_vert_impl` (a distinct unported kernel; guarded with an explicit
  `error stop`). With `.false.` the explicit/implicit split is trivial (`w_e=w`, `w_i=0`) AFTER `vert_vel_ale`, but
  `impl_vert_visc_ale` runs BEFORE `compute_Wvel_split` so it still consumes the PRESCRIBED `w_i` ≠ 0 (vertical
  momentum advection IS exercised). The w-split itself is gated at M2.7. Porting `adv_tra_vert_impl` +
  `use_wsplit=.true.` is a scoped follow-up.
- **Physical prescribed UV (0.50/0.40 m/s), NOT the viscosity gate's 2.0/1.5.** The strong-current stress test
  drives the single-step elevation `eta_n ~13 m`, past FESOM2's `±10 m` `check_blowup` guard → `oce_timestep_ale`
  STOPs ("eta_n become NaN or <-10,>10") AFTER all substep dumps fire (the gate still PASSES — dumps complete
  before the stop, verified by `hnode` matching). But an abnormal-termination + "NaN" log makes a permanent gate
  look fragile. Reducing UV to a physical 0.50/0.40 keeps `eta_n ∈ [-4.6, 1.9]` so the REAL step COMPLETES cleanly
  (the shim's own `stop` fires). The gate's correctness rests on the dump COMPARISON (all 13 substeps present +
  matched), not the exit code (`run_stepdump_pi.sh` tolerates the exit but checks dump existence + completeness),
  so a blowup would NOT silently pass — but a clean run is the trustworthy default. Branch coverage for the strong
  flow (the γ0/γ1/γ2 viscosity selection) stays at M2.4, not here.

## L25 — Forcing READ (M2.10a): the FIRST netCDF I/O — the 1-rank async-netcdf infinite-recursion fix, per-node partition-independence, the noleap `julday=365·yyyy` (710820 not 2.4e6) two-stage time-interp cancellation (PASSED first gate run)

The first kernel that reads real files (CORE2 NCAR stubs `test/input/global/{u_10,v_10,q_10,ncar_rad,t_10,ncar_precip}.1948.nc`,
192×94, 5 records, ~360 KB each). All 8 atmospheric fields (`u_wind`/`v_wind`/`Tair`/`shum`/`shortwave`/`longwave`/
`prec_rain`/`prec_snow`) `max|Δ|=0` vs the REAL FESOM2 `sbc_do` on pi 1-rank. Built `src/io/mod_io_netcdf.F90`
(thin `use netcdf` wrapper) + `src/forcing/mod_forcing_read.F90` (the read+interp) + `vector_g2r` into
`mod_mesh_rotate` + driver `src/drivers/fesom_forcingdump.F90`; gate `tools/run_forcing_gate.sh`.

- **The L8 "1-rank forcing hang" is a REAL infinite recursion, now FIXED in the oracle.** Empirically: a 1-rank
  full run hangs in `g_sbf::getcoeffld → io_netcdf_workaround::next_io_rank → mpi_topology`. Root cause:
  `next_io_rank_helper` (its own TODO admits it) infinite-recurses when the only rank IS `SEQUENTIAL_IO_RANK`(0) —
  the async-NetCDF IO-rank selector. Patched `port2/fesom2/src/io_netcdf_workaround_module.F90`: at `partit%npes==1`
  short-circuit to sequential I/O on rank 0 (`async_netcdf_allowed=.false.`) — **VALUE-IDENTICAL** (only changes
  which rank reads+bcasts, not the bytes). After the patch the full model runs PAST forcing and through a clean
  first step (`FDBG step=1 uv=0.13 ssh=0.10 Smax=37.4`, no NaN); it then crashes at OUTPUT I/O
  (`io_meandata::init_nod2d_lists`, a SEPARATE 1-rank gather bug) — irrelevant, the prescribe-and-stop shim stops
  before output. This is an UNCOMMITTED FESOM2 working-tree change (like the dump shims).
- **Forcing is PARTITION-INDEPENDENT per node** (each node's value is a pure function of its own geo-coords + the
  global file data + the model time — NO cross-node accumulation), unlike `area`/FCT scatter. So a multi-rank
  oracle would byte-match too; the np=1 fix simply keeps the 1-rank anchor (consistent with every prior gate).
- **The oracle shim** (`port2/fesom2/src/fesom_forcing_dump.F90`) is wired right AFTER `forcing_setup` in
  `fesom_module.F90:315` (NOT end-of-`ocean_setup` — `sbc_ini` only runs inside `forcing_setup`, which is AFTER
  `ocean_setup` and gated on `use_ice=.true.`, true on pi). It pins the model time, drives the REAL `sbc_do`, maps
  `atmdata→arrays` exactly as `gen_forcing_couple::update_atm_forcing:681-694`, dumps (FADVHDMP), stops.
- **`julday` for `noleap`/`none`/`365_days` = `365·yyyy`** (the else branch, `gen_surface_forcing.F90:1887`), so the
  pi CORE2 time axis is `365·1948=710820`-scale — NOT the ~2.43e6 NR Julian-Day-Number of the JRA gregorian case
  (that's where the HANDOFF's "2.4e6 cancellation" figure comes from). Same cancellation PRINCIPLE, smaller magnitude.
- **The time interp is a TWO-STAGE form that must NOT be collapsed**: `coef_a=(data2-data1)/Δt`;
  `coef_b=data1-coef_a·nc_time(t_indx)`; `atmdata=rdate·coef_a+coef_b`. Both `rdate·coef_a` (~710820·a) and `coef_b`
  (~−710820·a) are large and nearly cancel to the O(1) physical value. Refactoring to the algebraically-equal
  `data1+coef_a·(rdate−nc_time)` rounds the large `coef_b` differently → drift. (FESOM2 getcoeffld:1015-1016 +
  data_timeinterp:1041.)
- **Two rdates**: the coefficients are built at the COLD-START rdate (`nc_sbc_ini:643` — NO half-step, clock-init
  day 1 / sec 0) and the data is evaluated at the per-step rdate (`sbc_do:1527` — WITH the `−dt/86400/2` half-step,
  pinned day 1 / sec 43200 / dt 2400). The gate driver replicates BOTH. (For pi both land in the same window
  → `t_indx=1`, so a single-rdate read would also match, but the two-phase form is faithful and general.)
- **Raw slices stay `real(4)`** (the on-disk type) and promote to WP inside the bilinear weight expression — exactly
  FESOM2 (`real(4) sbcdata · real(WP) wgt`); `real4→real8` is exact so it's byte-identical either way, but keeping
  `real(4)` reads the on-disk bytes verbatim. The periodic-lon halo (`nlon+2`, mirror cols 1↔nlon-1 / nlon↔2, then
  `ic_cyclic` ±360) makes lon monotonic for `binarysearch`; lat-flip only if descending (CORE2 ascends → no flip).
- **netCDF link via `nf-config`** in CMakeLists.txt (the SAME spack `netcdf-fortran-4.5.3` the oracle links → byte
  reads) + **`-Wl,-rpath`** to both netcdf-fortran/-c lib dirs so binaries run without `LD_LIBRARY_PATH` (the spack
  modules set compile paths but not runtime). `use netcdf` (the F90 `.mod`) is compatible across intel 2021.5.0
  (netcdf build) ↔ 2022.0.1 (FESOM3) — no need for the F77 `netcdf.inc`.
- **The dump drivers do NOT call `set_partition`**, so `partit%myDim_nod2D` stays 0 — loop forcing over `mesh%nod2D`
  (== `myDim+eDim` at 1-rank), the same convention every FESOM3 kernel uses (`oce_pressure_bv.F90:72`
  `do node=1,mesh%nod2D`). The driver sets `frc%nnod=mesh%nod2D` and the read routines loop `1,frc%nnod`.
- **SCOPE**: M2.10a is the READ only. `heat_flux`/`water_flux` come from the air-sea **obudget** in
  `ice_thermo_oce.F90` (+ `oce_fluxes`), which is **M3 (ice/thermo)** per the plan ("M3 → thermo → oce_fluxes") —
  even for open water. So M2.10b (bulk transfer coeffs `Cd/Ch/Ce` + wind `stress_surf`) + M2.10c (SW penetration
  `sw_3d`) finish M2.10's producible fields; `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt` stay prescribed
  until M3. `gen_bulk_formulae.F90` computes ONLY the transfer coefficients (NOT stress/heat/evap).
- **M2.10b (bulk + wind stress, 17 fields, PASSED first gate run).** `ncar_ocean_fluxes_mode` (the LIVE NCAR routine:
  Large&Yeager 2004 + **Large-2009 drag** [the `u10**6` term + the 33 m/s → `2.34e-3` cap, NOT the commented-out
  L-Y2004 6a form], `n_itts=5`, a 3-height Monin-Obukhov stability iteration) → `Cd`/`Ch`/`Ce` — gated against the REAL
  routine. Then `stress_atmoce=Cd·(ρ_air·|Δu|)·Δu` (`Δu=u_wind−(1−Swind)·u_w`, `Swind=0`) and the node→elem
  `stress_surf(elem)=sum(stress_node_surf(elnodes))/3` (`a_ice=0` ⇒ `stress_node_surf=stress_atmoce`). **Byte traps —
  transcribe the un-suffixed default-real literals VERBATIM** (the L16 family): `inc_ratio=1.0e-4`, `inv_rhoair=1./1.3`,
  `tmelt=273.15`, `rhoair=1.3` (MOD_ICE type-defaults / the gen_bulk local) — a `_WP` suffix would round differently;
  also `(ustar*ustar)` not `ustar**2`, `atan(1.0_WP)` kept a runtime call, `test=abs(cd−cd_prev)/(cd+1.0e-8_WP)`.
  **SST + surface ocean velocity are PRESCRIBED** (a fresh `type(t_ice)` dummy supplies the thermo TYPE-DEFAULTS —
  `ice%thermo` is a non-allocatable component so `inv_rhoair`/`tmelt`/`rhoair` auto-initialize; only `srfoce_temp/u/v`
  need allocation). The wind stress + node→elem are inlined in BOTH sides (trivial formulas, L9-transitive — `Cd` is the
  gated substance). Non-vacuous: `Cd∈[5.5e-5,1.7e-2]`, the wide Tair(−44..31)−SST(−1..19) range fires both
  stability branches. `elem2D_nodes(1:3,elem)` slice (FESOM3 MAX_NV=4; FESOM2's is (3,·)) — the L15 trap. The bulk
  consumes the LIVE M2.10a `u_wind`/`v_wind`/`Tair`/`shum` (a real read→bulk assembly).
- **M2.10c (shortwave penetration, 20 fields total, PASSED first gate run).** `cal_shortwave_rad`
  (`oce_shortwave_pene.F90`, Morel&Antoine 1994 / Sweeney 2005): `swsurf=(1−albw)·shortwave·0.54`; `heat_flux+=swsurf`
  (the visible band is REMOVED from the +upward non-solar `heat_flux` and redeposited as `sw_3d`); chl floor `0.02`;
  the v1/v2/sc1/sc2 polynomials in `c=log10(chl)`; the two-exponential `sw_3d(k)=swsurf·(v1·exp(zbar/sc1)+v2·exp(zbar/
  sc2))` over `zbar_3d_n` (the per-node ALE depth, pressure-gate-proven) with the `aux<1e-5`/`k==nzmax` cutoff;
  `swsurf/=vcpw` (W/m²→K·m/s, `vcpw=4.2e6` — exactly representable so `_WP` is moot). Consumes the LIVE M2.10a
  `shortwave`; `chl`/`heat_flux`/`a_ice=0` prescribed. **`albw=0.066` forced to `0.066_WP` on BOTH sides** (the shim
  sets `ice%thermo%albw=0.066_WP`) so the un-suffixed default-real ambiguity can't bite. The shim's dummy ice needs
  `ice%data(1)%values` (= `a_ice`) allocated (the routine pointer-assigns it); the REAL `cal_shortwave_rad` is driven
  via an explicit interface. Non-vacuous: the chl floor fires on 476 polar nodes; `sw_3d` decays to 0 over 48 levels.
- **Build/CMake recap (M2.10):** netCDF added to FESOM3 via `nf-config` (include+`--flibs`) + `-Wl,-rpath` to the
  netcdf-fortran/-c lib dirs (self-contained binaries; no `LD_LIBRARY_PATH`). `use netcdf` (the F90 `.mod`) works
  across the intel 2021.5.0 (netcdf) ↔ 2022.0.1 (FESOM3) minor-version gap. `src/io/` + `src/forcing/` were already in
  `FESOM3_LIB_DIRS` (auto-globbed once created). One forcing gate (`run_forcing_gate.sh`) now covers all 20 fields
  (read→bulk→SW pene, a real producer chain on the LIVE read); `pressure_diff.py` is the generic comparator (its
  "PRESSURE/EOS…" trailer is cosmetic).
