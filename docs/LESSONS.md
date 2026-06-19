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
