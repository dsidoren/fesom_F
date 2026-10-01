# rvo_upwind: upwind-blended relative-vorticity reconstruction

## Overview

The vector-invariant momentum advection (`momadv_opt==1`, `src/oce/oce_dyn_vinv.F90`)
reconstructs relative vorticity at an element face as the plain average of its three
vertex values. This plan adds an **upwind-blended** reconstruction controlled by a single
coefficient:

```
zb_blended = zb + rvo_upwind * (zb_upwind - zb)
```

- `rvo_upwind = 0`: the current reconstruction, **bit-for-bit** (structurally — no new
  code executes on that path).
- `rvo_upwind = 1`: fully upwind — the face value comes entirely from the one or two
  edge-sharing neighbour elements that lie upstream, weighted by inflow strength.

Upwinding the advected vorticity adds the usual upwind damping to the vorticity /
enstrophy dynamics. It is **energy-neutral by construction**: `u·[(f+ζ)×u] ≡ 0` for *any*
ζ (the cross product is perpendicular to `u`), so the blend cannot pump or drain kinetic
energy through the Coriolis-vorticity term.

## Context (from discovery)

- Target: `src/oce/oce_dyn_vinv.F90` Block B2, where `zb = (1/3)·Σ vort(nz,elnodes)`.
- `mesh%elem_neighbors(MAX_NV, elem2D)` is **declared in `mod_mesh` but never allocated
  or filled** anywhere in the tree. `elem_edges` likewise — and stays unbuilt (see
  Decision 2).
- `edges` / `edge_tri` exist at multirank over `nEdgeL = myDim+eDim` edges with localized
  element ids (`mod_mesh_read.F90:381`).
- **Critical discovered constraint:** `edge_dxdy` / `edge_cross_dxdy` / `edge_len` exist
  for **owned edges only** (`nEdgeO`, `mod_mesh_areas.F90:320-321`), and an owned element
  can have one halo edge. Any design reading per-edge geometry therefore breaks at
  multirank. This killed the first design draft and forced Decision 3.
- `exchange_elem` exists (`mod_halo.F90:44`).
- No FESOM2 oracle for anything on this path (v2.7.3 `momadv_opt==1` aborts); the net is
  `test/test_vinv.F90` (V1–V7, all exact) plus the conservation gate (15 configs).

## Development Approach

- **testing approach**: **TDD (tests first)** — established preference for no-oracle
  schemes; the new assertions are written before the blend exists and are the
  specification.
- complete each task fully before moving to the next
- make small, focused changes
- **CRITICAL: every task MUST include new/updated tests** for code changes in that task
- **CRITICAL: all tests must pass before starting the next task** — no exceptions
- **CRITICAL: update this plan file when scope changes during implementation**
- run tests after each change
- backward compatibility is absolute: `momadv_opt==2` and `rvo_upwind=0` stay
  bit-identical

## Testing Strategy

- **unit tests**: extend `test/test_vinv.F90` (np 1 and 2). V1–V7 must stay green
  **untouched** — V2/V6 in particular remain valid for *any* `rvo_upwind` because a
  constant-ζ field has `zb_upwind == zb` (blend exactly inert), a property exploited
  deliberately by V9.
- **integration**: `tools/run_conserve_pi.sh` gains `momadv-vinv-upw`
  (`FESOM3_MOMADV_OPT=1 FESOM3_RVO_UPWIND=0.7`) at np=1 **and np=8** — `exchange_elem`
  is the new multirank risk, and np=8 is the rank count that exposes local/global array
  bugs (the gate's own header documents why).
- no e2e tests in this project.
- Commands: `source env/levante.dkrz.de/shell.intel`; `cd build_intel_dp && make -j16 &&
  ctest`; `bash tools/run_conserve_pi.sh`.

## Progress Tracking

- mark completed items with `[x]` immediately when done
- add newly discovered tasks with ➕ prefix
- document issues/blockers with ⚠️ prefix
- update plan if implementation deviates from original scope

## Solution Overview

### Decisions (validated in the brainstorm)

**Decision 1 — upwind criterion: normal flux, own-element velocity.** The prescription's
literal `edge_dxdy·u` is the *along-edge* (tangential) component, whose sign
distinguishes edge ends, not sides. The stated goal — "the one or two neighboring
elements that lie upwind" — is measured by the **inward-normal** component: `d > 0` ⟺
flow enters through that edge ⟺ that neighbour is upwind; a uniform flow gives exactly
1–2 inflow edges, and `d_up = d + |d|` with normalization works as prescribed. The
velocity in `d` is the element's own `UV(:,nz,elem)`.

**Decision 2 — cache `elem_neighbors` only; nothing else.** Per the user: the neighbour
across an edge is derivable from `edge_tri` ("the other entry"). What is *not* derivable
on the fly is the inverse map (element → its edges), so the neighbour ids are cached
**once** at init — 3 ints per owned element — by exactly that `edge_tri` check, swept
over the local edges. No `elem_edges`, no per-step `edge_tri` reads.

**Decision 3 — normals rebuilt locally from the element's own vertices.** Edge `k`
connects `elnodes(k)` and `elnodes(k+1)`, both local; the edge vector comes from
`coord_nod2D` with `trim_cyclic` and the element's own `elem_cos`. Of the two
perpendiculars, the **inward** one points toward the third vertex `elnodes(k+2)` — a
local, convention-free orientation (no `edge_tri` left/right branch; the class of sign
error that bit the circulation integral is structurally excluded). This also removes the
halo-edge geometry gap: no per-edge array is read at all.

**Decision 4 — dry/boundary neighbours are *excluded*, not zero-valued.** `zb_upwind` is
a pointwise **value** estimate, not a flux: by the project's own rule,
zero-with-full-weight is for flux-forming quantities; estimators drop dry contributors.
All excluded, or no inflow at all → `sumw = 0` → blend inert.

**Decision 5 — `omega_e` is a local allocatable, not state.** `(nl-1, nElemF)` inside
`momentum_adv_vinv`, allocated **only when `rvo_upwind > 0`**, filled for owned elements
as the same vertex average, then **one `exchange_elem`** so halo neighbour ids read valid
values. This sidesteps the owned-only `elem2D_nodes` problem (a halo neighbour's centre
value arrives by exchange; its vertex list is never needed). No `t_dyn_work` field, no
driver churn, no restart impact.

## Technical Details

### Adjacency build (once, at geometry setup)

`build_elem_adjacency(mesh, partit)` in `src/mesh/mod_mesh_areas.F90`, called from
`compute_geometry`:

```
allocate(mesh%elem_neighbors(3, nElemO));  elem_neighbors = -1     ! -1 = unresolved
do ed = 1, nEdgeL                       ! ALL local edges (owned + halo)
    (a, b) = mesh%edges(:, ed);  (e1, e2) = mesh%edge_tri(:, ed)
    for each of e1, e2 that is an OWNED element (1 <= e <= nElemO):
        k = the slot whose vertex pair {elnodes(k), elnodes(k+1)} == {a, b}
        elem_neighbors(k, e) = the OTHER element id   ! halo id ok; 0 if boundary (<=0)
end do
VERIFY: no owned element retains an unresolved (-1) slot -> error stop with element id.
```

The verification turns the partition assumption ("an owned element sees all three of its
edges in the local edge list") into a checked invariant, house style
(`assert_bottom_invariant`).

### The blend (Block B2, inside the level loop)

Per element, the three inward normals `nx(k), ny(k)` are built once **outside** the `nz`
loop (Decision 3). Then:

```fortran
zb = onethird*sum(vort(nz,elnodes))                       ! unchanged default
if (rvo_upwind > 0.0_WP) then
    sumw = 0 ;  zup = 0
    do k = 1, 3
        nb = mesh%elem_neighbors(k, elem)
        if (nb <= 0) cycle                                ! boundary
        if (nz < mesh%ulevels(nb) .or. nz > mesh%nlevels(nb)-1) cycle  ! dry here
        d = UV(1,nz,elem)*nx(k) + UV(2,nz,elem)*ny(k)     ! inflow through edge k
        w = d + abs(d)                                    ! 2d if inflow, 0 if outflow
        sumw = sumw + w
        zup  = zup  + w*omega_e(nz, nb)
    end do
    if (sumw > 0.0_WP) zb = zb + rvo_upwind*(zup/sumw - zb)
end if
```

Properties to state in the code comments: convex combination → max principle (no new
vorticity extrema); energy-neutral (`u·[(f+ζ)×u] ≡ 0` for any ζ); acts on
vorticity/enstrophy dynamics only.

### Parameter

`mod_param_phys`: `real(kind=WP) :: rvo_upwind = 0.0_WP`, documented range [0,1].
**Value-based** env `FESOM3_RVO_UPWIND` (like `FESOM3_MOMADV_OPT`, unlike the
presence-based switches); out-of-range → `error stop`, never a silent fallback.

## What Goes Where

- **Implementation Steps**: all code, tests, gate config, docs — this repo.
- **Post-Completion**: choosing a production `rvo_upwind` value and running the
  comparison — the user's call.

## Implementation Steps

### Task 1: Build and verify `elem_neighbors`

**Files:**
- Modify: `src/mesh/mod_mesh_areas.F90`
- Modify: `test/test_vinv.F90`

- [ ] add `build_elem_adjacency` (sweep over `1..nEdgeL` local edges; slot `k` matched by
      vertex pair; neighbour = other `edge_tri` entry, halo ids kept, boundary → 0)
- [ ] call it from `compute_geometry`; `-1`-init + post-sweep verification with
      `error stop` naming the first unresolved element
- [ ] check `mod_io_meshdiag`'s `elem_neighbors` references behave now the array is
      genuinely allocated
- [ ] write tests: every owned pi element has 3 resolved slots; each neighbour pair is
      mutual where both are owned (`elem ∈ neighbors(neighbors(elem))`); boundary count
      equals the count of `edge_tri(2,·) <= 0` edges
- [ ] run tests — 32/32 + new assertions green at np 1 and 2 before task 2

### Task 2: Write the upwind assertions V8–V11 (TDD — fail until Task 3)

**Files:**
- Modify: `test/test_vinv.F90`

- [ ] V8 bit-identity: `rvo_upwind=0` call of `momentum_adv_vinv` gives `UV_rhsAB`
      exactly `==` a reference call
- [ ] V9 constant-ζ inertness: rerun the V2 linear-shear exactness and the V6
      uniform-`u`/nonzero-`w` zero-tendency **with `rvo_upwind = 1.0`** — identical
      bounds (constant ζ ⟹ `zb_upwind == zb` ⟹ blend exactly inert)
- [ ] V10 orientation proof: uniform flow, `rvo_upwind=1`; the test recomputes normals +
      weights from public mesh data and asserts, for every owned element and edge with
      `w > 0`, that the neighbour centroid is upstream:
      `(centroid(nb) − centroid(elem))·u < 0` (owned neighbours only — halo centroids
      would need the owned-only `elem2D_nodes`)
- [ ] V11 end-to-end weighting: `u ∝ y²` (ζ ∝ y, **`y` from `coord_nod2D` in the rotated
      frame** — the V2 lesson); recover
      `δzb = (Δrhs_x·V − Δrhs_y·U)/((U²+V²)·elem_area)` between `rvo=1` and `rvo=0`
      calls (guard `U²+V² > tiny`); exactly one wet owned inflow neighbour ⟹
      `δzb == omega_e(nb) − zb₀` to round-off; two ⟹ max principle (blend within
      `[min,max]` of candidates)
- [ ] run tests — V8/V9 pass trivially against the stub-free current code only where
      inert; V10/V11 are EXPECTED TO FAIL (no blend exists). Record the failures as the
      specification; mark `[x] (fails until Task 3)`

### Task 3: Implement the blend in Block B2

**Files:**
- Modify: `src/oce/oce_dyn_vinv.F90`
- Modify: `src/params/mod_param_phys.F90`

- [ ] add `rvo_upwind = 0.0_WP` to `mod_param_phys` with the range/meaning comment
- [ ] `momentum_adv_vinv`: when `rvo_upwind > 0`, allocate local `omega_e(nl-1, nElemF)`,
      fill owned elements with the vertex average, `exchange_elem` under `is_multirank`
- [ ] per-element inward normals from own vertices (`coord_nod2D`, `trim_cyclic`,
      `elem_cos(elem)`, orientation toward `elnodes(k+2)`), outside the `nz` loop
- [ ] the blend exactly as in Technical Details, with the exclusion fallbacks and the
      max-principle / energy-neutrality comments
- [ ] verify V8–V11 now PASS at np 1 and 2; V1–V7 untouched and green
- [ ] run tests — full `ctest` green before task 4

### Task 4: Plumbing and the gate

**Files:**
- Modify: `src/drivers/fesom_conserve.F90`
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90`
- Modify: `tools/run_conserve_pi.sh`
- Modify: `work/job_levante` (commented-out example only; default OFF)

- [ ] `FESOM3_RVO_UPWIND` value-based read in `fesom_conserve` (beside
      `FESOM3_MOMADV_OPT`), validated to [0,1] with `error stop`
- [ ] same read in `fesom_lifecycle_native_mr`, rank-0 banner when `> 0` placed **at the
      assignment** (the momadv-banner lesson: never inside another option's branch)
- [ ] gate: `momadv-vinv-upw FESOM3_MOMADV_OPT=1 FESOM3_RVO_UPWIND=0.7` in the heredoc
      (np=1) and an explicit np=8 run of the same config
- [ ] `job_levante`: commented-out `# export FESOM3_RVO_UPWIND=...` with a one-line note
- [ ] write tests: the gate run IS the test — all previous configs plus the two new ones
      green
- [ ] run `bash tools/run_conserve_pi.sh` — fully green before task 5

### Task 5: Verify acceptance criteria
- [ ] `rvo_upwind=0` structurally inert (no allocation, no exchange, zb path identical)
- [ ] `rvo_upwind=1` fully upwind; intermediate values blend linearly (spot-check 0.5 in
      V11's recovered `δzb`)
- [ ] dry-neighbour, boundary and `u=0` fallbacks covered (V11 mesh has all three)
- [ ] run full test suite `ctest` and the full gate
- [ ] `momadv_opt==2` path untouched (unreachable) — suite green is the proof

### Task 6: [Final] Update documentation
- [ ] `docs/HANDOFF.md`: `rvo_upwind` / `FESOM3_RVO_UPWIND` under the momadv section
- [ ] `docs/LESSONS.md`: one entry — per-edge geometry is owned-edges-only, so
      element-local reconstruction (normals from own vertices, neighbour values by
      `exchange_elem`) is the MR-safe pattern; orientation chosen toward the third vertex
      needs no convention
- [ ] move this plan to `docs/plans/completed/`

## Post-Completion

*No checkboxes — informational.*

- **Choosing the production value**: `rvo_upwind` is a dissipation knob on the vorticity
  dynamics. A short pi/core2 pair (0 vs ~0.5) comparing eddy kinetic energy and
  enstrophy spectra is the natural first experiment; enabling it in `job_levante` is the
  user's call.
- Possible later refinement (explicitly out of scope now): edge-mean velocity in `d`
  (symmetric upwind criterion), and applying the same blend to the `(f+ζ)` average used
  in the Coriolis completion if PV-consistent upwinding is ever wanted.
