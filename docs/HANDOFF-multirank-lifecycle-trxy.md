# Handoff — the multi-rank free-running lifecycle `tr_xy` halo-exchange bug (✅ SOLVED 2026-06-22)

**Date:** 2026-06-22. **Status: SOLVED — the multi-rank free-running CORE2 lifecycle is now BYTE-EXACT vs FESOM2**
(`run_lifecycle_gate_multirank.sh 2 3` → 195 records, worst `|Δ|=0.000e+00`), no 1-rank/pi regression (ctest 13/13,
`run_step_gate_multirank.sh 2` and `run_step_gate.sh` → 65/65 `=0`).

## RESOLUTION (read this first)

It was **NOT a FESOM3 bug** — the exchange logic is correct. It was an **OpenMPI 4.1.2 `vader` (shared-memory) BTL
KNEM single-copy bug** that silently corrupts large messages (>~131 KB) on levante. FESOM3's `mod_halo::core_blk_r`
does a contiguous nonblocking `Isend/Irecv`; for the CORE2-scale `tr_xy` element-block exchange (94×636 = 478 KB —
the FIRST live F3 message above the threshold) KNEM delivered the first 134656 bytes correctly then filled the tail
with garbage. A cross-rank buffer trace pinned it precisely: packed `sbuf` (correct) → in-transit `rbuf` (corrupt from
byte 134656) → unpack `arr` (faithful) — the corruption is purely in the MPI transport, between pack and unpack.

- **Why it hid until the CORE2 lifecycle:** the bug needs a >~131 KB message. Every pi-mesh gate's halos stay well
  under it; FESOM2's `MPI_TYPE_INDEXED` exchange uses a different `vader` path that dodges KNEM; the 1-rank lifecycle
  has no exchange. CORE2 × free-running × multi-step was the first combination to send a large live F3 message.
- **Fix (committed in `env.sh`, levante branch):** `export OMPI_MCA_btl_vader_single_copy_mechanism=none` — forces the
  byte-faithful copy-in/copy-out path. Proof it's the transport, not us: the SAME binary/buffers give `total-diff=0`
  under `--mca btl self,tcp` AND under `single_copy_mechanism none`; only default `vader`+KNEM corrupts. (CMA is
  unavailable here — `kernel.yama.ptrace_scope=3` blocks the ptrace it needs — so `vader` falls back to the buggy KNEM.)
- **Both prior "leading hypotheses" were WRONG:** `myList_elem2D` (owned AND halo) is read VERBATIM from the same
  `dist_N/my_list*.out` as FESOM2 → identical ordering, so there was no eXDim-ordering ghost; and the manual pack vs
  persistent-`MPI_TYPE_INDEXED` difference mattered only in that the typed path happens to dodge KNEM. The full lesson
  + meta-lessons are in `docs/LESSONS.md` L35.
- **The UV over-exchange fix (§0 item 1) is real and KEEP it** (eDim `exchange_elem`, matching FESOM2).
- **Remaining (cleanup, before final commit):** strip the diagnostic instrumentation listed in §6 (F3 + F2 oracle);
  none of it is needed now. `tools/trxy_diff.py` / `trxy_analyze.py` / `trxy_repro.sh` are throwaway diagnostics.

---

<details><summary>Original investigation notes (pre-resolution; kept for the record)</summary>

**Original status (now obsolete):** the multi-rank *free-running* lifecycle is STABLE and matches FESOM2 to
`~1e-6` globally, but is **NOT byte-exact**. One residual bug remains, precisely localized but mechanism-elusive.
This doc is self-contained: read it cold next session and resume.

---

## 0. TL;DR

- We built the **production-validation gate**: a multi-rank **free-running** reduced-M2 lifecycle on CORE2 dist_2
  (dt=1800), the FIRST test of *multi-rank × free-running × many-steps* (every prior MR gate was single-step
  prescribe-and-stop, so latent bugs were masked).
- It immediately found **two** latent multi-rank bugs that single-step prescription hid:
  1. **UV halo blow-up — FOUND + FIXED.** `update_vel` exchanged UV over `com_elem2D_full` (eDim+eXDim);
     FESOM2 uses `com_elem2D` (eDim). Wrong UV halo → viscosity → blow-up (eta→2708 by step 3). Fix: eDim
     exchange. Model now stable, matches FESOM2's physical values to all printed digits.
  2. **`tr_xy` halo exchange byte-drift — OPEN (this doc).** `tr_xy` (the element tracer-gradient) exchanged over
     `com_elem2D_full` produces a **wrong halo at step 2+** (owned byte-exact). Feeds the MUSCL → `del_ttf` → T/S
     drift (~0.03 at boundary nodes) → amplifies through the dynamics to gross step-3 `ssh_rhs` (~1e6) at the
     central-Arctic nodes (the free-surface CG damps it, so the model stays globally stable at ~1e-6).

The bug is **provably impossible** by every verified fact (symmetric com, byte-exact owned source, correct gid
alignment, clean exchange of a gid pattern) yet reproducibly corrupts the *gradient* data at step 2. Points at a
subtle compiler/MPI-async or eXDim-second-layer subtlety we could not isolate by inspection.

---

## 1. Where we are (context)

- **Tag `m2-mvp`** (commit `908b25d`) = whole-model multi-rank byte-match vs FESOM2, **single-step
  prescribe-and-stop**, pi dist_2/dist_8, `max|Δ|=0`. That milestone is solid and committed.
- The **production validation** (optional secondary validation per the plan) is the multi-rank analog of the
  1-rank CORE2 lifecycle (`fesom_lifecycle` / `run_lifecycle_gate_core2.sh`, which is byte-exact at 1-rank, L29).
- ⚠️ The `com_elem2D_full` exchange bug is **latent in m2-mvp too** — m2-mvp's single-step gates prescribe the
  halo, so the buggy exchange is never exercised on live multi-step data. The free-running lifecycle is what
  exposed it. So this is a real correctness gap in the multi-rank port, not just a lifecycle-harness issue.

## 2. What was BUILT this session (all uncommitted, in the working tree)

| File | What |
|---|---|
| `src/oce/oce_initial_state.F90` | `do_ic3d`/`extrap_nod3D`/`getcoeffld`/`nc_ic3d_ini` lifted to MR (optional `partit`, owned/owned+halo bounds, `extrap` initial/per-iter/final `exchange_nod` + `allreduce_max(glob_max)`). assumed-shape `arr(:,:)`/`values(:,:)` so whole-array ops span the LOCAL size. |
| `src/oce/oce_pressure_bv.F90` | `insitu2pot` optional `partit` (owned+halo loop). |
| `src/infra/mod_halo.F90` | **NEW** `allreduce_max` (scalar, MPI_MAX). **`core_blk_r` rewritten** to flat rank-1 buffers + explicit indexing (was rank-3 buffers — equivalent, did NOT fix the bug). **`BLK_TAG=2`** distinct tag for `core_blk_r` (did NOT fix the bug). |
| `src/drivers/fesom_lifecycle_mr.F90` | **NEW** MR free-running lifecycle driver (= `fesom_stepfull_mr` scaffolding [remap + local sizes + gid dumps] + `fesom_lifecycle`'s cold-start + `do_ic3d(…,partit)` IC + N-step `step_oce(…,partit)` loop). |
| `src/oce/oce_ale.F90` | **THE UV FIX** — `update_vel` `exchange_elem_full(UV)` → `exchange_elem(UV)` (eDim, matching FESOM2). |
| `tools/run_lifecycle_gate_multirank.sh` | **NEW** MR lifecycle gate (oracle at NP + `fesom_lifecycle_mr` at NP, per-rank gid-keyed `dump_diff --glob`). |
| `tools/run_lifecycle_core2.sh` | added 4th arg `[nranks]` (default 1) so the oracle runs at NP>1 on CORE2 dist_NP. |

**Real fixes/features to KEEP:** the `do_ic3d` MR lift, the UV eDim fix, `fesom_lifecycle_mr`, the gates,
`allreduce_max`, the flat `core_blk_r` (cleaner), `BLK_TAG`. **Diagnostic scaffolding to strip** once the bug is
solved (see §6).

## 3. The OPEN bug — `tr_xy` (`com_elem2D_full`) halo exchange corrupts at step 2

### 3.1 What `tr_xy` is and how it flows
`tr_xy(2, nl-1, nElem)` = the element-wise horizontal gradient of the tracer (`∂T/∂x`, `∂T/∂y`). In
`init_tracers_AB` (FESOM2 `oce_tracer_mod.F90`; F3 `src/oce/oce_tracer_mod.F90`):
1. `tracer_gradient_elements` computes `tr_xy` over **OWNED** elements (`do elem=1, myDim_elem2D`).
2. **exchange the halo** — FESOM2 `exchange_elem_begin/end(tr_xy)` (before `fill_up_dn_grad`); F3
   `exchange_elem_full(tr_xy, partit)`.
3. `fill_up_dn_grad` reads `tr_xy` at the **up/downwind triangles** of each edge → these reach the **eXDim**
   (second halo layer), so the FULL halo (eDim+eXDim) must be valid → `edge_up_dn_grad`.
4. The MUSCL flux (`adv_tra_hor_mfct`) reads `edge_up_dn_grad` → `del_ttf` → tracer update.

So a wrong `tr_xy` halo → wrong `edge_up_dn_grad` at boundary edges → wrong MUSCL flux → T/S drift at
partition-boundary nodes.

### 3.2 The DEFINITIVE evidence (the `tr_xy` element dump)
We dumped `tr_xy(1,1,e)` for **every local element** (gid-keyed) from **both** codes, each `init_tracers_AB`
call (2-step run → records: step1-tr1, step1-tr2, step2-tr1, step2-tr2), and compared per-rank:

```
rank0: myDim=121861  eDim=285  eXDim=351  (nElemF=122497, halo=636)
  step1-tr1: OWNED diff=0  HALO diff=0          <- byte-exact
  step2-tr1: OWNED diff=0  HALO diff=456/636    <- HALO drifts (eDim 206/285 + eXDim 250/351)
rank1: myDim=123360 eDim=274  eXDim=331  (halo=605)
  step2-tr1: OWNED diff=0  HALO diff=426/605
```
Actual values (rank0, step2, a halo element owned by rank1):
```
  gid=76020  F2=+2.890428e-08   F3=+3.497329e-11   (same gid; element owned by r1)
```
- **F2's halo holds the real gradient (~1e-6..1e-8); F3's holds wrong values (~1e-11, a few ~1e-3).**
- The element's OWNER (rank1) has byte-exact owned `tr_xy` (=2.89e-8 on both F3 and F2). r1 SENDS it; the com is
  symmetric and aligned; yet r0's halo RECEIVES 3.5e-11.
- **F2's halo CHANGES every step (636/636 differ step1→step2)** → FESOM2 fills the full halo (eDim+eXDim) each
  step (so FESOM2's exchange is NOT eDim-only here; it covers eXDim).

### 3.3 The downstream chain (why the model is still ~1e-6 stable)
First gate divergence: **step 2, substep 15 (tracers)**, T/S ~0.03 at boundary nodes. By step 3 it amplifies
through density→pgf→vel_rhs→UV→ssh_rhs to gross `ssh_rhs ~1e6` at the naturally-huge-`ssh_rhs` **central-Arctic**
nodes (gid 100137: F2=5.8e4 vs F3=−1.2e6) — exactly the C-port's marginal-2Δx region. The free-surface CG damps
it (rtol scales with |b|), so `eta_n`/`d_eta` stay globally physical and the worst gated `|Δ|` is only ~1.4e-6.

### 3.4 The PARADOX — everything verified says it can't happen
- `core_blk_r` logic correct: rank-3-buffer AND flat-rank-1 rewrite give the SAME (wrong) result.
- The clean **gid-exchange test** (`GIDCHK 3Dfull` in `fesom_lifecycle_mr`) PASSES: set `g3(:,:,owned)=gid`,
  halo=−1, `exchange_elem_full`, check `g3(halo)==myList(halo)` → **WRONG=0** (all 636 halo get the correct gid).
  So slist/rlist mapping is correct, the exchange transports a clean pattern correctly.
- com is symmetric + constant (`COMCHK`): r0 recvtot=636/sendtot=605, r1 recvtot=605/sendtot=636, every call.
- owned `tr_xy` byte-exact (the SEND reads byte-exact owned via slist).
- step 1 byte-exact (same exchange, same com, different data).
- So: byte-exact source + correct mapping + symmetric com + full count → halo MUST be byte-exact. It isn't, only
  for the gradient data, only at step 2+. Logically impossible by these facts.

### 3.5 RULED OUT (don't re-investigate without new evidence)
- ❌ `core_blk_r` implementation (rank-3 buffer + flat rank-1 rewrite both fail identically).
- ❌ MPI tag cross-match with the CG (`BLK_TAG=2` unique tag — no change).
- ❌ com_elem2D_full corruption between steps (`COMCHK`: counts constant + symmetric).
- ❌ All other tracer inputs: `values`/`valuesAB`/`valuesold` (owned+halo, byte-exact), `UV`/`Kv`/`w`/`Ki`
  (byte-exact full-column at step 2), the FCT limiter (`fct_plus/minus` exchanged), vertical diffusion TDMA
  (per-column, no halo). Stale-halo probes (`stale_halo_max_nod/elem`) on all carried fields = 0.
- ❌ The eDim path: UV (eDim `exchange_elem`), UVnode (`com_nod2D` via `core_blk_r`), viscosity `U_c`
  (`com_elem2D` via `core_3D_r`) — all byte-exact at step 2.
- ❌ Same-partition: F2 and F3 read the SAME CORE2 `MeshPath`/dist_2; `dump_diff` is gid-keyed.

## 4. Leading HYPOTHESES for next session (untested or partly-tested)

1. **eXDim element ORDERING / second-layer construction differs F3 vs F2.** The gid-test passes for a clean
   pattern, but maybe F3's `com_elem2D_full` rlist/slist (or `myList_elem2D` eXDim part) lists the eXDim
   elements in a different ORDER than FESOM2, and the gid-test happens to be order-insensitive in a way the
   gradient data isn't. **TEST:** dump `partit%com_elem2D_full%{rlist,slist,rptr,sptr}` AND
   `myList_elem2D(myDim+1:nElemF)` from BOTH F3 and FESOM2 (instrument FESOM2's `init_tracers_AB` or
   `gen_modules_partitioning`), diff per rank. This is the most likely culprit and the cleanest test.
   *Note:* the M2.12a geometry gate only verified OWNED entries — the eXDim halo numbering was NEVER gated.
2. **FESOM2 PERSISTENT `tr_xy` + `MPI_TYPE_INDEXED` vs F3 LOCAL + manual pack.** FESOM2's `tr_xy` is a module
   array allocated ONCE (`oce_modules.F90:247`, `oce_setup_step.F90:927`) and exchanges via precompiled
   `MPI_TYPE_INDEXED` datatypes; F3 allocates `tr_xy` per-step (local, in `init_tracers_AB`) and packs manually.
   Maybe the eXDim layer needs the persistent array, OR FESOM2's typed exchange fills eXDim correctly where the
   manual pack doesn't. **TEST:** (a) make F3's `tr_xy` a persistent module/`tracers%work` array allocated once;
   (b) compare F3's `exchange_elem_full` against FESOM2's actual element-exchange routine line-by-line — check
   whether FESOM2's element exchange is genuinely `com_elem2D_full` or a two-stage (`com_elem2D` then a second
   hop) fill.
3. **Two-stage eXDim fill.** The second halo layer (eXDim) elements may be owned by a rank that is NOT the
   direct (eDim) neighbour, requiring a TWO-hop exchange (fill eDim from owners, then eXDim from eDim). At 2
   ranks this collapses (only one neighbour), but verify FESOM2 doesn't do an extra exchange we're missing.
   The gid-test single exchange filled all 636 though, so this is lower-priority.
4. **Compiler async-MPI optimization.** `-O3` + non-blocking `MPI_Irecv` into `rbuf` + the compiler not modeling
   the async write → stale `rbuf` read in the unpack. *Counter-evidence:* `core_2D_r`/`core_3D_r` use the same
   Irecv+Waitall+unpack pattern and ARE byte-exact. **TEST anyway:** `VOLATILE` on `rbuf`, or blocking
   `MPI_Sendrecv`, or build `mod_halo.F90` at `-O0`, and re-run the `tr_xy` dump.
5. **A specific data value triggers it (step-2 gradient vs step-1).** The step-2 values look normal (3.5e-11 is
   not denormal), but check for NaN/Inf/denormal in `tr_xy` owned at step 2 that the pack/SVML mishandles.

**Strongest lead = #1 (eXDim ordering) then #2 (persistent + typed exchange).** Both are "the F3 halo machinery
differs from FESOM2's in the untested eXDim layer," which fits ALL the evidence (clean-pattern gid-test passes,
gradient-data fails, owned byte-exact, step-2-onset is just the first step the boundary gradient is non-trivial).

## 5. REPRODUCTION (the diagnostic harness is in the tree)

Builds: F3 `./configure.sh --compiler intel --precision dp --build`; oracle `make -C
/home/a/a270088/port2/fesom2/build fesom.x` (the oracle's `oce_tracer_mod.F90` + `fesom_dump_shim.F90` have the
dump instrumentation — keep them until the bug is solved).

```bash
cd /home/a/a270088/fesom3 && source env.sh intel
RUN=/scratch/a/a270088/trxy_test; rm -rf $RUN; mkdir -p $RUN
# (1) ORACLE: 2-rank CORE2 dist_2, 2 steps, tr_xy dump
export FESOM_TRXY_DUMP=$RUN/trxy_f2
bash tools/run_lifecycle_core2.sh $RUN $RUN/life_f2 2 2
# (2) FESOM3: same
export FESOM3_MESH_DIR=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2 \
       FESOM3_IC_FILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc \
       FESOM_DUMP_FILE=$RUN/life_f3 FESOM_DUMP_MAXSTEPS=2 FESOM3_NSTEPS=2 FESOM_TRXY_DUMP=$RUN/trxy_f3
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 2 build_intel_dp/bin/fesom_lifecycle_mr
# (3) compare (owned vs halo diffs; record index 2 = step2-tr1)
python3 tools/trxy_diff.py   # see §5.1
```
Whole-gate one-liner (worst |Δ| over all substeps, full-owned): `FESOM_DUMP_ALL=1 bash
tools/run_lifecycle_gate_multirank.sh 2 3` → expect `worst |delta| ~ 1.4e-6, first at step 3 substep 9` (or
step-2 substep-15 tracers if you compare full-column).

### 5.1 The comparison script (save as `tools/trxy_diff.py`)
```python
import struct
NELEMO={0:121861,1:123360}   # myDim_elem2D per rank (CORE2 dist_2); printed by [IDXCHK]
def rd(p):
    d=open(p,'rb').read(); r=[]; o=0
    while o+8<=len(d):
        t,ne=struct.unpack_from('<2i',d,o); o+=8
        a=[struct.unpack_from('<id',d,o+12*j) for j in range(ne)]; o+=12*ne; r.append(a)
    return r
for rank in (0,1):
    myo=NELEMO[rank]
    f2=rd(f'/scratch/a/a270088/trxy_test/trxy_f2.{rank:05d}')
    f3=rd(f'/scratch/a/a270088/trxy_test/trxy_f3.{rank:05d}')
    for k,lab in [(0,'step1-tr1'),(2,'step2-tr1')]:
        r2,r3=f2[k],f3[k]; ne=len(r2)
        own=sum(1 for i in range(myo) if r2[i][1]!=r3[i][1])
        hal=sum(1 for i in range(myo,ne) if r2[i][1]!=r3[i][1])
        print(f'rank{rank} {lab}: OWNED diff={own}  HALO diff={hal}/{ne-myo}')
```

## 6. Diagnostic instrumentation currently in the tree (strip when solved)
- **F3** `src/oce/oce_tracer_mod.F90`: the `FESOM_TRXY_DUMP` block + the `[COMCHK]` print (after
  `exchange_elem_full(tr_xy)`).
- **F2 oracle** `port2/fesom2/src/oce_tracer_mod.F90`: the `FESOM_TRXY_DUMP` block (after `exchange_elem_end`).
- **F3** `src/infra/mod_dump.F90` + **F2** `port2/fesom2/src/fesom_dump_shim.F90`: the `FESOM_DUMP_ALL` mode
  (dump EVERY owned node, full column, gid-keyed — for whole-field localization).
- **F3** `src/drivers/fesom_lifecycle_mr.F90`: `[IDXCHK]` (com index bounds), `[GIDCHK 2D]`/`[GIDCHK 3Dfull]`
  (clean gid exchange self-test), `[STALE]` (per-step stale-halo probes on all carried fields).
- **F3** `src/step/mod_step_oce.F90`: `[SDBG]` per-kernel max-abs prints. **F3** `src/oce/oce_ssh_solve.F90`:
  `[CGDBG]` CG iters/residual print.
These are all env-gated or cheap and harmless; keep them for the next session, strip before the final commit.

## 7. Decision log / why each "fix attempt" was rejected
- **eDim+zero for `tr_xy`** (compute owned, zero, exchange eDim only): made step-1 drift (4.5e-5) because
  FESOM2's `tr_xy(eXDim)` is NOT zero — it IS filled each step (F2 halo changes 636/636). Reverted. Confirms
  FESOM2 fills the full eXDim halo, so F3 must too (full exchange is the right intent; its implementation drifts).
- **flat `core_blk_r` rewrite**: equivalent result (kept anyway — cleaner). Rules out rank-3 MPI buffer.
- **`BLK_TAG=2`**: no change. Rules out tag cross-match.

## 8. No-regression note
The UV eDim fix touches the MR path of `update_vel`. Before committing, re-verify the m2-mvp gate suite
(`run_step_gate_multirank.sh 2/8`, `run_step_gate.sh`, `run_stepdyn_gate_multirank.sh 2/8`,
`run_advhor_gate_multirank.sh 8`, `run_pressure_gate.sh`, `ctest`) — the eDim UV exchange should keep them
`max|Δ|=0` (eDim is a subset of the previously-used full halo, and the owned dumps depend only on eDim).

</details>
