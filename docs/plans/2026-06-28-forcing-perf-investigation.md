# Forcing performance investigation — HANDOFF (fresh-session-ready)

**Created 2026-06-28. Self-contained: a fresh session with no memory of the prior chat can pick this up.**

## Mission

FESOM3's surface **forcing** is the ONLY component materially slower than FESOM2: **2.24× (+8.7 ms/step)**
at dist_512. Everything else (ocean, ice) is at parity or *faster* in F3. Find *why* forcing is 2× and
fix it **byte-exactly** (the project bar is `max|Δ|=0` vs the FESOM2 oracle). Closing forcing makes F3
faster than F2 overall.

## THREE operating rules (do not violate)

1. **MEASURE, DON'T GUESS.** The prior session guessed twice and was wrong both times (see "Ruled out").
   Every hypothesis must be settled by a measurement before any fix is written. Add a timer, run it, read
   the number — then act.
2. **The `compute` partition is DOWN** (Levante Sunday cooling maintenance, ~10:00–22:00; nodes in
   `drain*`). dist_512 jobs PEND with no start estimate. A few 4-node jobs *did* slip through on planned
   nodes around 11:00–11:20 today, so opportunistic dist_512 is possible but unreliable — check
   `sinfo -p compute` and `squeue --start`. The real-scale measurement may have to wait for the window to lift.
3. **Be gentle on `interactive`.** It's the only multi-rank partition up during maintenance, but
   `MaxNodes=1` (128 ranks max) and it's shared/contended. Keep jobs small and few; don't hammer it or
   the admins get annoyed. NOTE: single-node interactive **masks the forcing I/O cost** via the OS page
   cache — see "Why scale matters" — so it is NOT a faithful proxy for the dist_512 gap anyway.

## THE FINDING (corrected — fair dist_512 comparison, both ran 2026-06-28)

Same mesh + physics (CORE2 dist_512, zstar + cvmix_TKE + GM + Redi + ice, step_per_day=48 → dt=30min,
600 steps). **F2 = `/scratch/a/a270088/benchf2_d512_n600/run.log`**, **F3 = `/scratch/a/a270088/proff3_d512_n{100,600}/log`**.

| component | F3 ms/step | F2 ms/step | F3/F2 |
|-----------|-----------:|-----------:|------:|
| ocean (step_oce) | 28.66 | 30.24 | **0.95× (F3 faster)** |
| mix/pres | 2.35 | 4.57 | 0.51× |
| dyn ssh | 7.62 | 8.27 | 0.92× |
| dyn u,v,w | 4.31 | 4.38 | 0.98× |
| ice | 5.97 | 6.67 | 0.90× |
| tracer | 12.19 | 11.45 | 1.06× |
| GM/Redi | 2.18 | 1.00 | 2.18× (small: +1.2ms) |
| **forcing** | **15.74** | **7.03** | **2.24×  ← THE gap (+8.71 ms/step)** |
| LOOP total | 51.50 | 47.81 | 1.08× |

### ⚠️ UNITS GOTCHA (this caused a wrong conclusion — read carefully)

FESOM2's `___MODEL RUNTIME per task [seconds]___` block is **SECONDS for the whole run** (divide by
`nsteps` for ms/step). FESOM3's `mod_timer` report (`mean_ms/step`) is **already ms/step**. ALWAYS
convert F2 by `/nsteps*1000` before comparing. (F2 `runtime forcing: 4.2186` s / 600 = **7.03 ms/step**.)

### Where the forcing cost sits

Inside F3 forcing (15.74 ms/step), it is essentially ALL in **`TMR_FRC_SBC` (sbc/crossing) = 15.46
ms/step**; time-interp (0.018), bulk NCAR (0.26), wind/ice stress (0.004) are negligible. Forcing record
crossings fire **every ~6 steps** (JRA55 3-hourly, dt=30min; verified: 17 `REFRESH fld=1` events in 100
steps). So per-crossing cost ≈ 15.46×6 ≈ **93 ms/crossing for 8 fields** (F2 ≈ 42 ms/crossing).

## RULED OUT — do not re-investigate these

1. **Read LOCATION (rank-0 read + `MPI_Bcast` instead of all-ranks read).** Implemented, **byte-exact**
   (`max|Δ|=0`, JRA55 MR8 gate), but a PERF dead end: **neutral at dist_512** (15.97→15.75) and a **+60%
   REGRESSION at dist_128** (29.84→47.61). On a node the page cache makes the "redundant" all-ranks read
   nearly free; the ~800 KB slice broadcast is slow two-copy vader (KNEM single-copy is disabled in
   env.sh). Committed then **reverted** (`83615d7` → `d8ad992`, full reasoning in the revert message).
   ⇒ Because moving *where* the read happens changed nothing, the cost is **not the per-rank read transfer**.
2. **Build flags.** F3 and F2 Fortran flags are **identical** except F2 adds `-fPIC` (which only slows F2):
   `-O3 -r8 -i4 -fp-model precise -no-prec-div -fimf-use-svml -init=zero -fpe0 -no-prec-sqrt -ip`.
3. **"Broad slowdown" / `associate_mesh` aliasing.** A units misread briefly suggested F3 ocean was 2×
   slower and blamed F3 accessing `mesh%` directly in hot loops vs F2's local-pointer aliasing. **WRONG** —
   with correct units F3 ocean is at parity/faster. The `associate_mesh` idiom is a **dead end for perf**;
   do not pursue it.

## Leading hypothesis (UNCONFIRMED — must be measured first)

F2's `getcoeffld` (`/home/a/a270088/port2/fesom2/src/gen_surface_forcing.F90`, ~750–1020) does LESS work
per crossing than F3:
- **Persistent file handle**: F2 reads via `forcing_provider%get_forcingdata(...)` (an async provider that
  keeps the file open) — NO `nf90_open`/`nf90_close` per crossing. F3 does `nc_open_read` + `nc_close`
  **every crossing on every rank** (`src/forcing/mod_forcing_read.F90:309,318`).
- **Double-buffer slice cache**: F2 keeps `sbcdata_a`/`sbcdata_b` (per field, `SAVE`d, lines ~788–897) and
  on a crossing **reuses the previously-read `t_indx_p1` slice as the new `t_indx`** — so it reads only the
  ONE genuinely-new slice (the cached one is reused; line 926 "use the cache instead of bcast"). **F3
  re-reads BOTH brackets every crossing** (`mod_forcing_read.F90:309–318`).
- The spatial bilinear interp loop (`do ii=1,nnod`) is the same shape in both (F3 line 320; F2 961–989).

So the suspects, in order: **(a) per-crossing `nc_open`/`nc_close` overhead**, **(b) the double-read** (F3
reads 2 slices, F2 1). But note: the broadcast test moved opens 512→1 and reads 512→1 yet was neutral —
so do NOT assume (a)/(b) without the sub-timing below proving it.

## THE PLAN — measure first

1. **Sub-time inside `forcing_getcoeffld`** to split **open / read / interp**. Extend `mod_timer`
   (`src/infra/mod_timer.F90`) with e.g. `TMR_FRC_OPEN`, `TMR_FRC_READ`, `TMR_FRC_INTRP2` (children of
   `TMR_FRC_SBC` or `TMR_FORCING`; bump `NTIMER`, add to `TNAME`/`TPARENT`/`TSUBSET` — the registry is a
   pre-order tree, see the file header), OR just bracket with local `MPI_Wtime()` calls and print. Wrap:
   `nc_open_read` (line 309), the two `nc_get_slice_r4` reads (311,315), and the `do ii=1,frc%nnod` interp
   (320) in `forcing_getcoeffld`.
2. **Build**: `cmake --build build_intel_dp --target fesom_lifecycle_native_mr -- -j8` (after
   `source env.sh intel`). Incremental ~10s.
3. **Measure at dist_512** (where the gap lives) when `compute` allows — reuse
   `tools/run_prof_f3_dist512.sbatch` (two-point [100,600]). If compute is down, queue it and wait; a
   single-node interactive run will NOT show the real I/O cost (page cache).
4. **Read which sub-cost dominates**, THEN write the matching fix:
   - if `nc_open`/`nc_close`: hold a persistent handle per field (stash `ncid` in `t_ffile`, reopen only
     on year change — fname depends on year).
   - if the double-read: add F2's double-buffer cache (persist `sbc1`/`sbc2` per field in `t_atm_forcing`,
     reuse the old `t_indx_p1` as the new `t_indx`).
   - if interp: profile the loop (unlikely — nnod/rank is small at dist_512).
5. **Re-measure at dist_512** to confirm the win; target F2's ~7 ms/step.

## Byte-gating does NOT need the compute partition

The forcing byte-gate runs on the **login node** (mpirun), so you can verify correctness anytime even with
compute down:
```
tools/run_lifecycle_jra55_gate_multirank.sh 8 12 0 /scratch/a/a270088/<rundir>
# runs FESOM2 oracle + FESOM3 (NP=8, 12 steps = 2 crossings), then:
python3 tools/dump_diff.py <rundir>/lifef_f2 <rundir>/lifen_f3 --glob --ignore-substep=2
#   -> "worst |delta| = 0.000e+00" is the pass bar
```
NB the gate harness does `rm -rf` its rundir — don't put your own logs there before it runs.
`ctest` (in `build_intel_dp/`, 13/13) is the quick byte sanity but does NOT exercise the forcing read path.

## Why scale matters (don't be fooled by single-node)

The forcing read cost is a **multi-node** phenomenon: on ONE node, all 128 ranks share the OS page cache
(1 cold read + 127 cached) so all-ranks reading looks cheap. At dist_512 = 4 nodes the picture differs
(4 cold reads, more `nc_open` metadata pressure). So `interactive` (1 node) under-represents the gap.
Trust **dist_512** numbers; treat dist_128/interactive as a cheap smoke test only.

## File / data map

**F3 code:**
- `src/forcing/mod_forcing_read.F90` — `forcing_getcoeffld` (277–354): read block 309–318 (all-ranks
  `nc_open_read` + 2× `nc_get_slice_r4` + `nc_close`), interp loop 320+. `forcing_sbc_do` (394–437): the
  per-step entry; calls `getcoeffld` on crossings. `TMR_FRC_SBC` wraps `forcing_sbc_do` in the driver.
- `src/infra/mod_timer.F90` — timer registry (extend here for open/read/interp sub-timers).
- `src/drivers/fesom_lifecycle_native_mr.F90` — `compute_native_forcing` (~768): timer wrapping.

**F2 reference (the model to match):**
- `/home/a/a270088/port2/fesom2/src/gen_surface_forcing.F90` — `getcoeffld`: double-buffer cache 788–897,
  rank-0 read + bcast 902–951, persistent `forcing_provider%get_forcingdata` 915/936, interp 961–989.

**Data (today, dist_512, valid — KEEP):**
- F2 baseline: `/scratch/a/a270088/benchf2_d512_n600/run.log` (`runtime *` = seconds/600).
- F3 profile: `/scratch/a/a270088/proff3_d512_n{100,600}/log` (mod_timer `mean_ms/step`).
- dist_128 A/B + saved binaries: `/scratch/a/a270088/ab_bcast/{fesom_old,fesom_new,ab_*_n*/log}`,
  `slurm-25934633.out`. (`fesom_old` = all-ranks read; `fesom_new` = the reverted broadcast version.)

**Run scripts:** `tools/run_prof_f3_dist512.sbatch` (compute, two-point), `tools/run_prof_f3_dist128.sbatch`
(interactive), `tools/run_ab_forcing_bcast_dist128.sbatch` (two-binary A/B harness),
`tools/run_lifecycle_jra55_gate_multirank.sh` (login-node byte-gate).

**Build/env:** `source env.sh intel` then
`cmake --build build_intel_dp --target fesom_lifecycle_native_mr -- -j8`. env.sh sets
`OMPI_MCA_btl_vader_single_copy_mechanism=none` (KNEM — needed for all multi-rank; relevant if any
broadcast is ever retried — large msgs corrupt without it).

## Open questions for the user (FESOM expert) — ask early

1. Is **7 ms/step** for F2 forcing at dist_512 itself high, or expected?
2. Instinct on whether F3's 2× is the **per-crossing reopen**, the **double-read**, or the **interp** —
   to prioritize the sub-timing.

## Git state at handoff

`main`: `6158715` (mod_timer) → `4ae3718` (granular forcing sub-timers) → `83615d7` (broadcast fix) →
`d8ad992` (REVERT of the broadcast fix). Net forcing code = unchanged from pre-perf-work; the timers are
the only kept addition. All unpushed — the fix+revert pair can be squashed to a clean state if desired.
`tools/run_ab_forcing_bcast_dist128.sbatch` is untracked (commit if keeping).
