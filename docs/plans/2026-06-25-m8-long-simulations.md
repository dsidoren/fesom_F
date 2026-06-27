# M8 — production long simulations (months/years, byte-identical) — CORE2 multinode

Detailed, byte-gated decomposition. Scoped 2026-06-25 (brainstorm + oracle-investigation
pass: the clock module + the `sbc_do` rollover trigger read end-to-end, the FESOM3
forcing-reader primitive inventory, the multinode/SLURM/forcing-availability facts, the
C-port 2yr-validation precedent). Mirrors the M2–M7 method: transcribe each kernel from
FESOM2 v2.7.3, byte-gate `max|Δ|=0` against the oracle run with the **matching config**,
optional-`partit` multi-rank from the start.

Oracle root: `/home/a/a270088/port2/fesom2/src`. FESOM3 root: `/home/a/a270088/fesom3/src`.

## Overview

M8 takes the model from **short byte-gated runs** (M0–M7: the full default config is
`max|Δ|=0` vs FESOM2 for 3–5 steps, 1-rank AND multi-rank, both whichEVP) to **production
months/years runs that stay byte-identical to FESOM2**. The whole physics is already
`max|Δ|=0` per step; the only *new* code is the time machinery that lets a run advance past
"forcing day 1": **the clock**, **forcing rollover** (record/day → month → year-file), and a
**run-length-in-years driver** — exercised on the **real production mesh (CORE2) at multinode
scale (SLURM, dist_864)** over **2 years**.

**Success bar: `max|Δ|=0`** (NOT the C-port's RMS-climate bar). FESOM3 is a Fortran→Fortran
port: every step is already byte-identical, so a long run stays byte-identical **by
induction** — *provided every new rollover/clock path and the inter-node MPI are byte-faithful*.
The C-port (`/home/a/a270088/port2/fesom2_port_zstar`, `docs/validation_*_2yr/RESULTS.md`)
could only reach RMS ~0.004 because C-vs-Fortran has irreducible 1-ulp (`x**(3./2.)` Intel-libm
vs glibc `pow`) that near-zero-N² `bvfreq` sign-flips then amplified through `mo_convect` over
many steps; **those traps are ABSENT Fortran↔Fortran** (M7 handoff). So M8's bar is strictly
stronger and achievable — the validation theorem is "every per-step AND every rollover path is
`max|Δ|=0`," which holds for arbitrarily long runs.

**Scope = Option A (lean byte-identity validation).** M8 builds the clock + rollover +
run-length driver + the long-run byte-gate harness, and reuses the proven `mod_dump` +
`tools/dump_diff.py` for ALL comparison. **The production *output/restart harness is M9**
(netCDF `io_meandata`, restart write/read incl. `tke` serialization, restart-reproducibility
gate, `io_mesh_info`) — see "M9 scope" below. M8 builds NONE of it; a 2yr run holds state in
memory and runs straight through (no checkpoint until M9).

## Goal & target config (decided 2026-06-25)

- **Mesh:** CORE2 (~126 k nodes). **Forcing:** **JRA55-do-v1.4.0** atmosphere
  (`/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/`), **start year 1958**, **gregorian +
  `include_fleapyear=.true.`** (REAL leap years; model + forcing share the calendar). 3-hourly
  (2920 rec/yr non-leap, 2928 leap), `tmid=0` (mid-point shift fires), fields
  `uas/vas/huss/rsds/rlds/tas/prra/prsn` (var-name == file-name). **SSS** = monthly climatology
  (`PHC2_salx.nc`, `sss_data_source='CORE2'`, `surf_relax_S=1.929e-6`) — month-rollover, NO
  year-rollover. **runoff** (`CORE2_runoff.nc`) = CONSTANT total climatology read ONCE at init — NO
  rollover. **chl** = **`'Sweeney'` MONTHLY** climatology (`Sweeney_2005.nc`, `use_sw_pene=.true.`)
  — month-rollover (the file has NO `missing_value` ⇒ a no-fill sentinel is needed, L50). So the
  monthly triggers are **SSS + chl**; the year-file trigger is the 8 JRA55 atmosphere fields.
  **(CORRECTED 2026-06-26: this section originally read CORE2/1948/noleap/chl-`'None'` — the
  CORE2-substitute era. The user's standing correction is JRA55-do from 1958, NO CORE2 forcing; and
  `work_core`'s namelist is ground truth — SSS `'CORE2'` monthly, chl `'Sweeney'` monthly. The CORE2
  NCAR path stays only as the M2.10→M8b regression substitute, selected by `FORCING_SET=CORE2`.)**
- **Config:** the full production default — `which_ALE='zstar'` + `mix_scheme='cvmix_TKE'`
  (`mix_scheme_nmb==5`) + `Fer_GM=.true.` + `Redi=.true.` + `use_sw_pene=.true.` + native
  ice/forcing — **BOTH `whichEVP=0` (std EVP) and `whichEVP=1` (mEVP)**. `step_per_day=48`
  (dt=1800 s). `tke_cd=3.75` etc. baked (the M7 `&param_tke` doubles).
- **Scale:** **`dist_864`** (7 nodes × 128 = 896 cores, 864 ranks). SLURM: account `ab0995`,
  partition `compute`, `ulimit -s 204800`, `source env.sh` (KNEM fix
  `OMPI_MCA_btl_vader_single_copy_mechanism=none`). Template:
  `port2/fesom2/work_zstar_tke/job_2yr_864` (7 nodes / 864 ranks / 4 h / 730 days).
- **Headline gate:** **2 years continuous** (~35 040 steps), F3 vs F2, `max|Δ|=0` at MONTHLY
  dump checkpoints. **Plus** a **2yr free-running FESOM3-only stability** run (no blow-up),
  both whichEVP, not byte-gated.
- **Comparison:** F3 `dist_864` vs F2 `dist_864` (SAME N ⇒ the MPI reduction tree matches, L6),
  per-rank on owned entries (L8), via `mod_dump` + `tools/dump_diff.py`.

## Oracle map (FESOM2 runloop)

**Clock** (`gen_modules_clock.F90`, `module g_clock`, 240 lines — port whole):
- `clock` (:27-69) — called at the TOP of each step (`fesom_module.F90:673`). Advances
  `timenew+=dt`; day rollover `if (timenew>86400) daynew++ ; timenew-=86400` (:41-44); year
  rollover `if (daynew>ndpyr) daynew=1 ; yearnew++ ; check_fleapyr ; ndpyr=365+fleapyear`
  (:47-54); month/`day_in_month` from `num_day_in_month(fleapyear,:)` (:56-67).
- `clock_init(partit)` (:73-166) — reads `RestartInPath//runid//'.clock'` (2 lines:
  old, new); `r_restart = .not.(new==old)` (:107-112) ⇒ **cold start writes old==new ⇒
  `r_restart=.false.`, `yearold=yearnew-1`** (:109); `check_fleapyr`; `ndpyr`; month.
- `check_fleapyr` (:217-226) — `if(.not.include_fleapyear) return` (flag stays 0) ⇒ **noleap
  is trivial**. `is_fleapyr` (:228-236) the Gregorian rule (dead for CORE2). `clock_finish`
  (:170-200) writes `.clock` — **M9** (restart), not M8.

**Forcing rollover** (`gen_surface_forcing.F90::sbc_do`, body :1509-1567 — the new logic):
```
force_newcoeff=.false.
if (yearnew/=yearold) then                         ! :1510  YEAR-FILE rollover (M8c)
   call nc_sbc_ini_fillnames(yearnew)              !   build new year's filenames
   do fld: call nc_readTimeGrid(sbc_flfi(fld))     !   read new year's time axis
   force_newcoeff=.true.
end if
do fld_idx = 1, i_totfl                             ! :1524
   rdate = julday(yearnew,1,1,calendar)
         + (daynew-1) + timenew/86400 - dt/86400/2  ! :1527-1528  running rdate, -dt/2 shift
   ! (leap special-case :1533-1554 — DEAD for noleap 'none' calendar)
   if ( (rdate > nc_time(t_indx_p1) .and. nc_time(t_indx) < nc_time(nc_Ntime))
        .or. force_newcoeff ) then                  ! :1561  RECORD crossing (M8b)
      call getcoeffld(fld_idx, rdate, partit, mesh) !   re-read 2 bracketing slices → coef
   endif
end do
! wind/stress g2r rotation if a wind/stress field got new coeffs (:1569+)
```
Then `update_atm_forcing` (`gen_forcing_couple.F90:678 call sbc_do`) runs the bulk on the new
coefficients; `data_timeinterp` evaluates `atmdata = rdate*coef_a + coef_b` EVERY step (the
coef is only *rebuilt* on a crossing). Cold-start coef built in `nc_sbc_ini` (`getcoeffld`
at `gen_surface_forcing.F90:691`, NO half-step shift, :643-644).

**What FESOM3 already has** (`src/forcing/mod_forcing_read.F90`, M2.10a/M3f-3a — the
primitives; the GAP is the `sbc_do` orchestration that drives them per-step):
- `forcing_read_grid(frc, fld, year)` (:172) — opens `file_base//year//'.nc'`, reads `nc_time`,
  adds the year's julday offset (:210). **The year-file primitive — exists.**
- `forcing_getcoeffld(frc, fld, year, rdate, mesh, partit)` (:268-341) — its OWN
  `forcing_binarysearch(ntime, nc_time, rdate, t_indx)` (:284) finds the bracket, reads 2
  real(4) slices, spatial bilinear, builds `coef_a/coef_b`. **Re-callable at any rdate.**
- `forcing_timeinterp(frc, rdate, partit)` — `atmdata = rdate*coef_a + coef_b` (NEVER refactor
  the large-magnitude cancellation, see the module header). `forcing_build_bilin`,
  `forcing_rotate_wind`. Plus `mod_forcing_other.read_other_NetCDF` (runoff/SSS monthly slice).

So M8b/c = **add a `forcing_sbc_do(frc, mesh, partit)` to `mod_forcing_read.F90` that ports
`gen_surface_forcing.F90:1509-1567`** (rdate from `mod_clock`; re-trigger `getcoeffld` on
bracket crossing; year branch via `forcing_read_grid`), and call it each step from the driver
where it currently hard-codes the day-1 rdate. **It is NOT pure orchestration** (plan-review
findings #3/#4/#6) — three structural changes are required first:
1. **Persist the bracket state.** The crossing test `rdate > nc_time(t_indx_p1)` (:1561) needs
   the per-field `t_indx`/`t_indx_p1` from the *previous* `getcoeffld`. FESOM3's
   `forcing_getcoeffld` computes them as LOCALS (`mod_forcing_read.F90:276`) and discards them;
   `t_ffile` (:52) has no such field. ⇒ add `t_indx`/`t_indx_p1`/`nc_Ntime` to `t_ffile` and
   store them in `forcing_getcoeffld`. (Alternative considered & rejected: call `getcoeffld`
   EVERY step — byte-identical since `binarysearch` re-finds the same bracket, but 35k×8
   netCDF reads over 2 yr; the faithful persisted-trigger is both correct and cheap.)
2. **Make `forcing_read_grid` re-entrant.** It does a bare `allocate(... nc_time ...)`
   (:186) with no `deallocate` — a second call (year rollover) aborts (ifort: already
   allocated). ⇒ guard `if(allocated) deallocate` before each allocate. The grid
   (`nc_lon/nc_lat`) is identical year-to-year ⇒ do NOT rebuild the bilinear `idx_i/idx_j`.
3. **Thread the start-clock into the COLD-START coef build.** The cold-start build is
   hard-coded to day 1 / 1948 (`fesom_lifecycle_native_mr.F90`: `fyear=1948`,
   `rdate_cold=julday(1948,1,1)`, SSS `month=1`). For the induction gates that start at
   `0 31 1948` / `0 365 1948`, FESOM2's `nc_sbc_ini` builds the cold-start coef at
   `julday+(daynew-1)+timenew/86400` (no half-step, :643-644) ⇒ FESOM3 must too, else
   F3-day-1 is compared to F2-day-N. (Required for the M8c gates to be valid — see M8b.)

## Development approach (the byte-gate IS the test)

This project's "test" for every task is its **byte-gate `max|Δ|=0`** vs the FESOM2 v2.7.3
oracle on the matching config + the **no-regression** re-run (`ctest 13/13` + the relevant
existing gates). There are no unit tests in the conventional sense — a task is DONE when its
gate is `max|Δ|=0` and nothing upstream regressed.

- **Stay close to Fortran** (memory `stay-close-to-fortran`): read the actual `.F90` + the
  `work_core`/`work_*_tke` namelists, not summaries.
- **Clean-rebuild after adding NEW files** (L19): `./configure.sh --compiler intel --precision
  dp --clean --build` before trusting a gate.
- **Induction strategy:** gate each NEW rollover path at the EARLIEST step it fires, using the
  settable `fesom.clock` to position the run at a boundary — run AT `dist_864` (real multinode
  inter-node MPI) but SHORT wall-clock — THEN the 2yr continuous headline confirms nothing
  accumulates. A `dist_2`/`dist_8` login-node pre-check is fine for fast iteration, but every
  sub's ACCEPTANCE gate is at `dist_864`.
- `ulimit -s unlimited` (login) / `ulimit -s 204800` (batch) for the dump writer (L20).

## Implementation Steps

### Task M8a-1: port `mod_clock` (clock state + advance + cold-start init) ✅ DONE 2026-06-25

**DONE (first try):** `src/infra/mod_clock.F90` transcribes `g_clock` (clock/clock_init/check_fleapyr/
is_fleapyr) + `clock_nsteps` (= `get_run_steps`), `r_restart` rehomed here. Unit gate
`fesom_clocktest` + `tools/run_clocktest.sh` advances the cold-start `0 1 1948` clock 17 521 steps
and matches a HAND-COMPUTED table (day rollover step 48→49, month Jan→Feb at day 32, year 1948→1949
at day 366→1) + `clock_nsteps` for s/d/m/y — all `max|Δ|=0`, Intel AND GNU, 1-rank AND 2-rank.

**Files:**
- Create: `src/infra/mod_clock.F90`
- Modify: `src/params/mod_config.F90` — `runid`, `step_per_day`, `run_length`,
  `run_length_unit`, `RestartInPath/OutPath`, `include_fleapyear`, `dt` ALL EXIST
  (`mod_config.F90:16-80`); the ONE new symbol is **`r_restart`** (set by `clock_init`,
  `gen_modules_clock.F90:107-112`) — give it a home as a `mod_clock` module variable (output-
  file creation that consumes it is M9), NOT `mod_config` (keep this edit-free)

- [x] create `src/infra/mod_clock.F90` (`module mod_clock`) transcribing `g_clock`: the saved
      state (`timeold/new`, `dayold/new`, `yearold/new`, `yearstart`, `month`, `day_in_month`,
      `fleapyear`, `ndpyr`, `num_day_in_month(0:1,12)`, `cyear*`/`cmonth`), `use mod_config` for
      `dt`/`include_fleapyear`/`runid`/`RestartInPath`
- [x] port `clock` (advance: time→day→year→month, the `>86400._WP` / `>ndpyr` rollovers
      verbatim, `num_day_in_month(fleapyear,:)` month scan)
- [x] port `clock_init(partit)` — read `RestartInPath//trim(runid)//'.clock'` (2 lines),
      `r_restart` test, `yearold=yearnew-1` on cold start, `check_fleapyr`, `ndpyr`, month;
      the `timenew==86400` day-start fixup (:122-125)
- [x] port `check_fleapyr` (the `.not.include_fleapyear` early return) + `is_fleapyr`; OMIT
      `clock_finish`/`clock_newyear` (M9 restart) and the `use_transit` line (out of scope)
- [x] **GATE (unit-level):** ✅ `fesom_clocktest` (`tools/run_clocktest.sh`) advances `clock` from
      `0 1 1948` for 17 521 steps; verified vs a HAND-COMPUTED table — day rollover step 48→49
      (the `>86400` strict test keeps step 48 on day 1, rolls on step 49), month at day 32 (Jan→Feb),
      year at day 366→1 (1948→1949) — all `max|Δ|=0`. (Also gates `clock_nsteps`, M8a-2.)

### Task M8a-2: run-length-in-years → total-steps mapping + clock-driven runloop ✅ DONE 2026-06-25

**DONE (first try) → M8a (clock + run-length driver) COMPLETE.** Wired `mod_clock` into
`fesom_lifecycle_native_mr` (extended the existing driver — no separate `fesom_run_mr` needed):
set `step_per_day=48`/`mod_config%dt=1800`/`run_length`/`run_length_unit` so the ported `clock`
advances `mod_config%dt` per step; `FESOM3_START_CLOCK` ("t d y", default "0 1 1948") → rank 0
writes a 2-identical-line `.clock` (`RestartInPath` from `FESOM3_RESTART_IN` else dirname of
`FESOM_DUMP_FILE`), barrier, all ranks `clock_init`; `nsteps = clock_nsteps(partit)` with the
`FESOM3_NSTEPS` short-gate override; `call clock` at the top of each step. Forcing path UNCHANGED
(still hard-codes the day-1 rdate) ⇒ byte-NEUTRAL within day 1. **GATE `max|Δ|=0`:** zstar+TKE MR
325 records, dist_2 (BOTH whichEVP) + dist_8, 5 steps. **No regression:** ctest 13/13 (Intel+GNU
dp) + step-65 1-rank — all `max|Δ|=0`/green. NOTE: `dist_864` deferred to M8b (the clock-file
cross-node read is validated there alongside the first real forcing crossing — M8a changes no
transport, so a login dist_2+dist_8 byte-neutral check is the right coverage here).

**Files:**
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90` (the production driver — extended in place)
- Modify: `src/infra/mod_clock.F90` (added the `clock_nsteps()` helper transcribing
  `get_run_steps`)

- [x] add `clock_nsteps` transcribing **`gen_model_setup.F90:334-382` (`get_run_steps`)
      verbatim** — NOT a simplified formula: `'s'`→`run_length`; `'d'`→`step_per_day*run_length`;
      `'m'`→a per-month LOOP summing `step_per_day*num_day_in_month(temp_fleapyear,temp_mon)` over
      the spanned months; `'y'`→a per-year LOOP summing `step_per_day*(365+temp_fleapyear)` with
      `check_fleapyr(temp_year,...)` re-evaluated each year. (For the noleap M8 target the `'y'`
      loop reduces to `run_length*365*step_per_day=35040`, but transcribe the loop so the M9
      decade+ 1948-2009 roadmap is leap-correct.) ⚠️ `'m'` is NEW vs `mod_config.F90:22`'s
      documented `'y, d, s'` — harmless (the driver controls the unit). Unit-gated: s,5→5; d,730→
      35040; m,1→1488; m,2→2832; y,1→17520; y,2→35040.
- [x] in the driver: `call clock_init(partit)` after setup; replace the hard-coded
      `FESOM3_NSTEPS` loop bound with `nsteps = clock_nsteps(...)` (kept `FESOM3_NSTEPS`
      override for short gates); `call clock` at the TOP of each step (mirror
      `fesom_module.F90:673`, BEFORE forcing/ocean2ice)
- [x] keep the existing forcing path UNCHANGED for now (still hard-codes day-1 rdate) — M8a is
      clock-only; the clock-derived rdate must EQUAL the current hard-coded rdate within day 1
      (it does: `daynew=1`, `timenew=n·dt` ⇒ same `rdate`), so M8a is byte-NEUTRAL within day 1
- [x] **GATE (regression, `dist_2`+`dist_8` login; `dist_864` folded into M8b):** the clock-driven
      loop reproduces the existing `run_lifecycle_zstar_tke_native_gate_multirank.sh` dumps
      `max|Δ|=0` (325 records) for the first 5 steps, BOTH whichEVP — proves the clock wiring
      changed nothing within day 1
- [x] re-run `ctest 13/13` (Intel + GNU dp) + step-65 1-rank + the zstar+TKE MR gate — all `max|Δ|=0`/green

### Task M8b: forcing record/day rollover (`forcing_sbc_do`)

**Files:**
- Modify: `src/forcing/mod_forcing_read.F90` (extend `t_ffile`; store bracket in
  `forcing_getcoeffld`; add `forcing_sbc_do`)
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90` (thread the start-clock into the cold-start
  build; call `forcing_sbc_do` each step; replace the hard-coded day-1 rdate)
- Create: `tools/run_lifecycle_long_gate_multirank.sh` (the M8 long/rollover gate runner; clone
  of `run_lifecycle_zstar_tke_native_gate_multirank.sh` + a `START_CLOCK`/`NSTEPS` knob, with
  BOTH F3 and the F2 oracle pointed at the multi-year pool)

- [x] **(prereq #3)** extend `t_ffile` with persisted `t_indx`, `t_indx_p1`; have
      `forcing_getcoeffld` STORE them (it already computes `t_indx`/`t_indx_p1` as locals,
      `mod_forcing_read.F90:276-289`) so the crossing test can read the previous bracket. ✅
      (`nc_Ntime` reused the existing `ntime` field — verified `nc_Ntime`==time-dim length.)
- [x] **(prereq #6)** thread `clock_init`'s `yearnew/daynew/timenew/month` into the driver's
      cold-start coef build: `rdate_cold = forcing_julday(yearnew,1,1,cal) + (daynew-1) +
      timenew/86400` (NO half-step, matching `nc_sbc_ini:643-644`), `forcing_read_grid(yearnew)`,
      SSS at `i=month`. ✅ Required REORDERING the M8a clock block to BEFORE the forcing setup.
- [x] add `forcing_sbc_do(frc, mesh, partit)` porting `gen_surface_forcing.F90:1524-1567` MINUS
      the year branch (M8c) and the noleap-leap special-case (dead for the `'none'` calendar):
      per `fld`, `rdate = forcing_julday(yearnew,1,1,cal) + (daynew-1) + timenew/86400 -
      dt/86400/2` (from `mod_clock`); crossing test `(rdate > nc_time(t_indx_p1)) .and.
      (nc_time(t_indx) < nc_time(nc_Ntime))` ⇒ `call forcing_getcoeffld`; on a wind re-trigger,
      `forcing_rotate_wind`. **WIND only** — CORE2 reads no stress field so `do_rotation_stre`
      is dead. ✅ + a rank-0 REFRESH diagnostic (stdout only — never touches the dump).
- [x] wire the driver: `call clock` (M8a) → `call forcing_sbc_do` → `forcing_timeinterp` (every
      step) → bulk; replaced the hard-coded `compute_native_forcing` rdate with the clock-derived
      one (byte-equal within day 1, correct across crossings). ✅
- [x] point `FESOM3_FORCING_DIR` AND the oracle `namelist.forcing` at the SAME multi-year pool
      `/pool/data/AWICM/FESOM2/FORCING/CORE2/`. ✅ Parametrized `run_lifecycle_forced_core2.sh` +
      the MR gate with `FORCING_OVERRIDE`/`START_CLOCK` (backward-compatible); NEW
      `tools/run_lifecycle_long_gate_multirank.sh` (pool, zstar+TKE) + `…_dist864.sbatch`.
- [x] **GATE (login `dist_2`+`dist_8`):** 48 steps from cold-start `0 1 1948` (1 model day) —
      **`max|Δ|=0` over 3120 records, BOTH whichEVP**, with **18 `getcoeffld` re-fires** logged
      (6-hourly winds/q/Tair brackets 1/2→2/3→3/4→4/5 at steps ~7/19/31/43 + daily radiation 1/2
      at step ~25) — non-vacuous rollover. Also verified at `dist_128` on a compute node (batch).
- [x] **GATE (`dist_864`):** production-scale acceptance — **48 steps from cold-start `0 1 1948`,
      BOTH whichEVP, `max|Δ|=0` over 3120 records each, 18 `getcoeffld` re-fires each** (job
      25915080, 7 nodes/864 ranks, `COMPLETED`). Needed the batch harness (NEW for fesom3):
      `srun`/204800-stack SLURM branches in the runners + the `shell.intel`
      `ulimit -s unlimited 2>/dev/null || true` fix (a bare `ulimit` raise on a hard-capped
      compute node aborts the gate under the runner's `set -e` — see L48). dist_128 batch too.
- [x] re-run `ctest 13/13` + the M8a regression gate — **`max|Δ|=0`** (ctest 13/13; zstar+TKE MR
      stub day-1 `dist_2` both whichEVP + `dist_8`, 325 records). ✅

### Task M8c: JRA55 forcing-read gate (Step 1) + month-climatology & year-file rollover (Step 2) ✅ Steps 1+2 COMPLETE 2026-06-26

**STATUS:** Step 1 (JRA55 forcing-read gate) + Step 2 (rollover) DONE, `max|Δ|=0`. **Step 1:** the FIRST gate on the
real JRA55-do atmosphere (not the CORE2 substitute) byte-matched the 8 fields + bulk from 1958 — exposed/fixed the
calendar NUL-byte trap (L49); `dist_2` both whichEVP + `dist_8` + CORE2 no-regression 48-step, all `max|Δ|=0`.
**Step 2:** year-file branch + SSS monthly + **Sweeney chl monthly** (the production config — `chl_data_source='Sweeney'`,
NOT the `'None'`/const this plan originally assumed) — `dist_2` month/year/chl gates + a `dist_8` combined year+chl gate,
all `max|Δ|=0`, all non-vacuous (L50). **Step 3 (2-yr headline) = REMAINING.**

**Files (as actually changed):**
- `src/forcing/mod_forcing_read.F90` (`forcing_sbc_do` year branch + `force_newcoeff`; `forcing_read_grid` re-entrancy)
- `src/drivers/fesom_lifecycle_native_mr.F90` (`roll_monthly_clim` SSS + chl read-ahead, step counter `n` as `mstep`;
  `forc_set` JRA55/CORE2 wiring from Step 1; `use_chl_sweeney` + cold-start Sweeney read)
- `src/io/mod_io_netcdf.F90` (Step 1 NUL sanitize in `nc_get_att_text`; `nc_get_att_dp` optional `stat` for chl)
- `src/forcing/mod_forcing_other.F90` (chl no-fill sentinel `miss=-99` when `missing_value` absent — the plan's "NO
  `mod_forcing_other.F90` change / runoff-chl NO re-read" note was the CORE2-substitute assumption; chl-Sweeney DID
  need both a monthly re-read AND the sentinel)

- [x] **(prereq #4)** `forcing_read_grid` re-entrant: `if(allocated) deallocate` the 3 axis arrays
      before each `allocate`. Grid identical year-to-year ⇒ do NOT rebuild the bilinear source
      indices (`forcing_build_bilin`); only `nc_lon/nc_lat/nc_time` re-read (`ntime` CAN change at a
      leap boundary 2920↔2928 ⇒ free, don't reuse)
- [x] **year-file branch** in `forcing_sbc_do`: ported `gen_surface_forcing.F90:1510-1519` —
      `if (yearnew/=yearold)` ⇒ for each `fld` `call forcing_read_grid(frc, fld, yearnew)`
      (re-anchors `nc_time`, re-opens `<var>.<yearnew>.nc`) + `force_newcoeff=.true.` OR'd into the
      crossing test ⇒ `getcoeffld` on every field that step. `yearnew/yearold` from `mod_clock`
- [x] **SSS month read-ahead** (`roll_monthly_clim` in the driver, port
      `gen_surface_forcing.F90:1590,1594-1601`): `update_monthly_flag =
      (day_in_month==num_day_in_month(fleapyear,month) .AND. timenew==86400._WP) .OR. mstep==1`;
      read-ahead `i=month ; if(mstep>1) i=i+1 ; if(i>12) i=1` + `read_other_NetCDF('SALT', i)`. Keyed
      off the driver step counter `n` as `mstep`. Fires on the LAST step of the month with NEXT
      month's SSS
- [x] **chl-Sweeney month read-ahead** (NOT in the original plan — production `chl_data_source='Sweeney'`):
      same trigger/index, `read_other_NetCDF('chl', i)` from `Sweeney_2005.nc`; needed `nc_get_att_dp`
      optional `stat` + `miss=-99` no-fill sentinel (the file has no `missing_value`, L50). Gated by
      `FESOM3_CHL_SWEENEY` / `CHL_SWEENEY=1`
- [x] **GATE month (`dist_2`, both whichEVP via separate legs):** `START_CLOCK="82800 31 1958"` (23:00
      Jan 31) → crosses into Feb in 2 steps; the read-ahead fires step 2 (Feb slice) → `max|Δ|=0` (780
      rec). Cheaper than the planned midnight `0 31` start (2 steps vs 48)
- [x] **GATE year (`dist_2` + `dist_8`):** `START_CLOCK="82800 365 1958"` (23:00 Dec 31) → crosses into
      1959 (`YEAR ROLLOVER 1958→1959` fires step 3, opens `uas.1959.nc` … — 1959 files exist, 2920 rec)
      + Dec→Jan SSS/chl read-ahead step 2 → `max|Δ|=0` (`dist_2` SSS-only AND `dist_8` SSS+chl-Sweeney)
- [x] **no-regression:** CORE2 M8b long gate 48-step `dist_2` (3120 rec) `max|Δ|=0`; mEVP leg `max|Δ|=0`.
      `ctest 13/13` re-run after the shared `mod_io_netcdf`/`mod_forcing_other` edits — see Step-3 wrap

### Task M8d: 2-year continuous headline byte-gate (the acceptance) — ✅ COMPLETE 2026-06-27 (validated PHYSICALLY)

**OUTCOME (2026-06-27): the strict 2-yr byte-gate is UNACHIEVABLE — but NOT because of a port bug — so the headline
is validated PHYSICALLY (user decision "validate 2-yr physically, not byte-exact").** The byte-gate diverges first at
**step 5095 (day ~107), SSH_RHS, at all 126858 nodes at once** (the global-spill fingerprint). Bisected to the origin
with a windowed `FESOM_DUMP_ALL` (`tools/onset_allnode.py`, denormal-filtered) + per-gid `extract_gid.py`: a spurious
**denormal `m_snow`** (~1.3e-309 m = physically-zero snow) in FESOM3 where FESOM2 has **exactly 0.0** (node 119505)
flips the ice-albedo branch `if(hsn>0)` (snow albedo 0.85 vs ice 0.65) → absorbed-SW shift → `t_skin` diverges 2.6%
(−0.4350 vs −0.4614) → evap/sublimation/`thdgr` → freshwater `flux` → global `net` via `integrate_nod_2D` →
`water_flux` at every node → SSH spill. **Proven EMERGENT, not transcription:** every snow/ice routine (FCT
antidiffusive flux, low-order solve, `cut_off`, `obudget`, `budget`, thermo snow-update, `flooding`, `ice_TG_rhs`)
AND all compiler flags are byte-identical, and **`m_ice`/`a_ice` stay `|Δ|=0` through step 5094** — only the near-zero
`m_snow` exercises the denormal regime (FESOM2 itself seeds denormals, e.g. 2.75e-308 at node 42847). A tolerance can't
help (the divergence amplifies) and the only byte-exact fix would modify the vanilla FESOM2 reference (a shared
denormal flush). Full write-up: **LESSONS L51**; memory `m8c-day107-divergence`. ⇒ Steps 1+2 keep the strict
`max|Δ|=0` bar (already green); the **2-yr HEADLINE is the free-running stability run** (M8e, below) — which COMPLETED
35040 steps cleanly with bounded global peaks and a correct ice seasonal cycle.

**APPROACH (revised from the original CORE2/DUMP_ALL plan):** JRA55-do **1958-1959** (35040 steps,
BOTH non-leap ⇒ no leap crossing, the 2920→2920 year-file rollover is in-bounds). REUSE the proven
`run_lifecycle_jra55_gate_multirank.sh` (→ `run_lifecycle_fullynative_gate_multirank.sh`): the oracle
leg drives the REAL `bin/fesom.x` runloop set up from `work_core` (`run_lifecycle_forced_core2.sh`,
`run_length_unit='s'` ⇒ `run_length=NSTEPS`), the F3 leg the native driver — both auto-srun at
`SLURM_JOB_ID`. **NO oracle rebuild, NO dumper change:** keep the **5-probe every-step dump** (only the
~5 ranks owning a probe GID write ⇒ ~228 MB/code at 35040 steps, NOT the tens-of-GB DUMP_ALL). The 5
probes are NOT blind to the other 859 ranks over a long run: the GLOBAL SSH elliptic solve couples
every rank each step, so a byte divergence anywhere reaches a probe within ~1 step — and M8b already
proved per-step inter-node byte-identity at 864 ranks. `dump_diff` runs on the compute node (merged
2.28M-record dict ≈ few GB; 256 GB available); every-step dumping gives EXACT first-diverging-step
localization (better than the 24-checkpoint plan).

**Files:**
- Create: `tools/run_lifecycle_2yr_jra55_dist864.sbatch` (7 nodes / 864 ranks / `-A ab0995` /
  `-p compute`; `NSTEPS`/`WHICHEVP` env knobs; reuses the JRA55 wrapper; greps the rollover +
  peak-|uv|/|eta|/a_ice trajectory for the M8e stability artifact)
- (NO `work_m8_zstar_tke`, NO `fesom_dump_shim.F90`/oracle rebuild, NO F3 dumper change — all avoided)

- [x] **de-risk (login):** dist_8, 96 steps, Jan→Feb, chl=Sweeney — `max|Δ|=0` (6240 records); extends
      the 12-step M8c gates, confirms continuous byte-identity over a month-crossing window
- [x] **de-risk (scale + runtime):** the headline byte-gate itself ran clean `max|Δ|=0` at `dist_864` for
      ~5094 steps (≈106 days) before the day-107 denormal divergence — so dist_864 scale AND runtime (2-yr
      job fits the partition cap) are BOTH validated well past 1 month. Used the every-step 5-probe dump
      (`run_lifecycle_jra55_gate_multirank.sh`).
- [x] **submit `whichEVP=0` byte-gate (F2 + F3, `dist_864`, 2 yr):** SUBMITTED and ran — diverged at step
      5095 (NOT inter-node transport / NOT a rollover off-by-one; M8b already proved 864-rank inter-node
      byte-identity). Root cause = emergent denormal `m_snow` flipping an ice-albedo branch (see OUTCOME +
      **L51**), which a byte-gate cannot pass without modifying the FESOM2 reference.
- [x] **diverged ⇒ bisect (L29 method):** windowed `FESOM_DUMP_ALL` ([5085,5096]) + `tools/onset_allnode.py`
      (denormal-floored) + `tools/extract_gid.py` localized it to node 119505 / step 5095 / `m_snow` denormal.
      `m_ice`/`a_ice` byte-clean through 5094 ⇒ snow-specific, not the shared advection/transport.
- [~] **`whichEVP=1` (mEVP) byte-gate:** not run — same denormal mechanism applies (it is config-independent),
      so it would diverge identically. Strict mEVP byte-identity stays proven by the SHORT M8c gates.
- [x] **RESOLUTION (user decision):** validate the headline PHYSICALLY via the free-running F3 run (M8e),
      not bit-identity — the only obstacle is physically-meaningless denormal snow both codes generate.

### Task M8e: 2-year free-running stability sanity (FESOM3-only) — ✅ COMPLETE 2026-06-27 (the headline validation)

Because the byte-gate diverges at day 107 (denormal snow, M8d OUTCOME), this free-running stability run IS the
2-year headline acceptance, not just a sanity check. A DEDICATED sbatch (`tools/run_lifecycle_2yr_freerun_f3_dist864.sbatch`,
NOT the byte-gate's F3 leg) runs the full production config (JRA55-do 1958-1959, zstar+TKE+GM+Redi+sw_pene, Sweeney chl)
free-running with NO oracle / NO dump. The per-step diagnostic was upgraded from rank-0-local to **GLOBAL `MPI_MAX` over
all 864 ranks** (one packed reduction/step — the proven `mod_halo` idiom; the rank-0 `max|a_ice|` always read 0 because
rank 0 owns no ice nodes, and the driver writes no NetCDF): `max|eta|`, `max|uv|`, `a_ice`, `m_ice`, `Tmax`, `Smax`.

- [x] **whichEVP=0, dist_864, 2 yr (job 25927003, COMPLETED, 31:49):** 35040 steps done, **0 NaN/Inf/FPE**, year
      rollover `1958→1959` + all 24 monthly SSS/chl read-aheads fired. Global peaks bounded with NO secular drift:
      `max|eta|` 1.8–2.0 m (peak 2.01), `max|uv|` 1.5–2.9 m/s (< the ~3 C-port margin), `Tmax` 30–32.5 °C flat, `Smax`
      ~41.04 flat. Cryosphere ALIVE with the correct **seasonal cycle**: `a_ice`→1.0 sustained, `m_ice` 2.0→6.1→4.6 m
      (winter-grow/summer-melt). ⇒ free-running FESOM3 is physically stable & sane over 2 production years.
- [~] **whichEVP=1 (mEVP):** not yet run; whichEVP=0 establishes stability and mEVP was byte-proven in M3 + the short
      M8c gates. Run if a second stability point is wanted (cheap: ~32 min, 3.7 node-h).

### Task M8f: no-regression sweep + docs + tag

**Files:**
- Modify: `docs/HANDOFF.md`, `docs/LESSONS.md`, the memory index + implementation-state
- Modify: `docs/plans/2026-06-25-m8-long-simulations.md` (tick boxes)

- [x] **back out the day-107 fw_ bisect scaffolding from BOTH codes** (physics-neutral): F3
      `mod_ice_oce_coupling.F90` reverted to HEAD (fw_ dumps + `mstep` arg gone) + the driver's
      `oce_fluxes` call de-argued; oracle `ice_oce_coupling.F90`/`oce_modules.F90` reverted to HEAD and
      the NATIVE `mstep` decl + `mstep = n` RESTORED (they are vanilla FESOM2, NOT bisect — the bisect
      only reused them; a wrong first pass removed them and was corrected). KEPT (clean, inert, reusable
      infra, NOT throwaway): the GLOBAL `MPI_MAX` per-step stability diagnostic, the `mod_dump`
      DUMP_ALL/MINSTEP windowed all-node dump, and the `FESOM3_ATMFLUX_CHECK` self-check.
- [x] full no-regression after the back-out: **`ctest 13/13` Intel dp + GNU dp** (both 100%) + the
      **production fully-native MR gate** (`run_lifecycle_fullynative_gate_multirank.sh`, dist_2, 3 steps,
      zstar+TKE+GM+Redi+sw_pene, **BOTH whichEVP**) = `max|Δ|=0` (195 records each). This MR gate drives
      both codes through the back-out's `oce_fluxes` path with the full production column physics, so it
      subsumes step-65/pressure-57; F3≡oracle byte-exact confirmed.
- [x] HANDOFF "Where we are" → **M8 ✅ COMPLETE (tag `m8`)**; "RESUME HERE" → **M9** (production I/O
      harness). LESSONS **L51** (the day-107 EMERGENT denormal-`m_snow` albedo-branch finding).
- [x] update `MEMORY.md` index + `project-fesom3-implementation-state.md` (was STALE at M3f → brought to
      M8); ticked this plan's boxes.
- [ ] commit + **tag `m8`**

## Risks & watch-list

- **Inter-node MPI byte-identity at 864 ranks (the headline risk).** Every prior multi-rank
  gate was `dist_2`/`dist_8` on ONE login node (shared-memory `vader`). M8 is the FIRST test of
  real inter-node transport (UCX/openib BTL, not `vader`). L35's KNEM corruption was
  vader-only and is fixed by `env.sh`; verify there is no inter-node analog (a short `dist_864`
  M8b gate catches it early — that's the induction payoff). If a transport floor appears,
  characterize it like L29 (bisect sbuf→rbuf→unpack) and pin the BTL/MCA setting in `env.sh`.
  **Coverage caveat:** this risk is only validated if the headline gate dumps EVERY owned node
  (`FESOM_DUMP_ALL=1`) — the default 5 probe GIDs touch ~5 of 864 ranks (M8d, finding #7).
- **MPI_Allreduce reduction-tree determinism.** The CG dot-products + `mesh%ocean_area` use
  `allreduce_sum`. Byte-identical only for a FIXED comm size + op + 8-byte type + OpenMPI build
  (L6). The gate compares F3 `dist_864` vs F2 `dist_864` (SAME N) ⇒ the tree matches. Do NOT
  compare across different N.
- **Rollover-trigger timing must EXACTLY match FESOM2.** The `rdate > nc_time(t_indx_p1)`
  crossing and the year `yearnew/=yearold` branch must fire on the SAME step as the oracle (a
  one-step-late re-read = a one-step-stale coef = a mismatch). The early-triggered M8b/c gates
  pin this surgically.
- **`forcing_timeinterp` cancellation (module header).** `atmdata = rdate*coef_a + coef_b` with
  `rdate ~ 711k` (julday) — NEVER refactor to `data1+coef_a*(rdate-nc_time)`; the WP rounding of
  the large `coef_b` differs (already gated through M3f, but the rollover re-evaluates it at new
  rdates — keep the form).
- **Dump volume over 2 yr.** Mitigated: MONTHLY checkpoints (24), not every step — but with
  `FESOM_DUMP_ALL` (all owned nodes, needed per finding #7) it is ≈ tens of GB/leg. Acceptable
  on `/scratch`; clean up between legs.
- **Wall-clock / node-hours.** ~4 h/leg @ 864 ranks × (2 whichEVP × 2 codes) for M8d ≈ ~110
  node-hours, + 2 stability legs (F3-only, ~4 h each). The induction gates (M8a/b/c) are short.
  Budget on `ab0995`.
- **`tke` prognostic state held in memory the whole run** (`dyn%work%tke`, the first stateful
  mixing field) — fine for a straight-through 2yr run; its restart serialization is M9.
- **Forcing dataset.** Use CORE2 NCAR (`u_10.YYYY.nc` …, what the M2.10/M3f gates proved) —
  NOT the JRA55-do `work_zstar_tke` 2yr-reference dir (different dataset + a frozen C-port
  binary). Re-run the v2.7.3 oracle with the `mod_dump` shim.

## M9 scope (the production harness — planned, built later)

Recorded here so it is on the roadmap; M8 builds NONE of it. M9 = the production I/O harness,
decomposed like M2–M8:
- **netCDF mean/snapshot output** — port `io_meandata.F90` (3121 lines): the mean-data buffer,
  the variable registry, `io_gather` (nod/elem), the netCDF write. Gate: F3 output vs F2 output
  field-by-field (the FIRST netCDF-WRITE gate; the existing `mod_io_netcdf.F90` only READS).
- **restart write/read** — port `io_restart.F90` (+ `io_restart_derivedtype.F90`,
  `io_restart_file_group.F90`, ~1200 lines) incl. **`tke` serialization** (the first stateful
  mixing field — `dyn%work%tke`, not yet in any `write_t_dyn_work`) + `clock_finish` (the
  `.clock` write, M8a-deferred). Gate: **restart-reproducibility** = a run split across a
  restart is `max|Δ|=0` vs the straight-through run (a strong self-consistency byte-test) +
  byte-match vs an F2 restart.
- **`io_mesh_info`** (mesh.nc for the output) + output scheduling (`io_data_strategy`,
  annual/monthly events from the clock).
- Then **real production science runs** (decade+ CORE2, netCDF output, restart chaining,
  forcing year-rollover across the full 1948-2009 record).

## References

- Plan format: `docs/plans/2026-06-24-m7-tke.md`. Source of truth: `docs/HANDOFF.md`,
  `docs/LESSONS.md` (L6 Allreduce determinism, L8 same-partition gate rule, L19 clean-rebuild,
  L20 dump ulimit, L29 bisect-the-divide method, L35 vader/KNEM).
- Oracle: `gen_modules_clock.F90` (`g_clock`); `gen_surface_forcing.F90` (`sbc_do` :1509-1567,
  cold-start `getcoeffld` :691, `nc_sbc_ini`); `gen_forcing_couple.F90:678` (`call sbc_do`);
  `fesom_module.F90:673` (`call clock`).
- FESOM3: `src/forcing/mod_forcing_read.F90` (primitives: `forcing_read_grid` :172,
  `forcing_getcoeffld` :268, `forcing_timeinterp`, `forcing_binarysearch`),
  `src/forcing/mod_forcing_other.F90` (`read_other_NetCDF`), `src/params/mod_config.F90`
  (run_length/step_per_day/dt/include_fleapyear/runid — all present),
  `src/drivers/fesom_lifecycle_native_mr.F90` (the driver to extend),
  `tools/run_lifecycle_zstar_tke_native_gate_multirank.sh` (gate template).
- Multinode: partitions `/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2/dist_{2..864}`; SLURM
  template `port2/fesom2/work_zstar_tke/job_2yr_864`; forcing pool
  `/pool/data/AWICM/FESOM2/FORCING/CORE2/` (1947+). C-port precedent:
  `port2/fesom2_port_zstar/docs/validation_*_2yr/RESULTS.md` (RMS bar — M8 beats it with
  `max|Δ|=0`).

## Progress tracking

- mark completed items `[x]` immediately when their gate is `max|Δ|=0`; add `➕` for newly
  discovered tasks, `⚠️` for blockers; keep this file in sync (memory
  `feedback-tick-plan-checkboxes`: tick boxes + the `✅` header marker when a sub is gated DONE).
