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

- **Mesh:** CORE2 (~126 k nodes). **Forcing:** CORE2 NCAR interannual
  (`/pool/data/AWICM/FESOM2/FORCING/CORE2/`, years 1947+), **start year 1948**, **noleap**
  (`include_fleapyear=.false.` ⇒ `check_fleapyr→0` ⇒ 365-day years). **SSS** = monthly
  climatology (`PHC2_salx.nc`) — month-rollover, NO year-rollover. **runoff** for
  `runoff_data_source='CORE2'` is a CONSTANT total climatology read ONCE at init
  (`gen_surface_forcing.F90:1286-1289`; the monthly branch :1622 is `Dai09/JRA55`-only) — NO
  rollover at all. **chl** is constant (`chl_data_source='None'`, `chl_const=0.1`). So the only
  monthly trigger is SSS; the only year-file trigger is the 8 NCAR atmosphere fields.
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

### Task M8a-1: port `mod_clock` (clock state + advance + cold-start init)

**Files:**
- Create: `src/infra/mod_clock.F90`
- Modify: `src/params/mod_config.F90` — `runid`, `step_per_day`, `run_length`,
  `run_length_unit`, `RestartInPath/OutPath`, `include_fleapyear`, `dt` ALL EXIST
  (`mod_config.F90:16-80`); the ONE new symbol is **`r_restart`** (set by `clock_init`,
  `gen_modules_clock.F90:107-112`) — give it a home as a `mod_clock` module variable (output-
  file creation that consumes it is M9), NOT `mod_config` (keep this edit-free)

- [ ] create `src/infra/mod_clock.F90` (`module mod_clock`) transcribing `g_clock`: the saved
      state (`timeold/new`, `dayold/new`, `yearold/new`, `yearstart`, `month`, `day_in_month`,
      `fleapyear`, `ndpyr`, `num_day_in_month(0:1,12)`, `cyear*`/`cmonth`), `use mod_config` for
      `dt`/`include_fleapyear`/`runid`/`RestartInPath`
- [ ] port `clock` (advance: time→day→year→month, the `>86400._WP` / `>ndpyr` rollovers
      verbatim, `num_day_in_month(fleapyear,:)` month scan)
- [ ] port `clock_init(partit)` — read `RestartInPath//trim(runid)//'.clock'` (2 lines),
      `r_restart` test, `yearold=yearnew-1` on cold start, `check_fleapyr`, `ndpyr`, month;
      the `timenew==86400` day-start fixup (:122-125)
- [ ] port `check_fleapyr` (the `.not.include_fleapyear` early return) + `is_fleapyr`; OMIT
      `clock_finish`/`clock_newyear` (M9 restart) and the `use_transit` line (out of scope)
- [ ] **GATE (unit-level):** a tiny standalone driver advances `clock` `step_per_day` times
      from `0 1 1948` and prints `timenew/daynew/yearnew/month/day_in_month`; verify against a
      hand-computed table (day rollover at step 48, month at day 32, year at day 366→1) — proves
      the calendar arithmetic before it drives forcing

### Task M8a-2: run-length-in-years → total-steps mapping + clock-driven runloop

**Files:**
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90` (the production driver to extend) — or
  Create: `src/drivers/fesom_run_mr.F90` if a clean separation is preferred (decide in M8a-2)
- Modify: `src/infra/mod_clock.F90` (add a `clock_nsteps()` helper transcribing the oracle's
  `get_run_steps`)

- [ ] add `clock_nsteps` transcribing **`gen_model_setup.F90:334-382` (`get_run_steps`)
      verbatim** — NOT a simplified formula: `'s'`→`run_length`; `'d'`→`step_per_day*run_length`;
      `'m'`→a per-month LOOP summing `step_per_day*num_day_in_month(temp_fleapyear,temp_mon)` over
      the spanned months; `'y'`→a per-year LOOP summing `step_per_day*(365+temp_fleapyear)` with
      `check_fleapyr(temp_year,...)` re-evaluated each year. (For the noleap M8 target the `'y'`
      loop reduces to `run_length*365*step_per_day=35040`, but transcribe the loop so the M9
      decade+ 1948-2009 roadmap is leap-correct.) ⚠️ `'m'` is NEW vs `mod_config.F90:22`'s
      documented `'y, d, s'` — harmless (the driver controls the unit)
- [ ] in the driver: `call clock_init(partit)` after setup; replace the hard-coded
      `FESOM3_NSTEPS` loop bound with `nsteps = clock_nsteps(...)` (keep `FESOM3_NSTEPS`
      override for short gates); `call clock` at the TOP of each step (mirror
      `fesom_module.F90:673`, BEFORE forcing/ocean2ice)
- [ ] keep the existing forcing path UNCHANGED for now (still hard-codes day-1 rdate) — M8a is
      clock-only; the clock-derived rdate must EQUAL the current hard-coded rdate within day 1
      (it does: `daynew=1`, `timenew=n·dt` ⇒ same `rdate`), so M8a is byte-NEUTRAL within day 1
- [ ] **GATE (regression, `dist_864` + `dist_2` pre-check):** the clock-driven loop reproduces
      the existing `run_lifecycle_zstar_tke_native_gate_multirank.sh` dumps `max|Δ|=0` for the
      first 5 steps, both whichEVP — proves the clock wiring changed nothing within day 1
- [ ] re-run `ctest 13/13` + step-65 1-rank + the zstar+TKE MR gate — all `max|Δ|=0`

### Task M8b: forcing record/day rollover (`forcing_sbc_do`)

**Files:**
- Modify: `src/forcing/mod_forcing_read.F90` (extend `t_ffile`; store bracket in
  `forcing_getcoeffld`; add `forcing_sbc_do`)
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90` (thread the start-clock into the cold-start
  build; call `forcing_sbc_do` each step; replace the hard-coded day-1 rdate)
- Create: `tools/run_lifecycle_long_gate_multirank.sh` (the M8 long/rollover gate runner; clone
  of `run_lifecycle_zstar_tke_native_gate_multirank.sh` + a `START_CLOCK`/`NSTEPS` knob, with
  BOTH F3 and the F2 oracle pointed at the multi-year pool)

- [ ] **(prereq #3)** extend `t_ffile` with persisted `t_indx`, `t_indx_p1`, `nc_Ntime`; have
      `forcing_getcoeffld` STORE them (it already computes `t_indx`/`t_indx_p1` as locals,
      `mod_forcing_read.F90:276-289`) so the crossing test can read the previous bracket
- [ ] **(prereq #6)** thread `clock_init`'s `yearnew/daynew/timenew/month` into the driver's
      cold-start coef build: `rdate_cold = forcing_julday(yearnew,1,1,cal) + (daynew-1) +
      timenew/86400` (NO half-step, matching `nc_sbc_ini:643-644`), `forcing_read_grid(yearnew)`,
      SSS at `i=month` — so a non-day-1 start clock initializes forcing at the SAME point as F2
- [ ] add `forcing_sbc_do(frc, mesh, partit)` porting `gen_surface_forcing.F90:1524-1567` MINUS
      the year branch (M8c) and the noleap-leap special-case (dead for the `'none'` calendar):
      per `fld`, `rdate = forcing_julday(yearnew,1,1,cal) + (daynew-1) + timenew/86400 -
      dt/86400/2` (from `mod_clock`); crossing test `(rdate > nc_time(t_indx_p1)) .and.
      (nc_time(t_indx) < nc_time(nc_Ntime))` ⇒ `call forcing_getcoeffld`; on a wind re-trigger,
      `forcing_rotate_wind` (rotates `coef_a/coef_b`). **WIND only** — CORE2 reads no stress
      field (stress is computed from wind), so `do_rotation_stre` (:1578-1585) is dead
- [ ] wire the driver: `call clock` (M8a) → `call forcing_sbc_do` → `forcing_timeinterp` (every
      step) → bulk; delete the hard-coded `compute_native_forcing` rdate (keep the bulk/stress
      calls). NCAR time res: winds/T/q 6-hourly, radiation daily, precip monthly — the first
      6-hourly crossing fires ~step 12 (dt=1800)
- [ ] point `FESOM3_FORCING_DIR` AND the oracle `namelist.forcing` (`run_lifecycle_forced_core2.sh`
      currently hard-codes `STUB=$F2/test/input/global`, :21/:91) at the SAME multi-year pool
      `/pool/data/AWICM/FESOM2/FORCING/CORE2/` so both read identical bytes (1948/1949 present)
- [ ] **GATE (`dist_864`):** run ~48 steps from cold-start `0 1 1948` (spans several 6-hourly
      wind crossings + the daily radiation crossing); confirm `getcoeffld` actually re-fired
      (log it); byte-gate the NODE substep set `max|Δ|=0` vs F2 over all 48 steps, both whichEVP
- [ ] re-run `ctest 13/13` + the M8a regression gate — `max|Δ|=0`

### Task M8c: month-climatology + year-file rollover

**Files:**
- Modify: `src/forcing/mod_forcing_read.F90` (`forcing_sbc_do` year branch; `forcing_read_grid`
  re-entrancy)
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90` (SSS monthly read-ahead trigger; pass the
  step counter `n` as `mstep`)
- (NO `mod_forcing_other.F90` change expected — `read_other_NetCDF` already takes a month index,
  used at the cold-start `:513-517`; M8c just calls it at the read-ahead index. runoff/chl: NO
  re-read for CORE2, finding #2.)

- [ ] **(prereq #4)** make `forcing_read_grid` re-entrant: `if(allocated) deallocate` before each
      `allocate` (`:186` is bare ⇒ a 2nd call aborts). Grid is identical year-to-year ⇒ do NOT
      rebuild `idx_i/idx_j` (`forcing_build_bilin`); only `nc_time` re-bases to the new year
- [ ] **year-file branch** in `forcing_sbc_do`: port `gen_surface_forcing.F90:1510-1519` —
      `if (yearnew/=yearold)` ⇒ for each `fld` `call forcing_read_grid(frc, fld, yearnew)`
      (re-reads `nc_time` + new year offset, re-opens `u_10.<yearnew>.nc` …) + `force_newcoeff
      = .true.` ⇒ `getcoeffld` on every field that step
- [ ] **SSS month read-ahead** (finding #1, port `gen_surface_forcing.F90:1590,1594-1601`
      VERBATIM): `update_monthly_flag = (day_in_month==num_day_in_month(fleapyear,month) .AND.
      timenew==86400._WP) .OR. mstep==1`; if `surf_relax_S>0 .and. sss_data_source=='CORE2'
      .and. update_monthly_flag` ⇒ read-ahead `i=month ; if(mstep>1) i=i+1 ; if(i>12) i=1` and
      `read_other_NetCDF('SALT', month=i)`. Key off the driver step counter `n` as `mstep`
      (oracle `mstep` set at `fesom_module.F90:663`). **NOTE: fires on the LAST step of the
      month using NEXT month's SSS — NOT on the month-change step.** runoff/chl: NO re-read
- [ ] **GATE month (`dist_864`):** `fesom.clock = 0 31 1948` (day 31, late Jan) → run across into
      day 32 (Feb); the read-ahead fires on the LAST Jan step → byte-gate `max|Δ|=0` (SSS
      re-read + downstream), both whichEVP. (This gate specifically catches a wrong trigger
      phase — it would fail if SSS loads one step late)
- [ ] **GATE year (`dist_864`):** `fesom.clock = 0 365 1948` (Dec 31) → run across into day 1
      year 1949 (`yearnew` 1948→1949, opens `u_10.1949.nc` etc. — 1949 files exist in the pool)
      → byte-gate `max|Δ|=0`, both whichEVP
- [ ] re-run `ctest 13/13` + M8a/M8b gates — `max|Δ|=0`

### Task M8d: 2-year continuous headline byte-gate (the acceptance)

**Files:**
- Create: `port2/fesom2/work_m8_zstar_tke/` (the REAL oracle production work dir — see below;
  the existing `work_zstar_tke/` is JRA55-do + a frozen C-port binary, NOT reusable for a
  v2.7.3 NCAR byte-gate)
- Modify: `port2/fesom2/src/fesom_dump_shim.F90` (oracle: gate the per-substep recorder on a
  checkpoint-step list — env `FESOM_DUMP_CHECKPOINTS`; `FESOM_DUMP_ALL` already exists, :69),
  rebuild `bin/fesom.x`
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90` (F3: dump at the same checkpoint steps)
- Create: `tools/run_lifecycle_2yr_gate_multirank.sh` + `tools/job_m8_f3_2yr_864.sbatch` +
  `tools/job_m8_f2_2yr_864.sbatch` (SLURM submit + checkpoint dump-diff; SBATCH header from
  `port2/fesom2/work_zstar_tke/job_2yr_864`: 7 nodes / 864 ranks / `-A ab0995` / `-p compute`)

- [ ] **(finding #5)** build the oracle work dir from `work_core` (NOT `work_step`/`fesom_step_dump`
      — that is a one-shot prescribe-and-stop reduced-M2 shim): `which_ale='zstar'`,
      `mix_scheme='cvmix_TKE'` (+ `namelist.cvmix &param_tke`), `Fer_GM/Redi=.true.`,
      `use_sw_pene=.true.`, CORE2 NCAR pool, `include_fleapyear=.false.`, `run_length=730
      run_length_unit='d'` (= 2 yr noleap), `whichEVP` knob, fresh `fesom.clock = 0 1 1948`;
      drive the REAL `fesom.x` runloop (it dumps via `fesom_dump_shim.F90` inside
      `oce_timestep_ale`, gated by `FESOM_DUMP_MAXSTEPS`/`FESOM_DUMP_ALL`)
- [ ] add a checkpoint-step list (`FESOM_DUMP_CHECKPOINTS`) to BOTH dumpers (F2
      `fesom_dump_shim.F90` + F3 driver): dump at the end of each model MONTH (24 checkpoints
      over 2 yr), NOT every step
- [ ] **(finding #7)** set `FESOM_DUMP_ALL=1` on BOTH codes — dump EVERY owned node, not the 5
      hard-coded probe GIDs (`fesom_dump_shim.F90:61`); at `dist_864` the 5-GID default exercises
      only ~5 of 864 ranks, blind to a transport bug on the other 859 (the inter-node MPI risk
      this gate exists to catch). Volume ≈ full owned state × 24 ckpts × 864 ranks ≈ tens of GB
      per leg — defensible; confirm `dump_diff.py` iterates 864 per-rank files × 24 checkpoints
- [ ] write the SLURM batch (`ulimit -s 204800`, `source env.sh`); F3 driver `run_length=2
      run_length_unit='y'`; both codes on the SAME CORE2 NCAR pool, years 1948-1949
- [ ] submit `whichEVP=0`: F2 + F3, `dist_864`, 2 yr; byte-gate all 24 monthly checkpoints
      `max|Δ|=0` (a divergence localizes to a month, bisectable)
- [ ] submit `whichEVP=1` (mEVP): same; byte-gate `max|Δ|=0`
- [ ] if any checkpoint diverges: bisect to the first differing checkpoint → step → scalar (L29
      method); characterize (likely inter-node transport or a rollover-trigger off-by-one) + fix

### Task M8e: 2-year free-running stability sanity (FESOM3-only)

**Files:**
- Create: `tools/job_m8_stability_2yr_864.sbatch`

- [ ] run F3-only `dist_864`, 2 yr, `whichEVP=0`: assert no blow-up (the `check_blowup` ±
      guard), log peak `|uv|` / `|eta|` / `a_ice` per month; confirm a sane margin (C-port
      peak-uv ≈ 2.9 over 5 yr ⇒ expect < ~3 over 2 yr)
- [ ] run F3-only `dist_864`, 2 yr, `whichEVP=1`: same
- [ ] record the peak-`|uv|` trajectory in the gate log (confidence artifact, not a byte-gate)

### Task M8f: no-regression sweep + docs + tag

**Files:**
- Modify: `docs/HANDOFF.md`, `docs/LESSONS.md`, the memory index + implementation-state
- Modify: `docs/plans/2026-06-25-m8-long-simulations.md` (tick boxes)

- [ ] full no-regression: `ctest 13/13` (Intel + GNU dp) + step-65 1-rank + pressure-57 +
      the zstar+TKE MR gate + the M8a/b/c gates — all `max|Δ|=0`/green
- [ ] HANDOFF "Where we are" + "Next task" → M8 DONE, M9 next; new LESSONS (inter-node MPI
      byte-identity at 864 ranks; the rollover-trigger; any bisection finding)
- [ ] update `MEMORY.md` index + `project-fesom3-implementation-state.md` (the index is
      currently STALE at M3f — bring it to M8); tick this plan's boxes (memory
      `feedback-tick-plan-checkboxes`)
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
