#!/usr/bin/env bash
# M8b LONG (forcing record/day rollover) PRODUCTION byte-gate (CORE2 dist_<NP>). Identical to the
# M7e zstar+TKE production gate (which_ALE='zstar' + cvmix_TKE + GM + Redi + sw_pene + native
# ice/forcing) EXCEPT both codes read the FULL-YEAR CORE2 pool instead of the day-1 stub, so a
# multi-step run crosses the 6-hourly wind records AND the daily radiation record. This is what
# exercises the NEW forcing_sbc_do path: every field whose model rdate passes its persisted bracket
# end (t_indx_p1) re-fires forcing_getcoeffld (winds also re-rotate g2r). Both sides cold-start from
# START_CLOCK (default '0 1 1948') and run NSTEPS steps; 48 steps = exactly 1 model day (48*1800s),
# spanning wind crossings at steps ~12/24/36 and the day boundary at step 48, all within January so
# the monthly SSS/CHL/runoff read fires only at mstep==1 (cold start) — pure atmosphere-rollover test.
#
#   tools/run_lifecycle_long_gate_multirank.sh [np] [nsteps] [whichEVP] [run_dir]
#   START_CLOCK="0 d 1948" env overrides the cold start (stays within 1948 for M8b).
set -euo pipefail
F3=/home/a/a270088/fesom3
CORE2POOL=/pool/data/AWICM/FESOM2/FORCING/CORE2   # full-year 6-hourly winds (1460) + daily rad (365)
NP="${1:-2}"
NSTEPS="${2:-48}"
WHICHEVP="${3:-0}"
RUN="${4:-/scratch/a/a270088/lifecycle_long_mr${NP}}"
START_CLOCK="${START_CLOCK:-0 1 1948}"

FORCING_OVERRIDE="$CORE2POOL" START_CLOCK="$START_CLOCK" \
WHICH_ALE=zstar MIX_TKE=1 SW_PENE=1 FER_GM=1 REDI=1 \
    bash "$F3/tools/run_lifecycle_fullynative_gate_multirank.sh" "$NP" "$NSTEPS" "$WHICHEVP" "$RUN"
