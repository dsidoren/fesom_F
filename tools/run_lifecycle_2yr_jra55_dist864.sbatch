#!/bin/bash
# M8c Step 3 — 2-YEAR PRODUCTION HEADLINE byte-gate at dist_864 (7 nodes x 128 = 896 cores, 864 ranks).
# The FULL production config on the REAL forcing: JRA55-do-v1.4.0 atmosphere (gregorian + leap,
# 3-hourly, start 1958) + zstar + cvmix_TKE + Fer_GM + Redi + sw_pene + Sweeney monthly chl, native
# ice/forcing. F2 oracle vs F3, byte-gate the 5-probe x 13-substep NODE set max|delta|=0 over EVERY
# step of a 2-model-year run (1958-1959 = 35040 steps @ step_per_day=48; both years non-leap, so NO
# leap crossing — the 2920->2920 year-file rollover is in-bounds, see L50). Exercises the WHOLE
# rollover stack continuously: ~5840 atmosphere 3-hourly crossings, 24 month-end SSS+chl read-aheads,
# and the 1958->1959 year-file rollover.
#
# 5-probe dump (NOT FESOM_DUMP_ALL): only the ~5 ranks owning probe GIDs write (~228 MB/code at 35040
# steps); the global SSH elliptic solve couples every rank, so a byte divergence on ANY of the 864
# ranks reaches a probe within ~1 step (M8b already proved per-step inter-node byte-identity at 864).
# dump_diff runs on the compute node (the merged 2.28M-record dict ~few GB; 256 GB available).
#
#   sbatch tools/run_lifecycle_2yr_jra55_dist864.sbatch              # 2 yr (35040), whichEVP=0
#   NSTEPS=1488 sbatch --time=01:00:00 tools/...                     # 1-month DE-RISK leg
#   WHICHEVP=1 sbatch tools/run_lifecycle_2yr_jra55_dist864.sbatch   # mEVP leg
#SBATCH --job-name=m8c_2yr864
#SBATCH -p compute
#SBATCH --nodes=7
#SBATCH --ntasks=864
#SBATCH --time=08:00:00
#SBATCH -A ab0995
#SBATCH -o /scratch/a/a270088/m8c_2yr864.%j.out
#SBATCH -e /scratch/a/a270088/m8c_2yr864.%j.err
set -uo pipefail
F3=/home/a/a270088/fesom3
NSTEPS="${NSTEPS:-35040}"     # 2 model years (1958+1959, both non-leap) at 48 steps/day
WHICHEVP="${WHICHEVP:-0}"     # 0 = std EVP, 1 = mEVP (submit a 2nd job for the other)
START_CLOCK="${START_CLOCK:-0 1 1958}"   # cold start (override for diagnostics, e.g. "0 106 1958")
TAG="${TAG:-}"                # optional run-dir suffix for diagnostic variants
RUN="/scratch/a/a270088/m8c_2yr864_evp${WHICHEVP}_n${NSTEPS}${TAG}"
# Do NOT pre-cap the stack (L48): env.sh's `ulimit -s unlimited` would then fail under the runner's
# set -e; the gate runners set the 204800-KB dump-writer stack per leg AFTER sourcing env.sh.
source "$F3/env.sh" intel >/dev/null 2>&1   # KNEM fix (L35), platform echo

date
echo "################################################################"
echo "## M8c Step 3 2-YEAR HEADLINE: $NSTEPS steps, dist_864, whichEVP=$WHICHEVP"
echo "##   JRA55-do 1958 + zstar+TKE+GM+Redi+sw_pene + Sweeney chl"
echo "################################################################"
rc=0
START_CLOCK="$START_CLOCK" CHL_SWEENEY=1 \
    bash "$F3/tools/run_lifecycle_jra55_gate_multirank.sh" 864 "$NSTEPS" "$WHICHEVP" "$RUN" || rc=1

echo "## rollover events (rank 0):"
grep -aE "YEAR ROLLOVER|roll_monthly_clim:" "$RUN/run_f3.log" 2>/dev/null | sort | uniq -c | head -40 \
    || echo "   (none logged — rollover did NOT fire!)"
echo "## peak |uv|/|eta|/a_ice trajectory (rank 0, every ~month):"
grep -aE "fesom_lifecycle_native_mr: step" "$RUN/run_f3.log" 2>/dev/null \
    | awk 'NR==1 || NR%1488==0' | tail -30 || true
date
echo "Exit rc=$rc"
exit $rc
