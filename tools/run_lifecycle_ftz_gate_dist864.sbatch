#!/bin/bash
# M8c Step-3 FIX TEST — does forcing flush-to-zero (FTZ) in FESOM3 cure the day-107 byte-divergence?
# The 2-yr headline first diverges at step 5095 because a DENORMAL m_snow (~1e-309) survives in FESOM3
# where FESOM2 has exactly 0.0, flipping the if(hsn>0) ice-albedo branch -> a global SSH spill that the
# 5 probes DO see (it hits every node). This runs the standard 5-probe JRA55 gate to just past the onset
# with FESOM3_FTZ=1 (F3 forces abrupt underflow). Both legs print their FTZ-PROBE state at startup.
#   - if dump_diff reports max|delta|=0 through 5096  => candidate #1 confirmed, FTZ was the cause.
#   - if it still diverges at ~5095                   => FTZ is not (the whole) cause; bisect the code seed.
#
#   sbatch tools/run_lifecycle_ftz_gate_dist864.sbatch
#   NSTEPS=5100 FTZ=1 sbatch tools/...
#SBATCH --job-name=m8c_ftz_gate
#SBATCH -p compute
#SBATCH --nodes=7
#SBATCH --ntasks=864
#SBATCH --time=01:00:00
#SBATCH -A ab0995
#SBATCH -o /scratch/a/a270088/m8c_ftz_gate.%j.out
#SBATCH -e /scratch/a/a270088/m8c_ftz_gate.%j.err
set -uo pipefail
F3=/home/a/a270088/fesom3
NSTEPS="${NSTEPS:-5100}"            # run to just past the onset (step 5095)
WHICHEVP="${WHICHEVP:-0}"           # std EVP (the leg that diverged)
FTZ="${FTZ:-1}"                     # 1 => FESOM3 forces flush-to-zero (the candidate fix)
START_CLOCK="${START_CLOCK:-0 1 1958}"
RUN="/scratch/a/a270088/m8c_ftz_gate_evp${WHICHEVP}_ftz${FTZ}_n${NSTEPS}"

source "$F3/env.sh" intel >/dev/null 2>&1   # KNEM fix (L35)

# The candidate fix: only FESOM3 reads FESOM3_FTZ (the oracle ignores it). Standard 5-probe dump
# (NO FESOM_DUMP_ALL) — the global SSH spill at 5095 reaches the probes.
export FESOM3_FTZ="$FTZ"

date
echo "################################################################"
echo "## M8c FTZ GATE: NSTEPS=$NSTEPS  dist_864  whichEVP=$WHICHEVP  FESOM3_FTZ=$FTZ"
echo "##   JRA55-do 1958 + zstar+TKE+GM+Redi+sw_pene + Sweeney chl ; START_CLOCK='$START_CLOCK'"
echo "################################################################"

START_CLOCK="$START_CLOCK" CHL_SWEENEY=1 \
    bash "$F3/tools/run_lifecycle_jra55_gate_multirank.sh" 864 "$NSTEPS" "$WHICHEVP" "$RUN"
rc=$?

echo ""
echo "## ===== FTZ-PROBE states (both legs) ====="
grep -aE "FTZ-PROBE" "$RUN"/run_f3.log "$RUN"/run_f2.log 2>/dev/null || \
  grep -aRE "FTZ-PROBE" "$RUN" 2>/dev/null | head
echo ""
echo "## ===== gate verdict (dump_diff) ====="
grep -aE "records compared|worst|MATCH|MISMATCH|first diverging|delta" "$RUN"/*diff* "$RUN"/run_*.log 2>/dev/null | tail -15 || true
date
echo "Exit rc=$rc"
exit $rc
