#!/bin/bash
# M8b ACCEPTANCE gate at PRODUCTION SCALE (dist_864 = 7 nodes x 128 = 896 cores, 864 ranks).
# Runs the M8b long/rollover byte-gate (run_lifecycle_long_gate_multirank.sh) at 864 ranks for BOTH
# whichEVP=0 and whichEVP=1: 48 steps from cold-start '0 1 1948' (1 model day, spans the 6-hourly
# wind crossings ~steps 7/19/31/43 + the daily radiation crossing ~step 25), with BOTH the F2 oracle
# and F3 reading the full-year CORE2 pool. The gate runners auto-detect $SLURM_JOB_ID -> srun + the
# 204800-KB batch stack. SAME N on both sides => the MPI reduction tree matches (L6); byte-gate the
# per-rank NODE substep set max|delta|=0 over all 48 steps.
#
#   sbatch tools/run_lifecycle_long_gate_dist864.sbatch
#SBATCH --job-name=m8b_long864
#SBATCH -p compute
#SBATCH --nodes=7
#SBATCH --ntasks=864
#SBATCH --time=00:30:00
#SBATCH -A ab0995
#SBATCH -o /scratch/a/a270088/m8b_long864.%j.out
#SBATCH -e /scratch/a/a270088/m8b_long864.%j.err
set -uo pipefail
F3=/home/a/a270088/fesom3
NSTEPS="${NSTEPS:-48}"
# Do NOT cap the stack here: env.sh (sourced per leg by the gate runners) does `ulimit -s unlimited`,
# and pre-capping the HARD limit to 204800 makes that RAISE fail -> under the runner's `set -e` the
# whole gate aborts at `source env.sh`. The gate runners set the 204800 dump-writer stack per leg,
# AFTER sourcing env.sh (correct order). We only source env.sh here for the platform echo + KNEM env.
source "$F3/env.sh" intel >/dev/null 2>&1   # KNEM fix (OMPI_MCA_btl_vader_single_copy_mechanism=none, L35)

date
rc=0
for EVP in 0 1; do
    echo "################################################################"
    echo "## M8b dist_864 ACCEPTANCE: $NSTEPS steps, whichEVP=$EVP"
    echo "################################################################"
    bash "$F3/tools/run_lifecycle_long_gate_multirank.sh" 864 "$NSTEPS" "$EVP" \
        "/scratch/a/a270088/lifecycle_long_mr864_evp${EVP}" || rc=1
    echo "## REFRESH events (whichEVP=$EVP, rank 0):"
    grep "forcing_sbc_do: REFRESH" "/scratch/a/a270088/lifecycle_long_mr864_evp${EVP}/run_f3.log" 2>/dev/null \
        | sort | uniq -c || echo "   (none logged — rollover did NOT fire!)"
done
date
echo "Exit rc=$rc"
exit $rc
