#!/bin/bash
# M8c Step-3 BISECT — windowed DUMP_ALL localization of the day-107 (step-5095) divergence.
# The 2-year headline diverges first at step 5095 SSH_RHS, lighting up all 5 probes at once —
# the fingerprint of a GLOBAL freshwater `net` (integrate_nod_2D) spilling into water_flux at
# every node. The 5-probe set structurally CANNOT see the origin. This run dumps EVERY owned
# node (FESOM_DUMP_ALL) over a narrow step WINDOW [MINSTEP, NSTEPS] on BOTH codes, then
# tools/onset_allnode.py answers: does the first diverging substep diverge at ONE node
# (=> local origin, inspect that node) or ~ALL nodes (=> global scalar, instrument oce_fluxes)?
#
# Window default [5085, 5096] (~12 steps). All-node node dumps only (no elem stream): ~19
# field-records/node/step x 126858 nodes x 12 steps ~ 4 GB/code; onset_allnode streams F3 so the
# resident set is ~one side (~12 GB). The run still SIMULATES all NSTEPS from cold start (~10-15
# min/leg), the window only bounds what is WRITTEN.
#
#   sbatch tools/run_lifecycle_bisect_dumpall_dist864.sbatch
#   MINSTEP=5090 NSTEPS=5095 sbatch tools/...    # tighter window
#SBATCH --job-name=m8c_bisect864
#SBATCH -p compute
#SBATCH --nodes=7
#SBATCH --ntasks=864
#SBATCH --time=01:30:00
#SBATCH -A ab0995
#SBATCH -o /scratch/a/a270088/m8c_bisect864.%j.out
#SBATCH -e /scratch/a/a270088/m8c_bisect864.%j.err
set -uo pipefail
F3=/home/a/a270088/fesom3
NSTEPS="${NSTEPS:-5096}"           # simulate to just past the onset (step 5095)
MINSTEP="${MINSTEP:-5085}"         # DUMP window start (>= a few clean steps before onset)
WHICHEVP="${WHICHEVP:-0}"          # std EVP (the leg that diverged at 5095)
START_CLOCK="${START_CLOCK:-0 1 1958}"
RUN="/scratch/a/a270088/m8c_bisect864_evp${WHICHEVP}_w${MINSTEP}_${NSTEPS}"

source "$F3/env.sh" intel >/dev/null 2>&1   # KNEM fix (L35); per-leg envs are re-sourced inside the harness

# All-owned-node windowed dump on BOTH legs. srun/mpirun forward the submit environment, and the
# harness sets FESOM_DUMP_FILE/MAXSTEPS per leg WITHOUT touching these two -> both codes inherit them.
export FESOM_DUMP_ALL=1
export FESOM_DUMP_MINSTEP="$MINSTEP"

date
echo "################################################################"
echo "## M8c BISECT DUMP_ALL: NSTEPS=$NSTEPS  window=[$MINSTEP,$NSTEPS]  dist_864  whichEVP=$WHICHEVP"
echo "##   JRA55-do 1958 + zstar+TKE+GM+Redi+sw_pene + Sweeney chl ; START_CLOCK='$START_CLOCK'"
echo "################################################################"

# Run both legs. The harness ends with a 5-probe-style dump_diff over the ALL-node data; that may be
# slow/heavy but the DUMPS are written before it, so we don't let its exit status abort the localizer.
START_CLOCK="$START_CLOCK" CHL_SWEENEY=1 \
    bash "$F3/tools/run_lifecycle_jra55_gate_multirank.sh" 864 "$NSTEPS" "$WHICHEVP" "$RUN" || \
    echo "## (harness dump_diff returned non-zero — proceeding to onset_allnode on the written dumps)"

echo ""
echo "## ===== per-node origin localization (tools/onset_allnode.py) ====="
python3 "$F3/tools/onset_allnode.py" "$RUN" "$MINSTEP" "$NSTEPS"
echo ""
echo "## dump sizes (sanity):"
du -sh "$RUN"/lifef_f2.* 2>/dev/null | tail -1
du -sh "$RUN"/lifen_f3.* 2>/dev/null | tail -1
date
