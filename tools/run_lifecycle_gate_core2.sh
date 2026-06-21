#!/usr/bin/env bash
# M2.11c multi-step lifecycle byte-gate (CORE2 1-rank, unforced): run the REAL FESOM2
# multi-step lifecycle (run_lifecycle_core2.sh -> the built-in per-substep dump_shim) and
# FESOM3's lifecycle driver (fesom_lifecycle -> mod_step_oce::step_oce per step), compare
# the 13 NODE substeps x 5 probes over N steps (dump_diff.py --ignore-substep=2).
#
# RESULT (see LESSONS L29): the WHOLE multi-step lifecycle is byte-identical (max|delta|=0)
# on CORE2 — 195 records (13 node substeps x 5 probes x 3 steps), worst |delta| = 0,
# INCLUDING the free-surface CG d_eta (substep 9) and everything downstream (eta_n / w /
# hbar / T / S). The earlier "CG reproducibility floor" (old L28) was a BUG, not a floor:
# ssh_solve_preconditioner's off-diagonal divide auto-vectorised (divpd) in FESOM3 vs
# scalar (divsd) in the oracle, drifting the un-gated pr_values by ~1 ULP and seeding the
# CG residual at ~iter 6 (CORE2 136 iters; pi's 37 stayed below the bit). One !DIR$ NOVECTOR
# in oce_ssh_solve.F90 forces the oracle's scalar divsd -> pr_values + d_eta byte-identical.
#
#   tools/run_lifecycle_gate_core2.sh [run_dir] [nsteps]
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
RUN="${1:-/scratch/a/a270088/lifecycle_core2}"
NSTEPS="${2:-3}"

echo "[1/3] FESOM2 oracle lifecycle ($NSTEPS steps, 1-rank CORE2, unforced)"
bash "$F3/tools/run_lifecycle_core2.sh" "$RUN" "$RUN/life_f2" "$NSTEPS" | tail -3

echo "[2/3] FESOM3 lifecycle ($NSTEPS steps, 1-rank CORE2)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE"
export FESOM_DUMP_FILE="$RUN/life_f3" FESOM_DUMP_MAXSTEPS="$NSTEPS" FESOM3_NSTEPS="$NSTEPS"
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_lifecycle" 2>&1 | grep -E 'nod2D|IC T|step|done'

echo "[3/3] compare (13 NODE substeps x 5 probes x $NSTEPS steps; SW_AB substep 2 ignored)"
python3 "$F3/tools/dump_diff.py" "$RUN/life_f2.00000" "$RUN/life_f3.00000" --ignore-substep=2
