#!/usr/bin/env bash
# M4c GM-bolus lifecycle byte-gate (CORE2 1-rank, unforced, Fer_GM=T / Redi=F): run the REAL
# FESOM2 multi-step lifecycle with Gent-McWilliams bolus advection ON (work_core GM config,
# Redi kept off) and FESOM3's lifecycle driver with FESOM3_FER_GM=1, compare the 13 NODE
# substeps x 5 probes over N steps. Isolates the GM tracer-transport loop end-to-end:
# producers (sw_alpha_beta -> sigma_xy) -> init_Redi_GM -> fer_solve_Gamma -> fer_gamma2vel
# -> vert_vel_ale fer_w -> solve_tracers_ale bolus add/advect/subtract, on the do_ic3d phc3.0
# state (real density gradients -> non-vacuous bolus). The GM-off gate is run_lifecycle_gate_core2.sh.
#
#   tools/run_lifecycle_gm_gate_core2.sh [run_dir] [nsteps]
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
RUN="${1:-/scratch/a/a270088/lifecycle_gm_core2}"
NSTEPS="${2:-3}"

echo "[1/3] FESOM2 oracle lifecycle ($NSTEPS steps, 1-rank CORE2, unforced, Fer_GM=ON)"
FER_GM=1 bash "$F3/tools/run_lifecycle_core2.sh" "$RUN" "$RUN/life_f2" "$NSTEPS" | tail -3

echo "[2/3] FESOM3 lifecycle ($NSTEPS steps, 1-rank CORE2, FESOM3_FER_GM=1)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE"
export FESOM_DUMP_FILE="$RUN/life_f3" FESOM_DUMP_MAXSTEPS="$NSTEPS" FESOM3_NSTEPS="$NSTEPS"
export FESOM3_FER_GM=1
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_lifecycle" 2>&1 | grep -E 'nod2D|IC T|Fer_GM|step|done'

echo "[3/3] compare (13 NODE substeps x 5 probes x $NSTEPS steps; SW_AB substep 2 ignored)"
python3 "$F3/tools/dump_diff.py" "$RUN/life_f2.00000" "$RUN/life_f3.00000" --ignore-substep=2
