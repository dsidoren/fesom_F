#!/usr/bin/env bash
# M2.11c CORE2 dynamics-kernel byte-gate (the M2.1-M2.6 per-kernel pressure gate, on the
# 40x CORE2 mesh): dump FESOM2's REAL pressure_bv/pgf/vel_rhs/viscosity/ivertvisc/SSH-stiff/
# CG fields (1-rank CORE2, reduced-M2, prescribed analytic state) + FESOM3's transcribed
# kernels (1-rank CORE2, FESOM3_STEP_PER_DAY=48 so the SSH stiffness dt=1800 matches), and
# compare. Localizes the M2.11c-1 lifecycle CG 1-ULP divergence (ssh_stiff_diag/ssh_Aeta vs
# d_eta).
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
RUN="${1:-/scratch/a/a270088/pressuredump_core2}"

echo "[1/3] FESOM2 oracle pressure dump (1-rank CORE2, reduced-M2)"
bash "$F3/tools/run_pressuredump_core2.sh" "$RUN" "$RUN/pressure_f2.bin" >/dev/null

echo "[2/3] FESOM3 pressure dump (1-rank CORE2, step_per_day=48)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH"
export FESOM3_STEP_PER_DAY=48
export FESOM3_PRESSURE_OUT="$RUN/pressure_f3.bin"
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_pressuredump" >/dev/null 2>&1

echo "[3/3] compare"
python3 "$F3/tools/pressure_diff.py" "$RUN/pressure_f2.bin" "$RUN/pressure_f3.bin"
