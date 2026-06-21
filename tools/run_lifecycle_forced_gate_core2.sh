#!/usr/bin/env bash
# M2.11c-2 FORCED multi-step lifecycle byte-gate (CORE2 1-rank): the REAL forced FESOM2
# lifecycle (use_ice=.true., real CORE2 forcing + ice EVP + oce_fluxes) vs FESOM3's
# lifecycle driver PRESCRIBING the oracle's per-step air-sea fluxes (the M3 gap). Compares
# the 13 NODE substeps x 5 probes over N steps (--ignore-substep=2).
#
# Like the unforced gate (run_lifecycle_gate_core2.sh + LESSONS L28): every substep that
# does NOT depend on the free-surface CG is byte-identical (max|delta|=0) — now with REAL
# wind stress + heat/freshwater/salt fluxes prescribed from the oracle. The CG d_eta and
# its downstream are at the iterative-solver reproducibility floor (see L28). The companion
# run_pressure_gate_core2.sh proves the whole dynamical core max|delta|=0 on CORE2.
#
#   tools/run_lifecycle_forced_gate_core2.sh [run_dir] [nsteps]
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
RUN="${1:-/scratch/a/a270088/lifecycle_forced_core2}"
NSTEPS="${2:-3}"

echo "[1/3] FESOM2 oracle FORCED lifecycle ($NSTEPS steps, 1-rank CORE2, use_ice)"
bash "$F3/tools/run_lifecycle_forced_core2.sh" "$RUN" "$RUN/lifef_f2" "$RUN/flux_f2" "$NSTEPS" | tail -4

echo "[2/3] FESOM3 lifecycle ($NSTEPS steps) PRESCRIBING the oracle's per-step fluxes"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE"
export FESOM3_FLUX_FILE="$RUN/flux_f2.00000"
export FESOM_DUMP_FILE="$RUN/lifef_f3" FESOM_DUMP_MAXSTEPS="$NSTEPS" FESOM3_NSTEPS="$NSTEPS"
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_lifecycle" 2>&1 | grep -E 'nod2D|FORCED|step|done'

echo "[3/3] compare (13 NODE substeps x 5 probes x $NSTEPS steps; SW_AB substep 2 ignored)"
python3 "$F3/tools/dump_diff.py" "$RUN/lifef_f2.00000" "$RUN/lifef_f3.00000" --ignore-substep=2
