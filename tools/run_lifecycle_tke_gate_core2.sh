#!/usr/bin/env bash
# M7b TKE-mixing lifecycle byte-gate (CORE2 1-rank, UNFORCED, mix_scheme='cvmix_TKE').
# Tests the prognostic tke recurrence + the Av/Kv producer integration + mo_convect-after-TKE.
# Unforced ⇒ stress_node_surf=0 ⇒ forc_tke_surf=0 ⇒ the surface wind term is NOT yet exercised
# (that is M7c); this isolates shear/buoyancy production + dissipation + the step->step recurrence.
#
#   tools/run_lifecycle_tke_gate_core2.sh [run_dir] [nsteps]
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
RUN="${1:-/scratch/a/a270088/lifecycle_tke_core2}"
NSTEPS="${2:-3}"

echo "[1/3] FESOM2 oracle lifecycle ($NSTEPS steps, 1-rank CORE2, unforced, mix_scheme=cvmix_TKE)"
MIX_TKE=1 bash "$F3/tools/run_lifecycle_core2.sh" "$RUN" "$RUN/life_f2" "$NSTEPS" | tail -3

echo "[2/3] FESOM3 lifecycle ($NSTEPS steps, 1-rank CORE2, FESOM3_MIX_TKE=1)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE"
export FESOM_DUMP_FILE="$RUN/life_f3" FESOM_DUMP_MAXSTEPS="$NSTEPS" FESOM3_NSTEPS="$NSTEPS"
export FESOM3_MIX_TKE=1
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_lifecycle" 2>&1 | grep -E 'nod2D|IC T|TKE|step|done'

echo "[3/3] compare (13 NODE substeps x 5 probes x $NSTEPS steps; SW_AB substep 2 ignored)"
python3 "$F3/tools/dump_diff.py" "$RUN/life_f2.00000" "$RUN/life_f3.00000" --ignore-substep=2
