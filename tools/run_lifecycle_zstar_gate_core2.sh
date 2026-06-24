#!/usr/bin/env bash
# M6a-2 zstar ALE byte-gate (CORE2 1-rank, UNFORCED): the same multi-step lifecycle as
# run_lifecycle_gate_core2.sh but with which_ALE='zstar' instead of 'linfs'. This exercises
# the production free-surface treatment — the Shchepetkin density-Jacobian PGF (M6a-1), the
# per-step SSH-stiffness 2nd-term update (update_stiff_mat_ale), and the vert_vel_ale /
# update_thickness_ale proportional-stretch thickness commit — on the reduced-M2 column
# physics (PP, no GM/Redi). UNFORCED => water_flux=0, so the freshwater path (M6a-3) is
# dormant and this isolates the ALE machinery. SSH evolves from the dynamics divergence, so
# hbar/hnode now EVOLVE (linfs kept them fixed) — a real test of the thickness commit.
#
#   tools/run_lifecycle_zstar_gate_core2.sh [run_dir] [nsteps]
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
RUN="${1:-/scratch/a/a270088/lifecycle_zstar_core2}"
NSTEPS="${2:-3}"

echo "[1/3] FESOM2 oracle lifecycle ($NSTEPS steps, 1-rank CORE2, unforced, which_ALE=zstar)"
WHICH_ALE=zstar bash "$F3/tools/run_lifecycle_core2.sh" "$RUN" "$RUN/life_f2" "$NSTEPS" | tail -3

echo "[2/3] FESOM3 lifecycle ($NSTEPS steps, 1-rank CORE2, FESOM3_WHICH_ALE=zstar)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE"
export FESOM_DUMP_FILE="$RUN/life_f3" FESOM_DUMP_MAXSTEPS="$NSTEPS" FESOM3_NSTEPS="$NSTEPS"
export FESOM3_WHICH_ALE=zstar
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_lifecycle" 2>&1 | grep -E 'nod2D|IC T|step|done'

echo "[3/3] compare (13 NODE substeps x 5 probes x $NSTEPS steps; SW_AB substep 2 ignored)"
python3 "$F3/tools/dump_diff.py" "$RUN/life_f2.00000" "$RUN/life_f3.00000" --ignore-substep=2
