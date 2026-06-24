#!/usr/bin/env bash
# M4d GM+Redi lifecycle byte-gate (CORE2 1-rank, unforced, Fer_GM=T + Redi=T) — the production
# tracer physics. The oracle keeps work_core Fer_GM/Redi=.true. (FER_GM=1 REDI=1); FESOM3
# enables both (FESOM3_FER_GM=1 FESOM3_REDI=1). Gates the 13 NODE substeps x 5 probes x N steps.
# Isolates the Redi terms (K13/K23 diff_part_hor_redi + K31/K32 diff_ver_part_redi_expl + K33
# diff_ver_part_impl_ale) on top of the proven M4c GM bolus: with the GM-only gate already green
# (run_lifecycle_gm_gate_core2.sh), a pass here pins the Redi isopycnal diffusion. NB a Redi-only
# (Fer_GM=F) gate is vacuous — work_core K_hor=0 makes Ki come solely from the Redi.and.Fer_GM
# coupling (Ki=max(fer_scal*Redi_Kmax,K_GM_min)).
#   tools/run_lifecycle_redi_gate_core2.sh [run_dir] [nsteps]
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
RUN="${1:-/scratch/a/a270088/lifecycle_redi_core2}"
NSTEPS="${2:-3}"

echo "[1/3] FESOM2 oracle lifecycle ($NSTEPS steps, 1-rank CORE2, unforced, Fer_GM+Redi=ON)"
FER_GM=1 REDI=1 bash "$F3/tools/run_lifecycle_core2.sh" "$RUN" "$RUN/life_f2" "$NSTEPS" | tail -3

echo "[2/3] FESOM3 lifecycle ($NSTEPS steps, 1-rank CORE2, FESOM3_FER_GM=1 FESOM3_REDI=1)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE"
export FESOM_DUMP_FILE="$RUN/life_f3" FESOM_DUMP_MAXSTEPS="$NSTEPS" FESOM3_NSTEPS="$NSTEPS"
export FESOM3_FER_GM=1 FESOM3_REDI=1
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_lifecycle" 2>&1 | grep -E 'nod2D|IC T|Redi|step|done'

echo "[3/3] compare (13 NODE substeps x 5 probes x $NSTEPS steps; SW_AB substep 2 ignored)"
python3 "$F3/tools/dump_diff.py" "$RUN/life_f2.00000" "$RUN/life_f3.00000" --ignore-substep=2
