#!/usr/bin/env bash
# M5b KPP-mixing lifecycle byte-gate (CORE2 1-rank, UNFORCED, mix_scheme='KPP'): run the REAL
# FESOM2 multi-step lifecycle with the production K-Profile Parameterization ON (work_core KPP
# config; use_sw_pene off since unforced => sw_3d=0, ghats unused in the TDMA) and FESOM3's
# lifecycle driver with FESOM3_MIX_KPP=1, compare the 13 NODE substeps x 5 probes over N steps.
# Exercises the KPP Kv/Av integration end-to-end: oce_mixing_kpp_driver (prestep dVsq/ustar/Bo
# + dbsfc(pressure_bv) -> ri_iwmix -> bldepth -> blmix_kpp -> enhance -> combine -> node->elem
# viscAE average) -> Kv=Kv_double(:,:,1) -> mo_convect; Av stays element-based so
# impl_vert_visc_ale is UNCHANGED. Runs on the do_ic3d phc3.0 IC (real density structure ->
# non-vacuous boundary layer). Optional FER_GM=1/REDI=1 layer GM+Redi on top (production combo).
# The PP (no-KPP) gate is run_lifecycle_gate_core2.sh.
#
#   tools/run_lifecycle_kpp_gate_core2.sh [run_dir] [nsteps]      env: FER_GM=0/1 REDI=0/1
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
RUN="${1:-/scratch/a/a270088/lifecycle_kpp_core2}"
NSTEPS="${2:-3}"
FER_GM="${FER_GM:-0}"
REDI="${REDI:-0}"

echo "[1/3] FESOM2 oracle lifecycle ($NSTEPS steps, 1-rank CORE2, unforced, KPP=ON FER_GM=$FER_GM REDI=$REDI)"
MIX_KPP=1 FER_GM="$FER_GM" REDI="$REDI" \
    bash "$F3/tools/run_lifecycle_core2.sh" "$RUN" "$RUN/life_f2" "$NSTEPS" | tail -3

echo "[2/3] FESOM3 lifecycle ($NSTEPS steps, 1-rank CORE2, FESOM3_MIX_KPP=1 FER_GM=$FER_GM REDI=$REDI)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE"
export FESOM_DUMP_FILE="$RUN/life_f3" FESOM_DUMP_MAXSTEPS="$NSTEPS" FESOM3_NSTEPS="$NSTEPS"
export FESOM3_MIX_KPP=1
if [ "$FER_GM" = "1" ]; then export FESOM3_FER_GM=1; fi
if [ "$REDI"   = "1" ]; then export FESOM3_REDI=1;   fi
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_lifecycle" 2>&1 | grep -E 'nod2D|IC T|KPP|Fer_GM|step|done'

echo "[3/3] compare (13 NODE substeps x 5 probes x $NSTEPS steps; SW_AB substep 2 ignored)"
python3 "$F3/tools/dump_diff.py" "$RUN/life_f2.00000" "$RUN/life_f3.00000" --ignore-substep=2
