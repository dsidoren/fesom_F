#!/bin/bash
# F3 per-component PROFILING at dist_128 (1 node — runs on the `interactive` partition during the
# compute-partition cooling maintenance). The F3-vs-F2 residual is compute-bound, so the single-node
# per-component split localizes it. Two-point [100,600] with FESOM3_TIMING=1:
#   steady-state ms/step per component = (mean600*600 - mean100*100)/500   (parsed from TIMING: lines)
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=00:15:00
#SBATCH -A ab0995
set -uo pipefail
F3=/home/a/a270088/fesom3
BIN="$F3/build_intel_dp/bin/fesom_lifecycle_native_mr"
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
source "$F3/env.sh" intel >/dev/null 2>&1
ulimit -s 204800 2>/dev/null || true
export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE"
export FESOM3_FORCING_DIR="$POOL"  FESOM3_FORCING=JRA55  FESOM3_START_CLOCK="0 1 1958"
export FESOM3_RUNOFF_FILE="$POOL/CORE2_runoff.nc" FESOM3_SSS_FILE="$POOL/PHC2_salx.nc" FESOM3_CHL_SWEENEY=1
export FESOM3_WHICH_ALE=zstar FESOM3_FER_GM=1 FESOM3_REDI=1 FESOM3_MIX_TKE=1 FESOM3_SW_PENE=1 FESOM3_WHICHEVP=0
export FESOM3_TIMING=1

date
for NS in 100 600; do
    export FESOM3_NSTEPS=$NS
    RUN="/scratch/a/a270088/proff3_d128_n${NS}"; mkdir -p "$RUN"
    echo "################ dist_128 NSTEPS=$NS ################"
    srun -l -n 128 "$BIN" > "$RUN/log" 2>&1
    echo "rc=$?"
    grep -aE "TIMING" "$RUN/log" | sed 's/^[ ]*[0-9]*:[ ]*//'
done
date
