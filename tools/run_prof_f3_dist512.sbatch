#!/bin/bash
# F3 per-component PROFILING at dist_512 (matches the F2 internal-timer breakdown scale).
# Two-point [100,600] with FESOM3_TIMING=1 so steady-state per-component = (tot600 - tot100)/500
# removes the warm-up (tot = mean_ms/step * nsteps, parsed from the TIMING: lines).
#SBATCH -p compute
#SBATCH --nodes=4
#SBATCH --ntasks=512
#SBATCH --time=00:20:00
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
export FESOM3_TIMING=1   # emit the per-component breakdown

date
for NS in 100 600; do
    export FESOM3_NSTEPS=$NS
    RUN="/scratch/a/a270088/proff3_d512_n${NS}"; mkdir -p "$RUN"
    echo "################ dist_512 NSTEPS=$NS ################"
    srun -l -n 512 "$BIN" > "$RUN/log" 2>&1
    echo "rc=$?"
    grep -aE "TIMING" "$RUN/log" | sed 's/^[ ]*[0-9]*:[ ]*//'
done
echo "################ steady-state per-component (tot600-tot100)/500 will be computed from TIMING lines ################"
date
