#!/bin/bash
# A/B of the forcing persistent-handle + double-buffer fix at dist_512 (4 nodes, `compute`
# partition during the compute cooling maintenance). Same binary, toggled by FESOM3_FORCING_PERSIST:
#   0 = legacy (nc_open + read BOTH brackets + nc_close every crossing, every rank)
#   1 = fix     (persistent handle, reopen only on year change; read only the 1 new slice/crossing)
# Two-point [100,600] per arm with FESOM3_TIMING=1 so steady-state per-component =
#   (mean600*600 - mean100*100)/500 removes the cold-start warm-up. The new TMR_FRC_OPEN /
# read+inflate / interp sub-timers under `sbc / crossing` show WHERE the crossing cost goes.
# NB single-node page-cache masks the absolute disk cost (see the investigation doc); the OPEN and
# READ(inflate) costs are CPU-bound so the A/B *mechanism* (open->~0, read halved) still shows here.
#SBATCH -p compute
#SBATCH --nodes=4
#SBATCH --ntasks=512
#SBATCH --time=00:25:00
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
for PERSIST in 0 1; do
    export FESOM3_FORCING_PERSIST=$PERSIST
    echo "######## ARM FESOM3_FORCING_PERSIST=$PERSIST ($([ $PERSIST = 0 ] && echo legacy || echo fix)) ########"
    for NS in 100 600; do
        export FESOM3_NSTEPS=$NS
        RUN="/scratch/a/a270088/ab_persist_d512_p${PERSIST}_n${NS}"; mkdir -p "$RUN"
        srun -l -n 512 "$BIN" > "$RUN/log" 2>&1
        echo "--- PERSIST=$PERSIST NSTEPS=$NS rc=$? ---"
        grep -aE "TIMING|forcing|sbc / crossing|open/reopen|read\+inflate|interp loop" "$RUN/log" | sed 's/^[ ]*[0-9]*:[ ]*//'
    done
done
date
echo "################ steady-state = (mean600*600 - mean100*100)/500 per row ################"
