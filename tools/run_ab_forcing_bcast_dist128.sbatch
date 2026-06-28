#!/bin/bash
# A/B: forcing rank-0-read+Bcast fix vs pre-fix all-ranks-read, dist_128 on `interactive`
# (1 node = 128 ranks; the only multi-rank partition available during the Sunday cooling
# maintenance). Two-point [100,600] -> steady-state ms/step = (mean600*600 - mean100*100)/500,
# parsed from the per-component TIMING lines. CAVEAT: on a SINGLE node the OS page cache masks
# part of the win (data reads warm after rank 0); the cross-node Lustre open/metadata contention
# that dominates at dist_512 is NOT exercised here, so the dist_128 delta is a LOWER BOUND on the
# real win. Both binaries are pre-built and pinned (deterministic build):
#   OLD = /scratch/a/a270088/ab_bcast/fesom_old   (all-ranks nc_open/read/close)
#   NEW = /scratch/a/a270088/ab_bcast/fesom_new   (if(mype==0) read + MPI_Bcast)
#SBATCH -p interactive
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=00:20:00
#SBATCH -A ab0995
set -uo pipefail
F3=/home/a/a270088/fesom3
SAVE=/scratch/a/a270088/ab_bcast
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
for TAG in old new; do
    BIN="$SAVE/fesom_$TAG"
    for NS in 100 600; do
        export FESOM3_NSTEPS=$NS
        RUN="$SAVE/ab_${TAG}_n${NS}"; mkdir -p "$RUN"; cd "$RUN"
        echo "################ $TAG  NSTEPS=$NS ################"
        srun -l -n 128 "$BIN" > "$RUN/log" 2>&1
        echo "rc=$?"
        grep -aE "TIMING|forcing|sbc / crossing|bulk NCAR|LOOP TOTAL" "$RUN/log" | sed 's/^[ ]*[0-9]*:[ ]*//'
    done
done
echo "################ steady-state = (mean600*600 - mean100*100)/500 per component; compare old vs new ################"
date
