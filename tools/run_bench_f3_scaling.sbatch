#!/bin/bash
# Clean FESOM3 speed/scaling benchmark — NO dump, production config (JRA55 1958, zstar+TKE+GM+Redi+sw_pene).
# Two-point timing per scale: run NSTEPS=100 and NSTEPS=600 from cold start; per-step = (T600-T100)/500
# cleanly removes the fixed setup cost (mesh read + IC + partition). The partition is auto-selected as
# dist_<SLURM_NTASKS>, so submit with --ntasks=<N> --nodes=<M> for each scaling point:
#   for d in 128 256 512 864; do
#     n=$(( (d+127)/128 ))   # nodes (128 cores/node)
#     sbatch --ntasks=$d --nodes=$n --job-name=benchf3_$d \
#            -o /scratch/a/a270088/benchf3_$d.%j.out tools/run_bench_f3_scaling.sbatch
#   done
#SBATCH -p compute
#SBATCH --time=00:25:00
#SBATCH -A ab0995
set -uo pipefail
F3=/home/a/a270088/fesom3
BIN="$F3/build_intel_dp/bin/fesom_lifecycle_native_mr"
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
NT="${SLURM_NTASKS:?must run under SLURM with --ntasks set}"

source "$F3/env.sh" intel >/dev/null 2>&1          # KNEM fix (L35)
ulimit -s 204800 2>/dev/null || true

# Full production config; NO FESOM_DUMP_FILE => the gid dumps are no-ops (clean run).
export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE"
export FESOM3_FORCING_DIR="$POOL"  FESOM3_FORCING=JRA55  FESOM3_START_CLOCK="0 1 1958"
export FESOM3_RUNOFF_FILE="$POOL/CORE2_runoff.nc" FESOM3_SSS_FILE="$POOL/PHC2_salx.nc" FESOM3_CHL_SWEENEY=1
export FESOM3_WHICH_ALE=zstar FESOM3_FER_GM=1 FESOM3_REDI=1 FESOM3_MIX_TKE=1 FESOM3_SW_PENE=1 FESOM3_WHICHEVP=0

date
echo "################################################################"
echo "## F3 SCALING BENCH: dist_$NT ($(( (NT+127)/128 )) nodes target), NO dump, two-point [100,600] steps"
echo "################################################################"
declare -A WALL
for NS in 100 600; do
    export FESOM3_NSTEPS=$NS
    RUN="/scratch/a/a270088/benchf3_d${NT}_n${NS}"; mkdir -p "$RUN"
    t0=$(date +%s.%N)
    srun -l -n "$NT" "$BIN" > "$RUN/log" 2>&1
    rc=$?
    t1=$(date +%s.%N)
    w=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f", b-a}')
    WALL[$NS]=$w
    done_line=$(grep -aE "done \($NS steps" "$RUN/log" | tail -1)
    echo "## dist_$NT NSTEPS=$NS  wall=${w}s  rc=$rc  [$done_line]"
done

t100="${WALL[100]}"; t600="${WALL[600]}"
awk -v t1="$t100" -v t6="$t600" -v nt="$NT" 'BEGIN{
    perstep=(t6-t1)/500.0;
    setup = t1 - 100*perstep;
    sps = (perstep>0)? 1.0/perstep : 0;
    printf "## ===== RESULT dist_%d =====\n", nt;
    printf "##   per-step (setup-removed) = %.2f ms\n", perstep*1000.0;
    printf "##   throughput              = %.2f steps/s  (= %.1f model-days/wall-hour @ 48 steps/day)\n", sps, sps*3600.0/48.0;
    printf "##   implied setup cost       = %.1f s\n", setup;
}'
date
