#!/bin/bash
# Clean FESOM2 (oracle) speed benchmark — MIRRORS tools/run_bench_f3_scaling.sbatch for an
# apples-to-apples compute-vs-compute comparison. Same mesh, same production config (JRA55 1958,
# zstar + cvmix_TKE + GM + Redi + sw_pene), same two-point method: NSTEPS=100 and 600 from cold
# start, per-step = (T600-T100)/500 (removes the fixed setup + the one-time IC write at step 0).
#
# CLEAN because: (a) NO dump shim — FESOM_DUMP_FILE is left UNSET, so the gid dump is a no-op,
# matching the F3 bench; (b) FESOM2 output is monthly means (namelist.io io_list freq '1,m'),
# and 600 steps = 12.5 days crosses NO month boundary, so zero NetCDF write fires in either point.
# => the timed delta is pure dynamics, directly comparable to F3.
#
# Partition auto-selected as dist_<SLURM_NTASKS> (FESOM2 reads dist_<npes>). Submit per scale:
#   for d in 512 864; do n=$(( (d+127)/128 ));
#     sbatch --ntasks=$d --nodes=$n --job-name=benchf2_$d \
#            -o /scratch/a/a270088/benchf2_$d.%j.out tools/run_bench_f2_scaling.sbatch; done
#SBATCH -p compute
#SBATCH --time=00:25:00
#SBATCH -A ab0995
set -uo pipefail
F2=/home/a/a270088/port2/fesom2
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
NT="${SLURM_NTASKS:?must run under SLURM with --ntasks set}"
START_CLOCK="0 1 1958"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1   # KNEM fix
ulimit -s 204800 2>/dev/null || true

setup_run () {   # $1 = run dir, $2 = nsteps ; production JRA55 config, NO dump
    local RUN="$1" NS="$2"
    rm -rf "$RUN"; mkdir -p "$RUN"
    cp "$F2"/work_core/namelist.* "$RUN"/
    cp "$F2"/work_core/namelist.forcing "$RUN"/namelist.forcing   # JRA55-do default
    cp "$F2"/work_tke_dump/namelist.cvmix "$RUN"/                 # &param_tke (tke_cd=3.75)
    sed -i "s/^whichEVP *=.*/whichEVP = 0/" "$RUN"/namelist.ice
    ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x
    printf '%s\n%s\n' "$START_CLOCK" "$START_CLOCK" > "$RUN"/fesom.clock
    python3 - "$RUN" "$POOL" "$NS" <<'PY'
import re,sys
run,pool,nsteps=sys.argv[1:4]
p=run+'/namelist.config'; s=open(p).read()
s=re.sub(r"ResultPath\s*=\s*'[^']*'","ResultPath       = './'",s,1)
s=re.sub(r"which_ALE\s*=\s*'zlevel'","which_ALE          = 'zstar'",s,1)   # production zstar
s=re.sub(r"run_length\s*=\s*\d+",f"run_length        = {nsteps}",s,1)
open(p,'w').write(s)
p=run+'/namelist.oce'; s=open(p).read()
s=re.sub(r"mix_scheme\s*=\s*'KPP'","mix_scheme         = 'cvmix_TKE'",s,1)  # production TKE
# Fer_GM/Redi already .true. in work_core; use_sw_pene already .true. — production keeps them.
open(p,'w').write(s)
p=run+'/namelist.forcing'; s=open(p).read()
s=re.sub(r"nm_runoff_file\s*=\s*'[^']*'",   f"nm_runoff_file ='{pool}/CORE2_runoff.nc'",s,1)
s=re.sub(r"nm_sss_data_file\s*=\s*'[^']*'", f"nm_sss_data_file ='{pool}/PHC2_salx.nc'",s,1)
open(p,'w').write(s)   # chl='Sweeney' kept (production)
PY
}

date
echo "################################################################"
echo "## F2 SCALING BENCH: dist_$NT ($(( (NT+127)/128 )) nodes), NO dump, monthly-output (no write <1mo), two-point [100,600]"
echo "################################################################"
declare -A WALL
for NS in 100 600; do
    RUN="/scratch/a/a270088/benchf2_d${NT}_n${NS}"
    setup_run "$RUN" "$NS"
    cd "$RUN"
    # NOTE: FESOM_DUMP_FILE intentionally UNSET => dump shim no-op (clean, matches F3 bench).
    t0=$(date +%s.%N)
    srun -l -n "$NT" ./fesom.x > "$RUN/run.log" 2>&1
    rc=$?
    t1=$(date +%s.%N)
    w=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f", b-a}')
    WALL[$NS]=$w
    last=$(grep -aiE "FESOM Run is finished|fesom.*finished|====" "$RUN/run.log" | tail -1)
    echo "## dist_$NT NSTEPS=$NS  wall=${w}s  rc=$rc  [$last]"
done

t100="${WALL[100]}"; t600="${WALL[600]}"
awk -v t1="$t100" -v t6="$t600" -v nt="$NT" 'BEGIN{
    perstep=(t6-t1)/500.0;
    setup = t1 - 100*perstep;
    sps = (perstep>0)? 1.0/perstep : 0;
    printf "## ===== RESULT dist_%d (FESOM2 oracle) =====\n", nt;
    printf "##   per-step (setup-removed) = %.2f ms\n", perstep*1000.0;
    printf "##   throughput              = %.2f steps/s  (= %.1f model-days/wall-hour @ 48 steps/day)\n", sps, sps*3600.0/48.0;
    printf "##   implied setup cost       = %.1f s\n", setup;
}'
date
