#!/usr/bin/env bash
# M5a-1 oracle: run the REAL FESOM2 on CORE2 1-rank with mix_scheme='KPP' so the
# KPP scheme's oce_mixing_kpp_init fires at array_setup and the FESOM_KPP_DUMP_DIR
# instrumentation writes kpp_init_rank0.txt + kpp_wscale_rank0.txt (the rank-
# independent Vtc/cg/deltaz/deltau + wmt/wst lookup tables + the wscale sweep).
# The init dump is emitted BEFORE the first timestep, so it lands even if a later
# (un-ported KPP path) step would crash. Reuses the proven forced-CORE2 setup
# (run_lifecycle_forced_core2.sh path-fixing) but KEEPS KPP (no PP downgrade).
#
#   tools/run_kppinit_oracle.sh [run_dir] [kpp_dump_dir]
set -euo pipefail
F2=/home/a/a270088/port2/fesom2
STUB="$F2/test/input/global"
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
RUN="${1:-/scratch/a/a270088/kppinit_oracle}"
KPPDIR="${2:-$RUN/kpp_oracle}"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1
rm -rf "$RUN"; mkdir -p "$RUN" "$KPPDIR"
cp "$F2"/work_core/namelist.* "$RUN"/
cp "$F2"/work_core/namelist.forcing.CORE2 "$RUN"/namelist.forcing
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x
printf '0 1 1948\n0 1 1948\n' > "$RUN"/fesom.clock

python3 - "$RUN" "$STUB" "$POOL" <<'PY'
import re,sys
run,stub,pool=sys.argv[1:4]
p=run+'/namelist.config'; s=open(p).read()
s=re.sub(r"ResultPath\s*=\s*'[^']*'","ResultPath       = './'",s,1)
s=re.sub(r"which_ALE\s*=\s*'zlevel'","which_ALE          = 'linfs'",s,1)   # proven-runnable reduced config
s=re.sub(r"yearnew\s*=\s*1958","yearnew = 1948",s,1)
s=re.sub(r"run_length\s*=\s*\d+","run_length        = 1",s,1)              # 1 step is plenty; init dumps at setup
s=re.sub(r"include_fleapyear\s*=\s*\.true\.","include_fleapyear = .false.",s,1)  # CORE noleap
open(p,'w').write(s)
p=run+'/namelist.forcing'; s=open(p).read()
s=s.replace("FORCING/CORE2/", stub+"/")
s=re.sub(r"nm_runoff_file\s*=\s*'[^']*'",   f"nm_runoff_file ='{pool}/CORE2_runoff.nc'",s,1)
s=re.sub(r"nm_sss_data_file\s*=\s*'[^']*'", f"nm_sss_data_file ='{pool}/PHC2_salx.nc'",s,1)
open(p,'w').write(s)
PY

cd "$RUN"
export FESOM_KPP_DUMP_DIR="$KPPDIR" FESOM_KPP_DUMP_STEP=1
ulimit -s unlimited
echo "run_kppinit_oracle: CORE2 1-rank, mix_scheme=KPP -> $KPPDIR"
timeout 600 mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 ./fesom.x > "$RUN/run.log" 2>&1 || true
if ls "$KPPDIR"/kpp_init_rank0.txt >/dev/null 2>&1; then
    echo "run_kppinit_oracle: OK -> $(ls "$KPPDIR")"
    wc -l "$KPPDIR"/kpp_init_rank0.txt "$KPPDIR"/kpp_wscale_rank0.txt
else
    echo "run_kppinit_oracle: MISSING init dump — see $RUN/run.log"; tail -40 "$RUN/run.log"; exit 1
fi
