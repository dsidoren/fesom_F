#!/usr/bin/env bash
# M7a oracle: run the REAL FESOM2 (cvmix_TKE) work_tke_dump config on CORE2 1-rank with
# FESOM_TKE_DUMP_DIR set, so the tke_dump_mod instrumentation (in the prebuilt dump-capable
# libfesom.so) writes the per-step per-gid INPUT + OUTPUT dumps:
#   tke_dump_s<step>_<tag>_rank0.txt  (tag in normstress/vshear2/bvfreq2/dztrr/tkeold +
#                                      tke/tkeav/tkekv + 10 budget diagnostics + kv/av)
# These are the byte-exact reference for the M7a-1 controlled-replay gate (fesom_tkereplay)
# and the M7a-2 native-driver gate (fesom_tkedump).
#
# Config = work_tke_dump verbatim (cvmix_TKE + dt=1800 + linfs + use_ice + use_sw_pene, JRA
# forcing, cold start 1958), 3 steps. Only ResultPath is repointed to the run dir. The oracle
# binary build/bin/fesom.x links build/lib64/libfesom.so (has cvmix_TKE + tke_dump_mod).
#
#   tools/run_tkedump_oracle.sh [run_dir] [tke_dump_dir]
set -euo pipefail
F2=/home/a/a270088/port2/fesom2
RUN="${1:-/scratch/a/a270088/tkedump_oracle}"
TKEDIR="${2:-$RUN/tke_oracle}"
NSTEPS="${NSTEPS:-3}"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1
rm -rf "$RUN"; mkdir -p "$RUN" "$TKEDIR"
cp "$F2"/work_tke_dump/namelist.* "$RUN"/
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x
printf '0 1 1958\n0 1 1958\n' > "$RUN"/fesom.clock

python3 - "$RUN" "$NSTEPS" <<'PY'
import re,sys
run,nsteps=sys.argv[1:3]
p=run+'/namelist.config'; s=open(p).read()
s=re.sub(r"ResultPath\s*=\s*'[^']*'","ResultPath       = './'",s,1)
s=re.sub(r"run_length\s*=\s*\d+",f"run_length        = {nsteps}",s,1)
open(p,'w').write(s)
PY

cd "$RUN"
export FESOM_TKE_DUMP_DIR="$TKEDIR" FESOM_TKE_DUMP_STEPS="$NSTEPS"
ulimit -s unlimited
echo "run_tkedump_oracle: CORE2 1-rank, mix_scheme=cvmix_TKE, $NSTEPS steps -> $TKEDIR"
timeout 600 mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 ./fesom.x > "$RUN/run.log" 2>&1 || true
if ls "$TKEDIR"/tke_dump_s1_tke_rank0.txt >/dev/null 2>&1; then
    echo "run_tkedump_oracle: OK -> $(ls "$TKEDIR" | wc -l) dump files"
    grep -A12 -i "initialise CVMIX_TKE" "$RUN/run.log" | head -16
else
    echo "run_tkedump_oracle: MISSING tke dumps — see $RUN/run.log"; tail -40 "$RUN/run.log"; exit 1
fi
