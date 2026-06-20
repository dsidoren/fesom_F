#!/usr/bin/env bash
# Run the instrumented FESOM2 (build/bin/fesom.x) on pi at SINGLE RANK to drive the
# REAL ocean timestep (oce_timestep_ale) on prescribed state for the FESOM3 M2.9b
# step-ASSEMBLY byte-gate. step_dump_write fires at the END of ocean_setup (env
# FESOM_STEP_DUMP, npes==1) — after init_ale / init_thickness_ale / arrays_init /
# init_stiff_mat_ale, BEFORE forcing init (which hangs on a 1-rank login-node run).
# It prescribes the clean state + forcing arrays + reduced-M2 config, calls
# compute_vel_nodes + oce_timestep_ale ONCE (whose built-in dump_shim_record_node
# emit the per-substep NODE dumps to FESOM_DUMP_FILE), and stops.
# Needs the hand-crafted dist_1 in the pi mesh dir (shared with the other gates).
#
#   tools/run_stepdump_pi.sh [run_dir] [dump_prefix]
set -euo pipefail

F2=/home/a/a270088/port2/fesom2
RUN="${1:-/scratch/a/a270088/stepdump_pi}"
DUMP="${2:-$RUN/step_f2}"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1

rm -rf "$RUN"; mkdir -p "$RUN"
cp "$F2"/work_pi/namelist.* "$RUN"/
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x     # fresh build (has the step shim)
printf '0 1 1948\n0 1 1948\n' > "$RUN"/fesom.clock

python3 - "$RUN" "$F2" <<'PY'
import re, sys
run, f2 = sys.argv[1], sys.argv[2]
p = run + '/namelist.config'
s = open(p).read()
s = re.sub(r"ResultPath\s*=\s*'[^']*'",  "ResultPath       = './'", s, count=1)
s = re.sub(r"ClimateDataPath\s*=\s*'[^']*'",
           f"ClimateDataPath  = '{f2}/test/input/global/'", s, count=1)
open(p, 'w').write(s)
PY

cd "$RUN"
export FESOM_STEP_DUMP=1                  # trigger step_dump_write
export FESOM_DUMP_FILE="$DUMP"            # oce_timestep_ale built-in node dumps
export FESOM_DUMP_MAXSTEPS=1
ulimit -s unlimited
echo "run_stepdump_pi: 1 rank -> ${DUMP}.<rank>"
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 ./fesom.x > "$RUN/run.log" 2>&1 || true
if ls "${DUMP}".* >/dev/null 2>&1; then
    echo "run_stepdump_pi: done -> ${DUMP}.* ($(stat -c%s "${DUMP}".00000) bytes)"
else
    echo "run_stepdump_pi: NO DUMP — see $RUN/run.log"; tail -30 "$RUN/run.log"; exit 1
fi
