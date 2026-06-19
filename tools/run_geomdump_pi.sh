#!/usr/bin/env bash
# Run the instrumented FESOM2 (fresh build, build/bin/fesom.x) on pi at SINGLE RANK
# to dump the fully-built mesh geometry for the FESOM3 geometry byte-gate
# (deferred M0.7 -> M1). geom_dump_write fires at mesh_setup (env FESOM_GEOM_DUMP,
# npes==1) and stops BEFORE forcing init, so this is fast and avoids the 1-rank
# forcing hang. Needs the hand-crafted dist_1 in the pi mesh dir (see tools/).
#
#   tools/run_geomdump_pi.sh [run_dir] [out_file]
set -euo pipefail

F2=/home/a/a270088/port2/fesom2
RUN="${1:-/scratch/a/a270088/geomdump_pi}"
OUT="${2:-$RUN/geom_f2.bin}"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1

rm -rf "$RUN"; mkdir -p "$RUN"
cp "$F2"/work_pi/namelist.* "$RUN"/
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x     # fresh build (has geom dump)
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
export FESOM_GEOM_DUMP="$OUT"
ulimit -s unlimited
echo "run_geomdump_pi: 1 rank -> $OUT"
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 ./fesom.x > "$RUN/run.log" 2>&1 || true
if [[ -f "$OUT" ]]; then
    echo "run_geomdump_pi: done -> $OUT ($(stat -c%s "$OUT") bytes)"
else
    echo "run_geomdump_pi: NO DUMP — see $RUN/run.log"; tail -15 "$RUN/run.log"; exit 1
fi
