#!/usr/bin/env bash
# M2.12a: run the instrumented FESOM2 on pi at NP ranks (NP>1, dist_<NP>) to dump
# the per-rank LOCAL OWNED mesh geometry for the multi-rank geometry byte-gate.
# geom_dump_write fires at mesh_setup (env FESOM_GEOM_DUMP, all ranks) and stops
# BEFORE forcing init, so this is fast and avoids the 1-rank forcing hang. Each rank
# writes <out>.<mype5> (1..myDim_* owned slices). Needs dist_<NP>/ in the pi mesh dir.
#
#   tools/run_geomdump_multirank.sh [run_dir] [out_prefix] [np]
set -euo pipefail

F2=/home/a/a270088/port2/fesom2
RUN="${1:-/scratch/a/a270088/geomdump_mr_pi}"
OUT="${2:-$RUN/geom_f2.bin}"
NP="${3:-2}"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1

rm -rf "$RUN"; mkdir -p "$RUN"
cp "$F2"/work_pi/namelist.* "$RUN"/
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x     # fresh build (has multi-rank geom dump)
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
echo "run_geomdump_multirank: $NP ranks (dist_$NP) -> ${OUT}.<rank>"
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n "$NP" ./fesom.x > "$RUN/run.log" 2>&1 || true
if ls "${OUT}".* >/dev/null 2>&1; then
    echo "run_geomdump_multirank: done -> $(ls "${OUT}".* | tr '\n' ' ')"
    for f in "${OUT}".*; do echo "    $(basename "$f"): $(stat -c%s "$f") bytes"; done
else
    echo "run_geomdump_multirank: NO DUMP — see $RUN/run.log"; tail -20 "$RUN/run.log"; exit 1
fi
