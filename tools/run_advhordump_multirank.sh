#!/usr/bin/env bash
# M2.12b: run the instrumented FESOM2 on pi at NP ranks (NP>1, dist_<NP>) to dump the
# per-rank LOCAL OWNED tracer-advection fields for the multi-rank advection byte-gate.
# advhor_dump_write fires at the END of ocean_setup (env FESOM_ADVHOR_DUMP, all ranks,
# after init_thickness_ale + muscl_adv_init) and STOPS before forcing — fast. Each rank
# writes <out>.<mype5> (1..myDim_* owned slices). Needs dist_<NP>/ in the pi mesh dir.
#
#   tools/run_advhordump_multirank.sh [run_dir] [out_prefix] [np]
set -euo pipefail

F2=/home/a/a270088/port2/fesom2
RUN="${1:-/scratch/a/a270088/advhordump_mr_pi}"
OUT="${2:-$RUN/advhor_f2.bin}"
NP="${3:-2}"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1

rm -rf "$RUN"; mkdir -p "$RUN"
cp "$F2"/work_pi/namelist.* "$RUN"/
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x     # fresh build (has multi-rank advhor dump)
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
export FESOM_ADVHOR_DUMP="$OUT"
ulimit -s unlimited
echo "run_advhordump_multirank: $NP ranks (dist_$NP) -> ${OUT}.<rank>"
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n "$NP" ./fesom.x > "$RUN/run.log" 2>&1 || true
if ls "${OUT}".* >/dev/null 2>&1; then
    echo "run_advhordump_multirank: done -> $(ls "${OUT}".* | tr '\n' ' ')"
    for f in "${OUT}".*; do echo "    $(basename "$f"): $(stat -c%s "$f") bytes"; done
else
    echo "run_advhordump_multirank: NO DUMP — see $RUN/run.log"; tail -25 "$RUN/run.log"; exit 1
fi
