#!/usr/bin/env bash
# Run the instrumented FESOM2 (fresh build, build/bin/fesom.x) on pi at SINGLE RANK
# to dump the surface-forcing atmospheric fields for the FESOM3 M2.10a forcing-read
# byte-gate. forcing_dump_write fires right AFTER forcing_setup (env FESOM_FORCING_DUMP,
# npes==1) — the REAL sbc_ini already ran (netCDF axes read + cold-start getcoeffld +
# g2r). The shim pins the model time, drives the REAL sbc_do, maps atmdata -> the
# physical node arrays, dumps, and stops.
#
# REQUIRES the 1-rank next_io_rank fix in io_netcdf_workaround_module.F90 (the async
# IO-rank selector infinite-recurses at npes==1 otherwise — the original L8 hang).
# Reads the CORE2 forcing stubs in $F2/test/input/global/ (u_10.1948.nc etc.).
#
#   tools/run_forcingdump_pi.sh [run_dir] [out_file]
set -euo pipefail

F2=/home/a/a270088/port2/fesom2
RUN="${1:-/scratch/a/a270088/forcingdump_pi}"
OUT="${2:-$RUN/forcing_f2.bin}"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1

rm -rf "$RUN"; mkdir -p "$RUN"
cp "$F2"/work_pi/namelist.* "$RUN"/
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x     # fresh build (has forcing dump + np=1 fix)
printf '0 1 1948\n0 1 1948\n' > "$RUN"/fesom.clock

python3 - "$RUN" "$F2" <<'PY'
import re, sys
run, f2 = sys.argv[1], sys.argv[2]
p = run + '/namelist.config'
s = open(p).read()
s = re.sub(r'run_length\s*=\s*\S+',      'run_length        = 1',   s, count=1)
s = re.sub(r"run_length_unit\s*=\s*\S+", "run_length_unit   = 's'", s, count=1)
s = re.sub(r"ResultPath\s*=\s*'[^']*'",  "ResultPath       = './'", s, count=1)
s = re.sub(r"ClimateDataPath\s*=\s*'[^']*'",
           f"ClimateDataPath  = '{f2}/test/input/global/'", s, count=1)
open(p, 'w').write(s)
PY

cd "$RUN"
export FESOM_FORCING_DUMP="$OUT"
ulimit -s unlimited
echo "run_forcingdump_pi: 1 rank -> $OUT"
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 ./fesom.x > "$RUN/run.log" 2>&1 || true
if [[ -f "$OUT" ]]; then
    echo "run_forcingdump_pi: done -> $OUT ($(stat -c%s "$OUT") bytes)"
else
    echo "run_forcingdump_pi: NO DUMP — see $RUN/run.log"; tail -25 "$RUN/run.log"; exit 1
fi
