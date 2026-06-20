#!/usr/bin/env bash
# Run the instrumented FESOM2 (build/bin/fesom.x -> libfesom.so) on CORE2 at SINGLE
# RANK to dump the fully-built mesh geometry for the FESOM3 M2.11a geometry byte-gate.
# geom_dump_write fires at the end of mesh_setup (env FESOM_GEOM_DUMP, npes==1) and
# MPI_FINALIZE+stops right there (oce_mesh.F90:206 -> fesom_geom_dump.F90:88), BEFORE
# ocean_setup/forcing — so this is fast and reads only the mesh (no IC/forcing).
# Needs the hand-crafted dist_1 in the CORE2 mesh dir (tools/make_dist1.py made it).
#
#   tools/run_geomdump_core2.sh [run_dir] [out_file]
set -euo pipefail

F2=/home/a/a270088/port2/fesom2
RUN="${1:-/scratch/a/a270088/geomdump_core2}"
OUT="${2:-$RUN/geom_f2.bin}"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1

rm -rf "$RUN"; mkdir -p "$RUN"
cp "$F2"/work_core/namelist.* "$RUN"/
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x     # loads build/lib64/libfesom.so (has geom dump)
printf '0 1 1958\n0 1 1958\n' > "$RUN"/fesom.clock   # work_core yearnew=1958; geom dump stops before forcing

python3 - "$RUN" <<'PY'
import re, sys
run = sys.argv[1]
p = run + '/namelist.config'
s = open(p).read()
s = re.sub(r"ResultPath\s*=\s*'[^']*'", "ResultPath       = './'", s, count=1)
open(p, 'w').write(s)
PY

cd "$RUN"
export FESOM_GEOM_DUMP="$OUT"
ulimit -s unlimited
echo "run_geomdump_core2: 1 rank -> $OUT"
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 ./fesom.x > "$RUN/run.log" 2>&1 || true
if [[ -f "$OUT" ]]; then
    echo "run_geomdump_core2: done -> $OUT ($(stat -c%s "$OUT") bytes)"
else
    echo "run_geomdump_core2: NO DUMP — see $RUN/run.log"; tail -25 "$RUN/run.log"; exit 1
fi
