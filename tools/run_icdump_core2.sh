#!/usr/bin/env bash
# Run the instrumented FESOM2 (build/bin/fesom.x -> libfesom.so) on CORE2 at SINGLE
# RANK to dump the live 3-D initial conditions (Tclim/Sclim) for the FESOM3 M2.11b
# do_ic3d byte-gate. ic_dump_write fires at the END of ocean_setup (env FESOM_IC_DUMP,
# npes==1), AFTER the REAL oce_initial_state::do_ic3d read phc3.0_winter.nc +
# interpolated + extrapolated + insitu2pot'd the tracers, then MPI_FINALIZE+stops
# BEFORE forcing init (the 1-rank login-node forcing hang, L8).
#
# Reduced-M2 namelist overrides on work_core: which_ALE 'zlevel'->'linfs' (so the
# oracle's init_ale Z_3d_n is the reference full-cell depth FESOM3 builds — same at
# init eta=0), mix_scheme 'KPP'->'PP', Fer_GM/Redi -> .false. (the proven M2 oracle
# config; KPP/GM/Redi are irrelevant to do_ic3d but this lightens ocean_setup). Needs
# the hand-crafted dist_1 in the CORE2 mesh dir (tools/make_dist1.py made it for M2.11a).
#
#   tools/run_icdump_core2.sh [run_dir] [out_file]
set -euo pipefail

F2=/home/a/a270088/port2/fesom2
RUN="${1:-/scratch/a/a270088/icdump_core2}"
OUT="${2:-$RUN/ic_f2.bin}"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1

rm -rf "$RUN"; mkdir -p "$RUN"
cp "$F2"/work_core/namelist.* "$RUN"/
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x     # loads build/lib64/libfesom.so (has IC dump shim)
printf '0 1 1958\n0 1 1958\n' > "$RUN"/fesom.clock   # work_core yearnew=1958; IC dump stops before forcing

python3 - "$RUN" <<'PY'
import re, sys
run = sys.argv[1]
# namelist.config: ResultPath -> ./, which_ALE zlevel -> linfs (reduced-M2; matches FESOM3 Z_3d_n)
p = run + '/namelist.config'
s = open(p).read()
s = re.sub(r"ResultPath\s*=\s*'[^']*'", "ResultPath       = './'", s, count=1)
s = re.sub(r"which_ALE\s*=\s*'zlevel'", "which_ALE          = 'linfs'", s, count=1)
open(p, 'w').write(s)
# namelist.oce: mix_scheme KPP -> PP, Fer_GM/Redi -> .false.
p = run + '/namelist.oce'
s = open(p).read()
s = re.sub(r"mix_scheme\s*=\s*'KPP'", "mix_scheme         = 'PP'", s, count=1)
s = re.sub(r"Fer_GM\s*=\s*\.true\.",  "Fer_GM             = .false.", s, count=1)
s = re.sub(r"Redi\s*=\s*\.true\.",    "Redi               = .false.", s, count=1)
open(p, 'w').write(s)
PY

cd "$RUN"
export FESOM_IC_DUMP="$OUT"
ulimit -s unlimited
echo "run_icdump_core2: 1 rank -> $OUT"
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 ./fesom.x > "$RUN/run.log" 2>&1 || true
if [[ -f "$OUT" ]]; then
    echo "run_icdump_core2: done -> $OUT ($(stat -c%s "$OUT") bytes)"
    grep -E 'gobal (max|min) init' "$RUN/run.log" || true
else
    echo "run_icdump_core2: NO DUMP — see $RUN/run.log"; tail -30 "$RUN/run.log"; exit 1
fi
