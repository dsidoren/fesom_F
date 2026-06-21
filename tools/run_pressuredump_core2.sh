#!/usr/bin/env bash
# M2.11c diagnostic / CORE2 dynamics-kernel gate: run the instrumented FESOM2 pressure
# shim (fesom_pressure_dump.F90, env FESOM_PRESSURE_DUMP, end-of-ocean_setup, npes==1) on
# the CORE2 mesh, reduced-M2 (linfs/PP/no-GM/no-Redi). It prescribes analytic T/S/UV/eta/
# w_e/w_i/Av/stress_surf and drives the REAL pressure_bv -> pgf -> compute_vel_rhs ->
# visc_filt_bidiff -> impl_vert_visc_ale -> init/compute_ssh_rhs_ale + solve_ssh_ale, dumping
# all 28+ fields (incl. ssh_stiff_diag / ssh_Aeta / ssh_rhs / d_eta). Used to LOCALIZE the
# M2.11c-1 lifecycle CG 1-ULP divergence (is it the stiffness matrix or the CG?) on CORE2.
# step_per_day=48 -> the oracle's init_stiff dt = 86400/48 = 1800 s (FESOM3 must match via
# FESOM3_STEP_PER_DAY=48).
#
#   tools/run_pressuredump_core2.sh [run_dir] [out_file]
set -euo pipefail

F2=/home/a/a270088/port2/fesom2
RUN="${1:-/scratch/a/a270088/pressuredump_core2}"
OUT="${2:-$RUN/pressure_f2.bin}"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1

rm -rf "$RUN"; mkdir -p "$RUN"
cp "$F2"/work_core/namelist.* "$RUN"/
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x
printf '0 1 1948\n0 1 1948\n' > "$RUN"/fesom.clock

python3 - "$RUN" <<'PY'
import re, sys
run = sys.argv[1]
p = run + '/namelist.config'; s = open(p).read()
s = re.sub(r"ResultPath\s*=\s*'[^']*'", "ResultPath       = './'", s, count=1)
s = re.sub(r"which_ALE\s*=\s*'zlevel'", "which_ALE          = 'linfs'", s, count=1)
open(p, 'w').write(s)
p = run + '/namelist.oce'; s = open(p).read()
s = re.sub(r"mix_scheme\s*=\s*'KPP'", "mix_scheme         = 'PP'", s, count=1)
s = re.sub(r"Fer_GM\s*=\s*\.true\.",  "Fer_GM             = .false.", s, count=1)
s = re.sub(r"Redi\s*=\s*\.true\.",    "Redi               = .false.", s, count=1)
open(p, 'w').write(s)
PY

cd "$RUN"
export FESOM_PRESSURE_DUMP="$OUT"
ulimit -s unlimited
echo "run_pressuredump_core2: 1 rank -> $OUT"
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 ./fesom.x > "$RUN/run.log" 2>&1 || true
if [[ -f "$OUT" ]]; then
    echo "run_pressuredump_core2: done -> $OUT ($(stat -c%s "$OUT") bytes)"
else
    echo "run_pressuredump_core2: NO DUMP — see $RUN/run.log"; tail -30 "$RUN/run.log"; exit 1
fi
