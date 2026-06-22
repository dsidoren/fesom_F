#!/usr/bin/env bash
# M2.12c-1 MULTI-RANK dynamics byte-gate. Runs the FESOM2 oracle (NP ranks, dist_<NP>, pi:
# the step-dump shim prescribes the analytic state + drives the REAL oce_timestep_ale, whose
# built-in dump_shim_record_node emits per-rank gid-keyed substep dumps) and the FESOM3
# dynamics chain UP TO compute_ssh_rhs_ale (fesom_stepdump_mr, NP ranks, dist_<NP>), then
# compares the per-rank PRE-SSH-solve probe fields — density/pressure/bvfreq (substep 1),
# Kv (4), ssh_rhs (8) — for max|delta|=0. The oracle's post-SSH substeps (9/11/12/13/15/16)
# and SW_AB (2) are ignored (FESOM3 stops before the un-lifted CG; those are M2.12c-2).
# Both share the same dist_<NP>/ myList, so each global probe id is owned by the same rank
# on both codes (the L8 same-partition rule) and dump_diff matches them by global id.
#
#   tools/run_stepdyn_gate_multirank.sh [np] [run_dir]
set -euo pipefail
F3=/home/a/a270088/fesom3
F2=/home/a/a270088/port2/fesom2
NP="${1:-2}"
RUN="${2:-/scratch/a/a270088/stepdyn_mr${NP}_pi}"
MESH="$F2/tests/data/MESHES/pi"

echo "[1/3] FESOM2 oracle: REAL oce_timestep_ale ($NP-rank pi, dist_$NP)"
source "$F3/env.sh" intel >/dev/null 2>&1
rm -rf "$RUN"; mkdir -p "$RUN"
cp "$F2"/work_pi/namelist.* "$RUN"/
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x
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
( cd "$RUN"
  export FESOM_STEP_DUMP=1 FESOM_DUMP_FILE="$RUN/step_f2" FESOM_DUMP_MAXSTEPS=1
  ulimit -s unlimited
  mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n "$NP" ./fesom.x > "$RUN/run_f2.log" 2>&1 || true )
if ! ls "$RUN"/step_f2.* >/dev/null 2>&1; then
    echo "  NO ORACLE DUMP — see $RUN/run_f2.log"; tail -30 "$RUN/run_f2.log"; exit 1
fi
echo "  oracle dumps: $(ls "$RUN"/step_f2.* | tr '\n' ' ')"

echo "[2/3] FESOM3: dynamics chain -> ssh_rhs ($NP-rank pi, dist_$NP)"
export FESOM3_MESH_DIR="$MESH" FESOM_DUMP_FILE="$RUN/step_f3" FESOM_DUMP_MAXSTEPS=1
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n "$NP" \
    "$F3/build_intel_dp/bin/fesom_stepdump_mr" > "$RUN/run_f3.log" 2>&1 || \
    { echo "  FESOM3 run failed"; tail -30 "$RUN/run_f3.log"; exit 1; }

echo "[3/3] compare per-rank (gid-keyed; substeps 1/4/8 pre-SSH + 9=d_eta from the c-2 CG)"
# ignore the oracle's post-SSH-solve substeps (FESOM3 stops at d_eta; ALE/tracers = c-3).
# substep 9 (d_eta) is now COMPARED — the M2.12c-2 free-surface CG byte-gate.
IGN="--ignore-substep=2 --ignore-substep=11 --ignore-substep=12 \
     --ignore-substep=13 --ignore-substep=15 --ignore-substep=16"
python3 "$F3/tools/dump_diff.py" "$RUN/step_f2" "$RUN/step_f3" --glob $IGN
