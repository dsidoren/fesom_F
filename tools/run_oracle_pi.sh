#!/usr/bin/env bash
# Run the instrumented FESOM2 v2.7.3 oracle on the pi mesh and produce a gid-keyed
# substep dump that tools/dump_diff.py compares against FESOM3 (byte-gate D8).
#
# PROVEN working 2026-06-19 (FESOM2 git SHA 9271ae92): a 2-step, 2-rank pi run takes
# ~0.25 s and writes <prefix>.<rank> dumps (density/pressure/bvfreq/sw_alpha-beta/Kv/
# ssh_rhs/d_eta/hbar/eta_n/w/T/S/hnode). The dump format is byte-identical to FESOM3
# src/infra/mod_dump.F90.
#
#   tools/run_oracle_pi.sh [run_dir] [nsteps] [nranks] [dump_prefix]
#
# NOTE: this uses the SHIPPED pi namelists (default mixing/GM). For the M2 byte-gate,
# switch to the reduced namelist (mix_scheme='PP', Fer_GM=.false., Redi=.false.,
# which_ALE='linfs', opt_visc=7) per docs/HANDOFF.md.
set -euo pipefail

F2=/home/a/a270088/port2/fesom2
RUN="${1:-/scratch/a/a270088/oracle_pi}"
NSTEPS="${2:-2}"
NRANKS="${3:-2}"
DUMP="${4:-$RUN/dump}"

# Toolchain matching the prebuilt fesom.x (Intel + OpenMPI + netCDF).
# shellcheck disable=SC1091
source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1

rm -rf "$RUN"; mkdir -p "$RUN"
cp "$F2"/work_pi/namelist.* "$RUN"/
ln -sf "$F2"/bin/fesom.x "$RUN"/fesom.x
printf '0 1 1948\n0 1 1948\n' > "$RUN"/fesom.clock   # fresh start, forcing year 1948

python3 - "$RUN" "$F2" "$NSTEPS" <<'PY'
import re, sys
run, f2, nsteps = sys.argv[1], sys.argv[2], sys.argv[3]
p = run + '/namelist.config'
s = open(p).read()
s = re.sub(r'run_length\s*=\s*\S+',      f'run_length        = {nsteps}', s, count=1)
s = re.sub(r"run_length_unit\s*=\s*\S+", "run_length_unit   = 's'",        s, count=1)
s = re.sub(r"ResultPath\s*=\s*'[^']*'",  "ResultPath       = './'",        s, count=1)
s = re.sub(r"ClimateDataPath\s*=\s*'[^']*'",
           f"ClimateDataPath  = '{f2}/test/input/global/'", s, count=1)
open(p, 'w').write(s)
PY

cd "$RUN"
export FESOM_DUMP_FILE="$DUMP" FESOM_DUMP_MAXSTEPS="$NSTEPS"
ulimit -s unlimited
echo "run_oracle_pi: $NRANKS ranks, $NSTEPS steps -> ${DUMP}.<rank>"
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n "$NRANKS" ./fesom.x > "$RUN/run.log" 2>&1
echo "run_oracle_pi: done. dumps:"
ls -la "${DUMP}".* 2>/dev/null || { echo "NO DUMP — see $RUN/run.log"; tail -5 "$RUN/run.log"; exit 1; }
