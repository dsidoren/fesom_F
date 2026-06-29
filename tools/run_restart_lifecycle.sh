#!/usr/bin/env bash
# Restart Stage 5 (Tasks 5.1/5.2/5.3) FUNCTIONAL 2-segment gate. Proves the lifecycle WIRING resumes
# cleanly (this is the FUNCTIONAL "does it resume?" gate — the BYTE-EXACT max|Δ|=0 gate is Task 6.1).
#
#   Segment 1 (cold): run K steps of the REAL fully-native CORE2 lifecycle with FESOM3_RESTART=<dir>
#       and a step cadence => writes a PERIODIC mid-run checkpoint + the last-step checkpoint, the
#       <runid>.clock (clock_finish after every write) and restart.latest under <dir>.
#   Segment 2 (resume): a FRESH process with FESOM3_RESTART_IN=<dir> (+FESOM3_RESTART=<dir>) over the
#       remaining steps => must DETECT restart (restart_mode + r_restart, the RESTART RUN banner),
#       restart_read the newest checkpoint, run the rest, and EXIT 0.
#
# Asserts: seg-1 emits >=2 checkpoint folders (periodic cadence fired mid-run) + restart.latest + a
# chained .clock; seg-2 prints "RESTART MODE"/"RESTART RUN"/"state restored" and its clock_init time
# equals seg-1's newest-checkpoint time (the .clock chained); seg-2 exits 0.
#
#   tools/run_restart_lifecycle.sh [run_dir] [NP]      # default NP=1
set -uo pipefail
F3="${F3:-/home/a/a270088/fesom3/.claude/worktrees/restart-checkpoint}"
BUILD="${BUILD:-$F3/build_intel_dp}"
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
NCAR=/home/a/a270088/port2/fesom2/test/input/global
RUN="${1:-/scratch/a/a270088/restart_lifecycle}"
NP="${2:-1}"
K=3            # segment-1 steps (checkpoints at step 2 [periodic] and step 3 [last])
CAD=2         # restart cadence in steps (FESOM3_RESTART_UNIT=s)
REMAIN=3      # segment-2 steps

source "$F3/env.sh" intel >/dev/null 2>&1
rm -rf "$RUN"; mkdir -p "$RUN"
CKDIR="$RUN/ckpt"; mkdir -p "$CKDIR"

export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE"
export FESOM3_FORCING_DIR="$NCAR"
export FESOM3_RUNOFF_FILE="$POOL/CORE2_runoff.nc"
export FESOM3_SSS_FILE="$POOL/PHC2_salx.nc"
export FESOM3_WHICHEVP=0
ulimit -s unlimited
MPIRUN=(mpirun --mca pml ob1 --mca btl self,vader --mca btl_vader_single_copy_mechanism none --oversubscribe)
BIN="$BUILD/bin/fesom_lifecycle_native_mr"

echo "############ Segment 1: COLD, $K steps, checkpoint every $CAD steps (NP=$NP) ############"
FESOM3_RESTART="$CKDIR" FESOM3_RESTART_IN="$CKDIR" \
FESOM3_RESTART_LENGTH="$CAD" FESOM3_RESTART_UNIT=s \
FESOM3_NSTEPS="$K" \
    "${MPIRUN[@]}" -n "$NP" "$BIN" 2>&1 | tee "$RUN/seg1.log" \
    | grep -E 'nod2D|INITIALISATION|RESTART|registered|checkpoint written|final checkpoint|done' || true
s1=${PIPESTATUS[0]}

echo
echo "---- segment 1 artifacts under $CKDIR ----"
ls -1 "$CKDIR" | sed 's/^/    /'
NFOLD=$(ls -1d "$CKDIR"/fesom.*/ 2>/dev/null | wc -l)
CLOCKF=$(ls -1 "$CKDIR"/*.clock 2>/dev/null | head -1)   # runid-named (e.g. test1.clock); discover it
S1LATEST=$(cat "$CKDIR/restart.latest" 2>/dev/null)
echo "    checkpoint folders: $NFOLD"
echo "    restart.latest -> $S1LATEST"
echo "    $(basename "$CLOCKF"):"; sed 's/^/        /' "$CLOCKF" 2>/dev/null
CKTIME=$(awk 'NR==2{printf "%d", $1}' "$CLOCKF" 2>/dev/null)        # new-clock sec-of-day (.clock line 2)
CKTAG=$(echo "$S1LATEST" | awk -F. '{printf "%d", $NF}')           # newest-checkpoint sec-of-day (folder tag)
echo "    chained new-clock sec-of-day (line 2) = $CKTIME ; newest-checkpoint tag = $CKTAG"

echo
echo "############ Segment 2: RESUME, $REMAIN steps, FRESH process (NP=$NP) ############"
FESOM3_RESTART="$CKDIR" FESOM3_RESTART_IN="$CKDIR" \
FESOM3_RESTART_LENGTH="$CAD" FESOM3_RESTART_UNIT=s \
FESOM3_NSTEPS="$REMAIN" \
    "${MPIRUN[@]}" -n "$NP" "$BIN" 2>&1 | tee "$RUN/seg2.log" \
    | grep -E 'nod2D|RESTART MODE|RESTART RUN|clock restarted|state restored|registered|checkpoint written|final checkpoint|done' || true
s2=${PIPESTATUS[0]}

echo
echo "############ Assertions ############"
rc=0
# (1) seg-1 wrote >=2 checkpoint folders (periodic mid-run cadence fired) + restart.latest
[ "$NFOLD" -ge 2 ] && echo "PASS  seg-1 wrote $NFOLD checkpoint folders (periodic + last-step)" \
                   || { echo "FAIL  seg-1 wrote $NFOLD folders (<2)"; rc=1; }
[ -s "$CKDIR/restart.latest" ] && echo "PASS  restart.latest present -> $(cat "$CKDIR/restart.latest")" \
                               || { echo "FAIL  restart.latest missing/empty"; rc=1; }
# (2) seg-2 detected restart + restored state
grep -q 'RESTART MODE'                "$RUN/seg2.log" && echo "PASS  seg-2 detected RESTART MODE (skipped cold .clock overwrite)" || { echo "FAIL  seg-2 no RESTART MODE"; rc=1; }
grep -q 'THIS IS A RESTART RUN'       "$RUN/seg2.log" && echo "PASS  seg-2 r_restart=.true. (RESTART RUN banner)"                || { echo "FAIL  seg-2 no RESTART RUN banner"; rc=1; }
grep -q 'state restored'              "$RUN/seg2.log" && echo "PASS  seg-2 restart_read restored the checkpoint"                 || { echo "FAIL  seg-2 did not restore"; rc=1; }
# (3) .clock chaining: seg-1 .clock line-2 sec-of-day == newest-checkpoint folder tag == seg-2 clock_init
S2TIME=$(grep 'clock restarted at time:' "$RUN/seg2.log" | head -1 | sed -E 's/.*time:[[:space:]]*([0-9.]+).*/\1/')
S2TI=$(printf "%.0f" "${S2TIME:-nan}" 2>/dev/null || echo nan)
if [ -n "$CKTIME" ] && [ "$CKTIME" = "$CKTAG" ] && [ "$S2TI" = "$CKTIME" ]; then
    echo "PASS  .clock chained: seg-1 .clock=$CKTIME == checkpoint tag=$CKTAG == seg-2 clock_init=$S2TI"
else
    echo "FAIL  .clock chain mismatch: seg-1 .clock=$CKTIME tag=$CKTAG seg-2 clock_init=$S2TI"; rc=1
fi
# (4) cold seg-1 must NOT be a restart run (sanity that detection is conditional)
grep -q 'INITIALISATION RUN' "$RUN/seg1.log" && echo "PASS  seg-1 was a cold INITIALISATION run (detection is conditional)" || { echo "FAIL  seg-1 was not cold"; rc=1; }
# (5) both segments exited 0
[ "$s1" = 0 ] && echo "PASS  segment 1 exit 0" || { echo "FAIL  segment 1 exit $s1"; rc=1; }
[ "$s2" = 0 ] && echo "PASS  segment 2 exit 0" || { echo "FAIL  segment 2 exit $s2"; rc=1; }

echo
[ "$rc" = 0 ] && echo "run_restart_lifecycle: FUNCTIONAL GATE GREEN (NP=$NP)  — clean 2-segment resume" \
             || echo "run_restart_lifecycle: FUNCTIONAL GATE FAILED (NP=$NP)"
exit "$rc"
