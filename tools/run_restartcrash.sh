#!/usr/bin/env bash
# Restart Stage 3 (Task 3.3) crash-safety GATE: atomic checkpoint finalize + restart.latest pointer +
# keep-N prune, proven WITHOUT a real crash via fesom_restartcrash (a real pi mesh, NO ocean — fast).
# The driver writes C1 atomically, injects a stray fesom.*.tmp/ and a finalized-but-unpointed
# fesom.<later>/ dir, asserts restart_resolve_latest STILL follows restart.latest to C1, then writes C2
# (keep-N=1) and asserts the pointer flipped, C1 was pruned, C2 + the protected later dir remain, and
# the stray .tmp is left untouched. Runs at np 1 AND 2 (atomic finalize is rank-0; verify it still works
# collectively at np=2). This script also dumps restart.latest + the dir listing as independent evidence.
#
#   tools/run_restartcrash.sh           # default: np 1 2
#   tools/run_restartcrash.sh 1         # just 1-rank
set -euo pipefail
F3="${F3:-/home/a/a270088/fesom3}"          # override to run from a git worktree checkout
BUILD="${BUILD:-$F3/build_intel_dp}"
PIMESH=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi
RUN="${RUN:-/scratch/a/a270088/restartcrash}"
NPS="${*:-1 2}"
source "$F3/env.sh" intel >/dev/null 2>&1

# Reconfigure (re-glob the new src/drivers file) then incremental build.
cmake "$BUILD" >/dev/null
cmake --build "$BUILD" --target fesom_restartcrash -j 4

export FESOM3_MESH_DIR="$PIMESH"
export FESOM3_CHUNK_HORIZ="${FESOM3_CHUNK_HORIZ:-1000}"
ulimit -s unlimited
MPIRUN=(mpirun --mca pml ob1 --mca btl self,vader --oversubscribe)

for NP in $NPS; do
    DIR="$RUN/np$NP"
    rm -rf "$DIR"; mkdir -p "$DIR"
    export FESOM3_RESTARTCRASH_DIR="$DIR"
    LOG="$DIR/.crash.log"
    echo "=== fesom_restartcrash np=$NP ==="
    "${MPIRUN[@]}" -n "$NP" "$BUILD/bin/fesom_restartcrash" >"$LOG" 2>&1 || true
    grep -E "PASS|FAIL|pointer|resolve|RESTARTCRASH|nod2D" "$LOG"
    echo "  --- restart.latest -> $(cat "$DIR/restart.latest")"
    echo "  --- restart dir (C1 pruned, C2 + later kept, .tmp ignored):"
    ls -1 "$DIR" | grep -v '^\.crash\.log$' | sed 's/^/        /'
    # hard-assert the gate line is present (driver error-stops nonzero on any failed assertion).
    grep -q "RESTARTCRASH OK" "$LOG"
done

echo "run_restartcrash: GATE GREEN (np: $NPS)"
