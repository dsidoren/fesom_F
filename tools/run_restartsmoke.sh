#!/usr/bin/env bash
# Restart Stage 3 (Task 3.2) GATE: the mod_io_restart WRITE path via fesom_restartsmoke (a real pi
# mesh, NO ocean — fast). It registers ONE node field (value(g)==g) as a live pointer and writes a
# checkpoint folder fesom.2000.001.03600/ holding eta_n.zarr (a single-variable single-entity SNAPSHOT
# store — NO time dim) + checkpoint.json. tools/zarr_diff.py --restart asserts the store opens in
# xarray (dims (nod2,), embedded finite lon/lat, _ARRAY_DIMENSIONS, value(g)==g) and that
# checkpoint.json parses. Run at np 1 AND 2 (the np=2 store must be byte-identical-in-shape +
# value-identical to np=1 — partition independence of the canonical write).
#
#   tools/run_restartsmoke.sh           # default: np 1 2
#   tools/run_restartsmoke.sh 1         # just 1-rank
set -euo pipefail
F3="${F3:-/home/a/a270088/fesom3}"          # override to run from a git worktree checkout
BUILD="${BUILD:-$F3/build_intel_dp}"
PY=/work/ab0995/a270088/mambaforge/bin/python3
PIMESH=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi
RUN="${RUN:-/scratch/a/a270088/restartsmoke}"
NPS="${*:-1 2}"
source "$F3/env.sh" intel >/dev/null 2>&1

# Reconfigure (re-glob the new src/io + src/drivers files) then incremental build. The explicit
# reconfigure is needed the first time a brand-new target/file appears.
cmake "$BUILD" >/dev/null
cmake --build "$BUILD" --target fesom_restartsmoke -j 4

export FESOM3_MESH_DIR="$PIMESH"
export FESOM3_CHUNK_HORIZ="${FESOM3_CHUNK_HORIZ:-1000}"   # small => multiple chunks + writer subset
ulimit -s unlimited
MPIRUN=(mpirun --mca pml ob1 --mca btl self,vader --oversubscribe)

for NP in $NPS; do
    DIR="$RUN/np$NP"
    rm -rf "$DIR"; mkdir -p "$DIR"
    export FESOM3_RESTARTSMOKE_DIR="$DIR"
    echo "=== fesom_restartsmoke np=$NP ==="
    "${MPIRUN[@]}" -n "$NP" "$BUILD/bin/fesom_restartsmoke" 2>&1 | grep -E "RESTARTSMOKE|nod2D"
    "$PY" "$F3/tools/zarr_diff.py" --restart "$DIR"
done

echo "run_restartsmoke: GATE GREEN (np: $NPS)"
echo "  manual ushow (needs a display): \$ushow $RUN/np1/fesom.2000.001.03600/eta_n.zarr -m fesom.mesh.diag.zarr"
