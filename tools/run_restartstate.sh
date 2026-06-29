#!/usr/bin/env bash
# Restart Stage 3 (Task 3.4) GATE: restart_register_state — the FULL oce+ice prognostic field set incl.
# EVP sigma and the real(MP) mesh%hbar/hnode — via fesom_restartstate (a real pi mesh, NO ocean). It
# allocates correctly-shaped REAL t_dyn/t_tracer/t_ice + mesh%hbar/hnode, fills owned slots with a
# partition-independent formula (g for 2-D, g+0.5*L for 3-D), registers every field via
# restart_register_state, and writes one checkpoint folder fesom.2000.001.03600/. tools/zarr_diff.py
# --restart-state asserts every expected store exists with the right entity x level-kind shape +
# _ARRAY_DIMENSIONS, embedded finite lon/lat, monotonic nz/nz1, finite data, value==formula (max|Δ|=0,
# so the MP->WP staging of hbar/hnode is proven lossless), and element stores ~2x the node count.
# Runs AB_order=2 (26 stores, tke on) AND AB_order=3 (30 stores: +urhs_AB3/vrhs_AB3/temp_M2/salt_M2),
# at np 1 AND 2, plus a np1-vs-np2 partition-independence compare (canonical write is np-independent).
#
#   tools/run_restartstate.sh           # default: np 1 2
#   tools/run_restartstate.sh 1         # just 1-rank
set -euo pipefail
F3="${F3:-/home/a/a270088/fesom3}"          # override to run from a git worktree checkout
BUILD="${BUILD:-$F3/build_intel_dp}"
PY=/work/ab0995/a270088/mambaforge/bin/python3
PIMESH=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi
RUN="${RUN:-/scratch/a/a270088/restartstate}"
NPS="${*:-1 2}"
source "$F3/env.sh" intel >/dev/null 2>&1

cmake "$BUILD" >/dev/null
cmake --build "$BUILD" --target fesom_restartstate -j 4

export FESOM3_MESH_DIR="$PIMESH"
export FESOM3_CHUNK_HORIZ="${FESOM3_CHUNK_HORIZ:-1000}"   # small => multiple chunks + writer subset
ulimit -s unlimited
MPIRUN=(mpirun --mca pml ob1 --mca btl self,vader --mca btl_vader_single_copy_mechanism none --oversubscribe)

for AB in 2 3; do
    [ "$AB" = 2 ] && NEXP=28 || NEXP=32
    export FESOM3_AB_ORDER="$AB"
    for NP in $NPS; do
        DIR="$RUN/ab${AB}_np$NP"
        rm -rf "$DIR"; mkdir -p "$DIR"
        export FESOM3_RESTARTSTATE_DIR="$DIR"
        echo "=== fesom_restartstate AB_order=$AB np=$NP (expect $NEXP stores) ==="
        "${MPIRUN[@]}" -n "$NP" "$BUILD/bin/fesom_restartstate" 2>&1 | grep -E "RESTARTSTATE|fields="
        "$PY" "$F3/tools/zarr_diff.py" --restart-state "$DIR" --ab-order "$AB" | tail -3
    done
done

# partition independence: AB2 np1 vs np2 checkpoint folders must be value-identical (max|Δ|=0)
if [[ " $NPS " == *" 1 "* && " $NPS " == *" 2 "* ]]; then
    echo "=== partition-independence: AB2 np1 vs np2 (max|Δ|=0) ==="
    "$PY" "$F3/tools/zarr_diff.py" --output-cmp \
        "$RUN/ab2_np1/fesom.2000.001.03600" "$RUN/ab2_np2/fesom.2000.001.03600" | tail -1
fi

echo "run_restartstate: GATE GREEN (np: $NPS)"
