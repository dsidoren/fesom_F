#!/usr/bin/env bash
# M9 Stage 2 GATE (Task 2.1 + 2.3): field output via mod_io_means. Runs fesom_outputsmoke (a real pi
# mesh, NO ocean — fast) which writes per-variable-per-year Zarr stores of two synthetic node fields
# whose value at (record k, CANONICAL node g) follows a partition-independent generator formula.
# tools/zarr_diff.py --output then asserts every value == the formula (writer self-consistency: proves
# registry + growing time dim + append + decomp + CF time + lon/lat embed), and --output-cmp asserts
# dist_2 == dist_8 == 1-rank (partition independence). No FESOM2 oracle needed (the model STATE is
# already byte-exact vs FESOM2 through M8; this gate proves the WRITER serializes it faithfully).
#
#   tools/run_output_gate.sh           # default: np 1,2,8 + partition-independence
#   tools/run_output_gate.sh 1         # just 1-rank
set -euo pipefail
F3=/home/a/a270088/fesom3
BUILD="${BUILD:-$F3/build_intel_dp}"
PY=/work/ab0995/a270088/mambaforge/bin/python3
PIMESH=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi
RUN="${RUN:-/scratch/a/a270088/outsmoke}"
NPS="${*:-1 2 8}"
source "$F3/env.sh" intel >/dev/null 2>&1

cmake "$BUILD" >/dev/null
cmake --build "$BUILD" --target fesom_outputsmoke -j 4 >/dev/null

export FESOM3_MESH_DIR="$PIMESH"
export FESOM3_OUTSMOKE_NREC="${FESOM3_OUTSMOKE_NREC:-5}"
export FESOM3_CHUNK_HORIZ="${FESOM3_CHUNK_HORIZ:-1000}"   # small => multiple chunks + writer subset
ulimit -s unlimited
MPIRUN=(mpirun --mca pml ob1 --mca btl self,vader --oversubscribe)

# ---- GEOGRAPHIC frame (default): scalars + fld_u/fld_v r2g-rotated (numpy reference) -------------
export FESOM3_VEC_FRAME=geographic
for NP in $NPS; do
    DIR="$RUN/np$NP"
    rm -rf "$DIR"; mkdir -p "$DIR"
    export FESOM3_OUTSMOKE_DIR="$DIR"
    echo "=== fesom_outputsmoke np=$NP (geographic) ==="
    "${MPIRUN[@]}" -n "$NP" "$BUILD/bin/fesom_outputsmoke" 2>&1 | grep -E "OUTPUTSMOKE|nod2D"
    "$PY" "$F3/tools/zarr_diff.py" --output "$DIR" --frame geographic
done

# partition-independence (geographic): compare consecutive rank counts (incl. the rotated vectors)
prev=""
for NP in $NPS; do
    if [ -n "$prev" ]; then
        "$PY" "$F3/tools/zarr_diff.py" --output-cmp "$RUN/np$prev" "$RUN/np$NP"
    fi
    prev="$NP"
done

# ---- NATIVE frame: fld_u/fld_v must equal the raw (unrotated) generator (max|Δ|=0) --------------
export FESOM3_VEC_FRAME=native
NNAT=""
for NP in $NPS; do [ "$NP" -le 2 ] && NNAT="$NNAT $NP"; done   # np 1[,2] suffice for native passthrough
for NP in $NNAT; do
    DIR="$RUN/np${NP}_native"
    rm -rf "$DIR"; mkdir -p "$DIR"
    export FESOM3_OUTSMOKE_DIR="$DIR"
    echo "=== fesom_outputsmoke np=$NP (native) ==="
    "${MPIRUN[@]}" -n "$NP" "$BUILD/bin/fesom_outputsmoke" 2>&1 | grep -E "OUTPUTSMOKE|nod2D"
    "$PY" "$F3/tools/zarr_diff.py" --output "$DIR" --frame native
done
# native partition-independence (raw vectors identical at any rank count)
if echo "$NNAT" | grep -qw 2; then
    "$PY" "$F3/tools/zarr_diff.py" --output-cmp "$RUN/np1_native" "$RUN/np2_native"
fi

echo "run_output_gate: GATE GREEN (np: $NPS; frames: geographic+native)"
echo "  manual ushow (needs a display): $F3/../ushow/ushow $RUN/np1/fld_u.fesom.2000.zarr"
