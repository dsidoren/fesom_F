#!/usr/bin/env bash
# M9 Stage 1 GATE (1-rank, pi): write fesom.mesh.diag.zarr via fesom_meshdiagdump and compare it
# vs the FESOM2 fesom.mesh.diag.nc oracle (test/output_pi), over the emitted subset, in canonical
# order — integer/connectivity EXACT, floats max|Δ|=0. Covers Task 1.2a (coords/connectivity/UGRID)
# + 1.2b (derived: areas, nlevels, nod_in_elem2D, edge_cross_dxdy, gradient_sca, zbar bottoms,
# nod_area). face_edges/face_links/gradient_vec are DEFERRED (FESOM3 never builds
# elem_edges/elem_neighbors/gradient_vec; nod_part/elem_part are partition descriptors -> excluded).
#
#   tools/run_meshdiag_gate.sh [nranks]   (nranks default 1; >1 = partition-independence, Task 1.3)
set -euo pipefail
F3=/home/a/a270088/fesom3
NP="${1:-1}"
BUILD="${BUILD:-$F3/build_intel_dp}"
PY=/work/ab0995/a270088/mambaforge/bin/python3
PIMESH=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi
REF=/home/a/a270088/port2/fesom2/test/output_pi/fesom.mesh.diag.nc
RUN="${RUN:-/scratch/a/a270088/meshdiag_pi}"
mkdir -p "$RUN"
source "$F3/env.sh" intel >/dev/null 2>&1

cmake "$BUILD" >/dev/null
cmake --build "$BUILD" --target fesom_meshdiagdump -j 4 >/dev/null

STORE="$RUN/fesom.mesh.diag.np${NP}.zarr"
rm -rf "$STORE"
export FESOM3_MESH_DIR="$PIMESH"
export FESOM3_MESHDIAG_OUT="$STORE"
export FESOM3_CHUNK_HORIZ="${FESOM3_CHUNK_HORIZ:-1000}"   # small => multiple chunks + writer subset
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n "$NP" \
    "$BUILD/bin/fesom_meshdiagdump" 2>&1 | grep -E "nod2D|MESHDIAGDUMP"

"$PY" "$F3/tools/zarr_diff.py" --meshdiag "$STORE" "$REF"

echo "run_meshdiag_gate: GATE GREEN (np=$NP)"
echo "  manual ushow (needs a display): $F3/../ushow/ushow $STORE -m $STORE"
