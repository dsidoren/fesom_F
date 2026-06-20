#!/usr/bin/env bash
# M2.11a CORE2 geometry byte-gate: dump FESOM2 geometry (1-rank CORE2), dump FESOM3
# geometry (1-rank CORE2), compare for max|delta|=0. Scales the proven pi geom gate
# 40x (nod2D=126858/elem2D=244659/edge2D=371644/nl=48) and EMPIRICALLY closes the L8
# CW-orientation-swap deferral (CORE2 needs ~244654/244659 enforce_cw_orientation
# swaps vs pi's 0 — the reorder path is exercised + gated for the first time).
#
# Prereqs (one-time):
#   - FESOM2 build/lib64/libfesom.so has the geom dump shim (oce_mesh.F90 -> geom_dump_write).
#   - CORE2 mesh has the hand-crafted dist_1/ (tools/make_dist1.py <core2_mesh_dir>).
#   - FESOM3 built: build_intel_dp/bin/fesom_geomdump.
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
RUN="${1:-/scratch/a/a270088/geomdump_core2}"

echo "[1/3] FESOM2 oracle geometry dump (1-rank CORE2)"
bash "$F3/tools/run_geomdump_core2.sh" "$RUN" "$RUN/geom_f2.bin" >/dev/null

echo "[2/3] FESOM3 geometry dump (1-rank CORE2)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH"
export FESOM3_GEOM_OUT="$RUN/geom_f3.bin"
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_geomdump" 2>&1 | grep -E 'CW swaps|nod2D'

echo "[3/3] compare"
python3 "$F3/tools/geom_diff.py" "$RUN/geom_f2.bin" "$RUN/geom_f3.bin"
