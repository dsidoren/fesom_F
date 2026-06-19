#!/usr/bin/env bash
# End-to-end M1 geometry byte-gate (closes deferred M0.7): dump FESOM2 geometry
# (1-rank pi), dump FESOM3 geometry (1-rank pi), and compare for max|delta|=0.
#
# Prereqs (one-time, see docs/HANDOFF.md):
#   - FESOM2 rebuilt with src/fesom_geom_dump.F90 wired into mesh_setup
#     (build/bin/fesom.x -> build/lib64/libfesom.so).
#   - pi mesh has the hand-crafted dist_1/ (tools/ generated it; rpart.out has the
#     npes + counts + identity node->contiguous mapping; com_info uses blank lines
#     for the zero-size halo arrays).
#   - FESOM3 built: build_intel_dp/bin/fesom_geomdump.
set -euo pipefail
F3=/home/a/a270088/fesom3
RUN="${1:-/scratch/a/a270088/geomdump_pi}"

echo "[1/3] FESOM2 oracle geometry dump (1-rank pi)"
bash "$F3/tools/run_geomdump_pi.sh" "$RUN" "$RUN/geom_f2.bin" >/dev/null

echo "[2/3] FESOM3 geometry dump (1-rank pi)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_GEOM_OUT="$RUN/geom_f3.bin"
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_geomdump" >/dev/null 2>&1

echo "[3/3] compare"
python3 "$F3/tools/geom_diff.py" "$RUN/geom_f2.bin" "$RUN/geom_f3.bin"
