#!/usr/bin/env bash
# M2.12a end-to-end MULTI-RANK geometry byte-gate: dump FESOM2 geometry (NP ranks,
# dist_<NP>, pi) and FESOM3 geometry (NP ranks, dist_<NP>, pi), and compare per-rank
# on OWNED entries (1..myDim_*) for max|delta|=0. Both share the same myList ordering
# (same dist_<NP>/), so local index i == the same global id and the per-node area
# accumulation order matches (the L8 same-partition rule). This validates the FESOM3
# local-mesh remap (mod_mesh_read.read_mesh_local) + partition-aware compute_geometry.
#
#   tools/run_geom_gate_multirank.sh [np] [run_dir]
set -euo pipefail
F3=/home/a/a270088/fesom3
NP="${1:-2}"
RUN="${2:-/scratch/a/a270088/geomdump_mr${NP}_pi}"
MESH=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi

echo "[1/3] FESOM2 oracle geometry dump ($NP-rank pi, dist_$NP)"
bash "$F3/tools/run_geomdump_multirank.sh" "$RUN" "$RUN/geom_f2.bin" "$NP" | tail -2

echo "[2/3] FESOM3 geometry dump ($NP-rank pi, dist_$NP)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$MESH" FESOM3_GEOM_OUT="$RUN/geom_f3.bin"
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n "$NP" \
    "$F3/build_intel_dp/bin/fesom_geomdump" >/dev/null 2>&1

echo "[3/3] compare per-rank (OWNED entries, max|delta|=0)"
fail=0
for r in $(seq 0 $((NP-1))); do
    rs=$(printf '%05d' "$r")
    out=$(python3 "$F3/tools/geom_diff.py" "$RUN/geom_f2.bin.$rs" "$RUN/geom_f3.bin.$rs" 2>&1)
    if echo "$out" | grep -q "PASS (max|delta|=0)"; then
        echo "  rank $rs: PASS (all fields max|delta|=0)"
    else
        echo "  rank $rs: FAIL"; echo "$out" | grep -iE "FAIL|differ" | head; fail=1
    fi
done
if [[ $fail -eq 0 ]]; then
    echo "MULTI-RANK GEOMETRY BYTE-GATE ($NP ranks): PASS — all ranks max|delta|=0"
else
    echo "MULTI-RANK GEOMETRY BYTE-GATE ($NP ranks): FAIL"; exit 1
fi
