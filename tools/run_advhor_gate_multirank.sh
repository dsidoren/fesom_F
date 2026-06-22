#!/usr/bin/env bash
# M2.12b end-to-end MULTI-RANK tracer-advection byte-gate: dump FESOM2 advection (NP
# ranks, dist_<NP>, pi) and FESOM3 advection (NP ranks, dist_<NP>, pi), and compare per
# rank on OWNED entries (1..myDim_*) for max|delta|=0. Both share the same myList order
# (same dist_<NP>/) so local index i == the same global id and the per-node/per-edge
# accumulation order matches (the L8 same-partition rule). Validates the multi-rank
# advection assembly (find_neighbors halo dance + exchange_elem_full(tr_xy/elem_area) +
# exchange_nod(fct_LO/fct_plus/fct_minus) + owned/halo loop bounds).
#
#   tools/run_advhor_gate_multirank.sh [np] [run_dir]
set -euo pipefail
F3=/home/a/a270088/fesom3
NP="${1:-2}"
RUN="${2:-/scratch/a/a270088/advhordump_mr${NP}_pi}"
MESH=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi

echo "[1/3] FESOM2 oracle advhor dump ($NP-rank pi, dist_$NP)"
bash "$F3/tools/run_advhordump_multirank.sh" "$RUN" "$RUN/advhor_f2.bin" "$NP" | tail -2

echo "[2/3] FESOM3 advhor dump ($NP-rank pi, dist_$NP)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$MESH" FESOM3_ADVHOR_OUT="$RUN/advhor_f3.bin"
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n "$NP" \
    "$F3/build_intel_dp/bin/fesom_advhordump_mr" >/dev/null 2>&1

echo "[3/3] compare per-rank (OWNED entries, max|delta|=0)"
fail=0
for r in $(seq 0 $((NP-1))); do
    rs=$(printf '%05d' "$r")
    if python3 "$F3/tools/advhor_diff.py" "$RUN/advhor_f2.bin.$rs" "$RUN/advhor_f3.bin.$rs" >/dev/null 2>&1; then
        echo "  rank $rs: PASS (all fields max|delta|=0)"
    else
        echo "  rank $rs: FAIL"
        python3 "$F3/tools/advhor_diff.py" "$RUN/advhor_f2.bin.$rs" "$RUN/advhor_f3.bin.$rs" 2>&1 | grep -iE "FAIL|MISSING|mismatch" | head
        fail=1
    fi
done
if [[ $fail -eq 0 ]]; then
    echo "MULTI-RANK ADVECTION BYTE-GATE ($NP ranks): PASS — all ranks max|delta|=0"
else
    echo "MULTI-RANK ADVECTION BYTE-GATE ($NP ranks): FAIL"; exit 1
fi
