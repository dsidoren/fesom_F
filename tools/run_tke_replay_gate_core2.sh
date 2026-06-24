#!/usr/bin/env bash
# M7a-1 CONTROLLED-REPLAY byte-gate (CORE2 1-rank): lock the TKE column algebra
# (integrate_tke + solve_tridiag) in ISOLATION on bit-identical inputs.
#
#   [1] run the FESOM2 oracle (cvmix_TKE work_tke_dump) -> per-step tke_dump INPUT+OUTPUT files
#   [2] build + run fesom_tkereplay: ingest the oracle INPUT dumps -> integrate_tke -> emit OUTPUTS
#   [3] diff replay OUTPUTS vs oracle OUTPUTS, all steps, all tags -> max|Δ|=0
#
# Steps 2-3 (tke_old≠0) are the real prognostic-solver test; step 1 is the degenerate floor.
#
#   tools/run_tke_replay_gate_core2.sh [run_dir] [nsteps]
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
RUN="${1:-/scratch/a/a270088/tke_replay_core2}"
NSTEPS="${2:-3}"
ORACLE="$RUN/tke_oracle"
REPLAY="$RUN/tke_replay"

echo "[1/3] FESOM2 oracle tke_dump ($NSTEPS steps, 1-rank CORE2, cvmix_TKE linfs)"
NSTEPS="$NSTEPS" bash "$F3/tools/run_tkedump_oracle.sh" "$RUN/oracle_run" "$ORACLE" | tail -4

echo "[2/3] build + run fesom_tkereplay (ingest oracle inputs -> integrate_tke -> emit outputs)"
source "$F3/env.sh" intel >/dev/null 2>&1
cmake --build "$F3/build_intel_dp" --target fesom_tkereplay -j 8 >/dev/null 2>&1
mkdir -p "$REPLAY"
export FESOM3_MESH_DIR="$COREMESH"
export FESOM3_TKE_IN_DIR="$ORACLE" FESOM3_TKE_OUT_DIR="$REPLAY"
export FESOM_TKE_DUMP_STEPS="$NSTEPS"
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_tkereplay" 2>&1 | grep -E 'nod2D|replayed|cannot|short'

echo "[3/3] compare replay vs oracle (max|Δ|=0; tke/tkeav/tkekv + lmix/pr/tbpr/tspr/tdif/tdis ...)"
# The replay emits the oracle's EXACT es24.16 format, so byte-exact values => byte-identical
# files (cmp, instant). The python differ (tke_dump_diff.py) localizes the worst gid/comp/Part
# if a file differs.
fail=0; pass=0
for s in $(seq 1 "$NSTEPS"); do
  for tag in tke tkeav tkekv lmix pr tbpr tspr tdif tdis twin tiwf tbck ttot; do
    f="tke_dump_s${s}_${tag}_rank0.txt"
    if cmp -s "$REPLAY/$f" "$ORACLE/$f"; then pass=$((pass+1))
    else
      fail=$((fail+1)); echo "  DIFFER: s$s $tag"
      python3 "$F3/tools/tke_dump_diff.py" "$REPLAY" "$ORACLE" --steps "$s" --tags "$tag" | sed -n '1p'
    fi
  done
done
echo "----------------------------------------------------------------------"
if [ "$fail" -eq 0 ]; then
    echo "M7a-1 replay gate: PASS — $pass/$pass tags byte-identical (max|Δ|=0), all $NSTEPS steps"
else
    echo "M7a-1 replay gate: FAIL — $fail tag(s) differ (see localization above)"; exit 1
fi
