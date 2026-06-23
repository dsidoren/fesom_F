#!/usr/bin/env bash
# M3b sea-ice EVP-dynamics byte-gate: run the FESOM2 oracle's REAL ocean2ice + EVP solve
# (120 subcycles) on the cold-start ice IC + a prescribed analytic ocean-forcing state
# (surface UV / hbar / stress_atmice), dump srfoce_u/v + uice/vice + sigma11/12/22 (1-rank
# CORE2); run FESOM3's transcribed mod_ice_dyn::ocean2ice + EVPdynamics_solve (1-rank
# CORE2) on the IDENTICAL prescribed state; compare for max|delta|=0. Both build the ice IC
# from the SAME do_ic3d surface T + the SAME ssh_stiff CSR mass matrix, and prescribe the
# SAME analytic forcing from the byte-identical rotated coords.
# Tests BOTH whichEVP=0 (standard EVP) AND whichEVP=1 (modified EVP / mEVP).
# Fields: ice_srfoce_u/v (ocean2ice), ice_uice/vice, ice_sigma11/12/22 (EVP).
#
# Prereqs:
#   - FESOM2 build/lib/libfesom.so has the evp dump shim (fesom_module.F90 -> evp_dump_write);
#     rebuild after editing src/fesom_ice_dump.F90 (re-run cmake if a NEW shim file was added).
#   - CORE2 mesh has the hand-crafted dist_1/ (tools/make_dist1.py, from M2.11a).
#   - FESOM3 built CLEAN (configure.sh): build_intel_dp/bin/fesom_evpdump.
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
RUN="${1:-/scratch/a/a270088/evpdump_core2}"

source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH"
export FESOM3_IC_FILE="$ICFILE"
ulimit -s unlimited

run_one () {
    local mode="$1" name="$2"
    echo "============================================================"
    echo "  M3b EVP gate — whichEVP=${mode} (${name})"
    echo "============================================================"
    echo "[1/3] FESOM2 oracle EVP dump (1-rank CORE2, use_ice, reduced-M2)"
    bash "$F3/tools/run_evpdump_core2.sh" "$RUN" "$RUN/evp_f2_${mode}.bin" "$mode"

    echo "[2/3] FESOM3 EVP dump (1-rank CORE2, whichEVP=${mode})"
    export FESOM3_WHICHEVP="$mode"
    export FESOM3_EVP_OUT="$RUN/evp_f3_${mode}.bin"
    mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
        "$F3/build_intel_dp/bin/fesom_evpdump" 2>&1 | grep -E 'nod2D|range|whichEVP|max\|uice|max\|sig|wrote'

    echo "[3/3] compare (whichEVP=${mode})"
    python3 "$F3/tools/pressure_diff.py" "$RUN/evp_f2_${mode}.bin" "$RUN/evp_f3_${mode}.bin"
}

run_one 0 "standard EVP"
run_one 1 "modified EVP (mEVP)"
