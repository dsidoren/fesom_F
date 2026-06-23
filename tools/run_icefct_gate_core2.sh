#!/usr/bin/env bash
# M3c sea-ice FCT-advection byte-gate: run the FESOM2 oracle's REAL ocean2ice + EVP solve
# (120 subcycles) + ice FCT advection (ice_TG_rhs + ice_fct_solve) on the cold-start ice
# tracers + a prescribed analytic ocean-forcing state (surface UV / hbar / stress_atmice),
# dump the TG rhs (rhs_a/m/ms) + post-advection a_ice/m_ice/m_snow (1-rank CORE2); run
# FESOM3's transcribed mod_ice_fct::ice_TG_rhs + ice_fct_solve (1-rank CORE2) on the
# IDENTICAL prescribed state + the M3b-proven EVP velocity; compare for max|delta|=0. Both
# build the ice IC from the SAME do_ic3d surface T + the SAME ssh_stiff CSR mass matrix, and
# prescribe the SAME analytic forcing from the byte-identical rotated coords.
# Tests BOTH whichEVP=0 (standard EVP) AND whichEVP=1 (modified EVP / mEVP) feeding the advection.
# Fields: ice_rhs_a/m/ms (ice_TG_rhs) + ice_a_ice/m_ice/m_snow (post ice_fct_solve).
#
# Prereqs:
#   - FESOM2 build/lib/libfesom.so has the ice_fct dump shim (fesom_module.F90 ->
#     ice_fct_dump_write); rebuild after editing src/fesom_ice_dump.F90 (re-run cmake only if
#     a NEW shim file was added — here it is the SAME file, so plain make suffices).
#   - CORE2 mesh has the hand-crafted dist_1/ (tools/make_dist1.py, from M2.11a).
#   - FESOM3 built CLEAN (configure.sh): build_intel_dp/bin/fesom_icefctdump.
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
RUN="${1:-/scratch/a/a270088/icefctdump_core2}"

source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH"
export FESOM3_IC_FILE="$ICFILE"
ulimit -s unlimited

run_one () {
    local mode="$1" name="$2"
    echo "============================================================"
    echo "  M3c ice FCT gate — whichEVP=${mode} (${name})"
    echo "============================================================"
    echo "[1/3] FESOM2 oracle FCT dump (1-rank CORE2, use_ice, reduced-M2)"
    bash "$F3/tools/run_icefctdump_core2.sh" "$RUN" "$RUN/icefct_f2_${mode}.bin" "$mode"

    echo "[2/3] FESOM3 FCT dump (1-rank CORE2, whichEVP=${mode})"
    export FESOM3_WHICHEVP="$mode"
    export FESOM3_FCT_OUT="$RUN/icefct_f3_${mode}.bin"
    mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
        "$F3/build_intel_dp/bin/fesom_icefctdump" 2>&1 | grep -E 'nod2D|range|whichEVP|max\|uice|post-adv|wrote'

    echo "[3/3] compare (whichEVP=${mode})"
    python3 "$F3/tools/pressure_diff.py" "$RUN/icefct_f2_${mode}.bin" "$RUN/icefct_f3_${mode}.bin"
}

run_one 0 "standard EVP"
run_one 1 "modified EVP (mEVP)"
