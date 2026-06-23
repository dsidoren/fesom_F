#!/usr/bin/env bash
# M3d sea-ice THERMODYNAMICS byte-gate: run the FESOM2 oracle's REAL ocean2ice + EVP solve
# (120 subcycles) + ice FCT advection (ice_TG_rhs + ice_fct_solve) + cut_off + thermodynamics
# on the cold-start ice tracers + a prescribed analytic ocean+atmospheric forcing (1-rank
# CORE2); run FESOM3's transcribed mod_ice_thermo::cut_off + thermodynamics (1-rank CORE2) on
# the IDENTICAL prescribed state + the M3b/M3c-proven EVP velocity + FCT advection; compare
# for max|delta|=0. Both build the ice IC from the SAME do_ic3d surface T + the SAME ssh_stiff
# CSR mass matrix, prescribe the SAME analytic forcing from the byte-identical rotated coords,
# and use the SAME &ice_therm namelist doubles + atmflux scalar config.
# Tests BOTH whichEVP=0 (standard EVP) AND whichEVP=1 (mEVP) feeding the advection upstream.
# Fields: ice_a_ice/m_ice/m_snow (post-thermo) + ice_flx_h/flx_fw + ice_t_skin.
#
# Prereqs:
#   - FESOM2 build/lib/libfesom.so has the ice_thermo dump shim (fesom_module.F90 ->
#     ice_thermo_dump_write); rebuild after editing src/fesom_ice_dump.F90 + fesom_module.F90.
#   - CORE2 mesh has the hand-crafted dist_1/ (tools/make_dist1.py, from M2.11a).
#   - FESOM3 built CLEAN (configure.sh): build_intel_dp/bin/fesom_icethermodump.
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
RUN="${1:-/scratch/a/a270088/icethermodump_core2}"

source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH"
export FESOM3_IC_FILE="$ICFILE"
ulimit -s unlimited

run_one () {
    local mode="$1" name="$2"
    echo "============================================================"
    echo "  M3d ice THERMO gate — whichEVP=${mode} (${name})"
    echo "============================================================"
    echo "[1/3] FESOM2 oracle thermo dump (1-rank CORE2, use_ice, reduced-M2)"
    bash "$F3/tools/run_icethermodump_core2.sh" "$RUN" "$RUN/icethermo_f2_${mode}.bin" "$mode"

    echo "[2/3] FESOM3 thermo dump (1-rank CORE2, whichEVP=${mode})"
    export FESOM3_WHICHEVP="$mode"
    export FESOM3_THERMO_OUT="$RUN/icethermo_f3_${mode}.bin"
    mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
        "$F3/build_intel_dp/bin/fesom_icethermodump" 2>&1 | grep -E 'nod2D|range|whichEVP|post-adv|post-thermo|wrote'

    echo "[3/3] compare (whichEVP=${mode})"
    python3 "$F3/tools/pressure_diff.py" "$RUN/icethermo_f2_${mode}.bin" "$RUN/icethermo_f3_${mode}.bin"
}

run_one 0 "standard EVP"
run_one 1 "modified EVP (mEVP)"
