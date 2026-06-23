#!/usr/bin/env bash
# M3e sea-ice -> ocean COUPLING-OUT byte-gate (THE PAYOFF): run the FESOM2 oracle's REAL
# ocean2ice + EVP solve (120 subcycles) + ice FCT advection + cut_off + thermodynamics +
# oce_fluxes_mom + oce_fluxes on the cold-start ice + a prescribed analytic ocean/atmosphere/
# atm-ocean-stress/SSS-climatology (1-rank CORE2); run FESOM3's transcribed
# mod_ice_oce_coupling::oce_fluxes_mom + oce_fluxes (1-rank CORE2) on the IDENTICAL prescribed
# state + the M3b/c/d-proven EVP + FCT + thermo; compare for max|delta|=0. Both build the ice
# IC from the SAME do_ic3d surface T/S + the SAME ssh_stiff CSR mass matrix, prescribe the SAME
# analytic forcing from the byte-identical rotated coords, and use the SAME reduced-M2/linfs +
# &ice_therm namelist doubles + surf_relax_S=1.929e-06. Tests BOTH whichEVP=0 (standard EVP)
# AND whichEVP=1 (mEVP) feeding the advection upstream.
# Fields: heat_flux/water_flux/virtual_salt/relax_salt (nod2D) + stress_surf (2,elem2D) —
# the proven M2.11c-2 fesom_flux_dump field set, now produced NATIVELY by FESOM3.
#
# Prereqs:
#   - FESOM2 build/lib/libfesom.so has the M3e flux dump shim (fesom_module.F90 ->
#     ice_flux_dump_write); rebuild after editing src/fesom_ice_dump.F90 + fesom_module.F90.
#   - CORE2 mesh has the hand-crafted dist_1/ (tools/make_dist1.py, from M2.11a).
#   - FESOM3 built CLEAN (configure.sh): build_intel_dp/bin/fesom_icefluxdump.
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
RUN="${1:-/scratch/a/a270088/icefluxdump_core2}"

source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH"
export FESOM3_IC_FILE="$ICFILE"
ulimit -s unlimited

run_one () {
    local mode="$1" name="$2"
    echo "============================================================"
    echo "  M3e ice -> ocean FLUX gate — whichEVP=${mode} (${name})"
    echo "============================================================"
    echo "[1/3] FESOM2 oracle flux dump (1-rank CORE2, use_ice, reduced-M2)"
    bash "$F3/tools/run_icefluxdump_core2.sh" "$RUN" "$RUN/iceflux_f2_${mode}.bin" "$mode"

    echo "[2/3] FESOM3 flux dump (1-rank CORE2, whichEVP=${mode})"
    export FESOM3_WHICHEVP="$mode"
    export FESOM3_FLUX_OUT="$RUN/iceflux_f3_${mode}.bin"
    mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
        "$F3/build_intel_dp/bin/fesom_icefluxdump" 2>&1 | grep -E 'nod2D|range|whichEVP|post-thermo|max\||wrote'

    echo "[3/3] compare (whichEVP=${mode})"
    python3 "$F3/tools/pressure_diff.py" "$RUN/iceflux_f2_${mode}.bin" "$RUN/iceflux_f3_${mode}.bin"
}

run_one 0 "standard EVP"
run_one 1 "modified EVP (mEVP)"
