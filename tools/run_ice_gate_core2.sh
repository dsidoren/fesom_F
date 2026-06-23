#!/usr/bin/env bash
# M3a sea-ice foundation byte-gate: dump the FESOM2 oracle's live cold-start ice IC
# (a_ice/m_ice/m_snow) + FCT mass matrix (1-rank CORE2, the REAL ice_setup), dump
# FESOM3's transcribed mod_ice_setup::ice_setup (1-rank CORE2), compare for max|delta|=0.
# Both build the ice IC from the SAME do_ic3d surface T (SST<0 sign test) and the SAME
# ssh_stiff CSR (mass matrix). Fields: ice_a_ice, ice_m_ice, ice_m_snow, ice_massmatrix.
#
# Prereqs:
#   - FESOM2 build/lib/libfesom.so has the ice dump shim (fesom_module.F90 -> ice_dump_write);
#     rebuild after adding src/fesom_ice_dump.F90 (re-run cmake first so the GLOB picks it up).
#   - CORE2 mesh has the hand-crafted dist_1/ (tools/make_dist1.py, from M2.11a).
#   - FESOM3 built CLEAN (configure.sh): build_intel_dp/bin/fesom_icedump.
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
RUN="${1:-/scratch/a/a270088/icedump_core2}"

echo "[1/3] FESOM2 oracle ice dump (1-rank CORE2, use_ice, reduced-M2)"
bash "$F3/tools/run_icedump_core2.sh" "$RUN" "$RUN/ice_f2.bin"

echo "[2/3] FESOM3 ice dump (1-rank CORE2)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH"
export FESOM3_IC_FILE="$ICFILE"
export FESOM3_ICE_OUT="$RUN/ice_f3.bin"
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_icedump" 2>&1 | grep -E 'nod2D|range|max a_ice|nza|wrote'

echo "[3/3] compare"
python3 "$F3/tools/pressure_diff.py" "$RUN/ice_f2.bin" "$RUN/ice_f3.bin"
