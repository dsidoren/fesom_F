#!/usr/bin/env bash
# M2.11b CORE2 initial-conditions (do_ic3d) byte-gate: dump the FESOM2 oracle's live
# Tclim/Sclim (1-rank CORE2, the REAL oce_initial_state), dump FESOM3's transcribed
# oce_initial_state::do_ic3d (1-rank CORE2), compare for max|delta|=0. Both read the
# SAME phc3.0_winter.nc climatology; the gate covers the netCDF 3-D read + bilinear +
# vertical-interp + extrapolation + insitu2pot chain. Fields: Z_3d_n (input), ic_temp
# (potential T), ic_salt.
#
# Prereqs (one-time):
#   - FESOM2 build/lib64/libfesom.so has the IC dump shim (oce_setup_step.F90 ->
#     ic_dump_write); rebuild after adding src/fesom_ic_dump.F90 (re-run cmake first
#     so the GLOB picks up the new file).
#   - CORE2 mesh has the hand-crafted dist_1/ (tools/make_dist1.py, from M2.11a).
#   - FESOM3 built CLEAN: build_intel_dp/bin/fesom_icdump.
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
RUN="${1:-/scratch/a/a270088/icdump_core2}"

echo "[1/3] FESOM2 oracle IC dump (1-rank CORE2, reduced-M2 linfs)"
bash "$F3/tools/run_icdump_core2.sh" "$RUN" "$RUN/ic_f2.bin"

echo "[2/3] FESOM3 IC dump (1-rank CORE2)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH"
export FESOM3_IC_FILE="$ICFILE"
export FESOM3_IC_OUT="$RUN/ic_f3.bin"
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_icdump" 2>&1 | grep -E 'nod2D|range|wrote'

echo "[3/3] compare"
python3 "$F3/tools/pressure_diff.py" "$RUN/ic_f2.bin" "$RUN/ic_f3.bin"
