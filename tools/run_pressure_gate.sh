#!/usr/bin/env bash
# End-to-end M2.1 EOS/pressure/N^2 operator byte-gate: dump FESOM2's density_m_rho0/
# hpressure/bvfreq (1-rank pi, real pressure_bv on a prescribed T/S), dump FESOM3's
# (1-rank pi, transcribed oce_pressure_bv on the identical prescribed T/S), and
# compare for max|delta|=0 on density_m_rho0 + hpressure + bvfreq (raw + smoothed)
# and every input (temp, salt, density_ref, zbar_3d_n, Z_3d_n, hnode).
#
# Prereqs (one-time, see docs/HANDOFF.md):
#   - FESOM2 rebuilt with src/fesom_pressure_dump.F90 wired into ocean_setup
#     (build/lib64/libfesom.so; fesom.x loads it at runtime).
#   - pi mesh has the hand-crafted dist_1/ (shared with the geometry/advection gates).
#   - FESOM3 built: build_intel_dp/bin/fesom_pressuredump.
set -euo pipefail
F3=/home/a/a270088/fesom3
RUN="${1:-/scratch/a/a270088/pressuredump_pi}"

echo "[1/3] FESOM2 oracle pressure dump (1-rank pi)"
bash "$F3/tools/run_pressuredump_pi.sh" "$RUN" "$RUN/pressure_f2.bin" >/dev/null

echo "[2/3] FESOM3 pressure dump (1-rank pi)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_PRESSURE_OUT="$RUN/pressure_f3.bin"
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_pressuredump" >/dev/null 2>&1

echo "[3/3] compare"
python3 "$F3/tools/pressure_diff.py" "$RUN/pressure_f2.bin" "$RUN/pressure_f3.bin"
