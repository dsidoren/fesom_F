#!/usr/bin/env bash
# End-to-end M2.10a forcing-READ operator byte-gate: dump FESOM2's u_wind/v_wind/Tair/
# shum/shortwave/longwave/prec_rain/prec_snow (1-rank pi, the REAL sbc_ini + sbc_do —
# netCDF read of the CORE2 stubs + spatial bilinear + linear time interp + g2r wind
# rotation + the update_atm_forcing mapping), dump FESOM3's (1-rank pi, transcribed
# mod_forcing_read on the identical CORE2 files + pinned model time), and compare for
# max|delta|=0 on all 8 fields. This is the FIRST netCDF I/O in the port.
#
# Prereqs (one-time, see docs/HANDOFF.md):
#   - FESOM2 rebuilt with src/fesom_forcing_dump.F90 (wired after forcing_setup in
#     fesom_module.F90) AND the np=1 fix in src/io_netcdf_workaround_module.F90
#     (build/lib64/libfesom.so; re-run cmake first since fesom_forcing_dump is NEW).
#   - pi mesh has the hand-crafted dist_1/ (shared with the other gates).
#   - FESOM3 built (with netCDF): build_intel_dp/bin/fesom_forcingdump.
set -euo pipefail
F3=/home/a/a270088/fesom3
RUN="${1:-/scratch/a/a270088/forcingdump_pi}"

echo "[1/3] FESOM2 oracle forcing dump (1-rank pi, REAL sbc_do)"
bash "$F3/tools/run_forcingdump_pi.sh" "$RUN" "$RUN/forcing_f2.bin" >/dev/null

echo "[2/3] FESOM3 forcing dump (1-rank pi, transcribed mod_forcing_read)"
source "$F3/env.sh" intel >/dev/null 2>&1
ulimit -s unlimited
export FESOM3_FORCING_OUT="$RUN/forcing_f3.bin"
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_forcingdump" >/dev/null 2>&1

echo "[3/3] compare"
python3 "$F3/tools/pressure_diff.py" "$RUN/forcing_f2.bin" "$RUN/forcing_f3.bin"
