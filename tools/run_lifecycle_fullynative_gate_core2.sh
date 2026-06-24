#!/usr/bin/env bash
# M3f-3 (COMPLETE) FULLY-NATIVE forced lifecycle byte-gate (CORE2 1-rank). The REAL forced
# FESOM2 lifecycle (use_ice=.true., real CORE2 forcing + ice EVP + oce_fluxes) vs FESOM3's
# FULLY NATIVE lifecycle: FESOM3 reads NOTHING from the oracle atmosphere dump — it computes
# the WHOLE air-sea forcing itself each step:
#   8 NCAR fields  (mod_forcing_read: shortwave/longwave/Tair/shum/prec_rain/prec_snow/wind)
#   NCAR bulk      (mod_forcing_bulk: Ch/Ce_atm_oce + stress_atmoce + wind-on-ice stress)
#   runoff + Ssurf (mod_forcing_other: read_other_NetCDF monthly climatology, M3f-3b)
# then ocean2ice -> ice_timestep -> oce_fluxes_mom -> oce_fluxes -> step_oce. NO FESOM3_ATMFLUX_FILE.
# Compares the 13 NODE substeps x 5 probes over N steps (--ignore-substep=2), like the M2.11c-2
# prescribed-flux gate. A per-step native-vs-oracle flux self-check (FESOM3_FLUX_FILE) localizes
# any coupling drift before the substep gate.
#
#   tools/run_lifecycle_fullynative_gate_core2.sh [run_dir] [nsteps] [whichEVP]
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
RUN="${1:-/scratch/a/a270088/lifecycle_fullynative_core2}"
NSTEPS="${2:-3}"
WHICHEVP="${3:-0}"
# M4e: FER_GM=1 / REDI=1 turn ON the work_core GM bolus / Redi isopycnal diffusion in BOTH the
# oracle and the FESOM3 native lifecycle. Default 0 = the proven M3f GM/Redi-off no-regression gate.
FER_GM="${FER_GM:-0}"
REDI="${REDI:-0}"

echo "[1/3] FESOM2 oracle FORCED lifecycle ($NSTEPS steps, 1-rank CORE2, use_ice; FER_GM=$FER_GM REDI=$REDI)"
FER_GM="$FER_GM" REDI="$REDI" \
    bash "$F3/tools/run_lifecycle_forced_core2.sh" "$RUN" "$RUN/lifef_f2" "$RUN/flux_f2" "$NSTEPS" "$RUN/atmflux_f2" "$WHICHEVP" | tail -4

echo "[2/3] FESOM3 FULLY NATIVE lifecycle ($NSTEPS steps, whichEVP=$WHICHEVP) — NO prescribed atmosphere"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE"
# NO FESOM3_ATMFLUX_FILE -> fully native; the whole atmosphere is computed in-driver.
export FESOM3_FORCING_DIR="/home/a/a270088/port2/fesom2/test/input/global"  # native NCAR read
export FESOM3_RUNOFF_FILE="$POOL/CORE2_runoff.nc"   # native runoff (M3f-3b)
export FESOM3_SSS_FILE="$POOL/PHC2_salx.nc"          # native SSS restoring (M3f-3b)
export FESOM3_FLUX_FILE="$RUN/flux_f2.00000"        # native-vs-oracle flux self-check
export FESOM3_WHICHEVP="$WHICHEVP"
if [ "$FER_GM" = 1 ]; then export FESOM3_FER_GM=1; fi   # M4e: GM bolus in the native lifecycle
if [ "$REDI" = 1 ]; then export FESOM3_REDI=1; fi       # M4e: Redi isopycnal diffusion
export FESOM_DUMP_FILE="$RUN/lifen_f3" FESOM_DUMP_MAXSTEPS="$NSTEPS" FESOM3_NSTEPS="$NSTEPS"
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_lifecycle_native" 2>&1 | grep -E 'nod2D|NATIVE|native runoff|ENABLED|selfcheck|step [0-9]|done'

echo "[3/3] compare (13 NODE substeps x 5 probes x $NSTEPS steps; SW_AB substep 2 ignored)"
python3 "$F3/tools/dump_diff.py" "$RUN/lifef_f2.00000" "$RUN/lifen_f3.00000" --ignore-substep=2
