#!/usr/bin/env bash
# M3f native-flux forced lifecycle byte-gate (CORE2 1-rank). The REAL forced FESOM2
# lifecycle (use_ice=.true., real CORE2 forcing + ice EVP + oce_fluxes) vs FESOM3's
# NATIVE-flux lifecycle: FESOM3 reads only the per-step POST-bulk atmospheric forcing
# (fesom_atmflux_dump: shortwave/longwave/Tair/shum/prec/runoff/wind/Ch-Ce/stress_atm*/
# Ssurf) and computes heat_flux/water_flux/virtual_salt/relax_salt/stress_surf ITSELF via
# ocean2ice + ice_timestep + oce_fluxes_mom + oce_fluxes, then steps the ocean on them.
# Compares the 13 NODE substeps x 5 probes over N steps (--ignore-substep=2), like the
# M2.11c-2 prescribed-flux gate. A per-step native-vs-oracle flux self-check (FESOM3_FLUX_FILE)
# localizes any coupling drift before the substep gate.
#
#   tools/run_lifecycle_native_gate_core2.sh [run_dir] [nsteps] [whichEVP]
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
RUN="${1:-/scratch/a/a270088/lifecycle_native_core2}"
NSTEPS="${2:-3}"
WHICHEVP="${3:-0}"

echo "[1/3] FESOM2 oracle FORCED lifecycle ($NSTEPS steps, 1-rank CORE2, use_ice) + atm-forcing dump"
bash "$F3/tools/run_lifecycle_forced_core2.sh" "$RUN" "$RUN/lifef_f2" "$RUN/flux_f2" "$NSTEPS" "$RUN/atmflux_f2" "$WHICHEVP" | tail -4

echo "[2/3] FESOM3 native-flux lifecycle ($NSTEPS steps, whichEVP=$WHICHEVP) — fluxes computed natively"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE"
export FESOM3_ATMFLUX_FILE="$RUN/atmflux_f2.00000"
export FESOM3_FLUX_FILE="$RUN/flux_f2.00000"        # native-vs-oracle flux self-check
export FESOM3_FORCING_DIR="/home/a/a270088/port2/fesom2/test/input/global"  # M3f-3 native forcing self-check
export FESOM3_WHICHEVP="$WHICHEVP"
export FESOM_DUMP_FILE="$RUN/lifen_f3" FESOM_DUMP_MAXSTEPS="$NSTEPS" FESOM3_NSTEPS="$NSTEPS"
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_lifecycle_native" 2>&1 | grep -E 'nod2D|prescribing|selfcheck|forcing step|NATIVE|step [0-9]|done'

echo "[3/3] compare (13 NODE substeps x 5 probes x $NSTEPS steps; SW_AB substep 2 ignored)"
python3 "$F3/tools/dump_diff.py" "$RUN/lifef_f2.00000" "$RUN/lifen_f3.00000" --ignore-substep=2
