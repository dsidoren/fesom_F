#!/usr/bin/env bash
# M6a-3 forced/native zstar byte-gate (CORE2 1-rank): the FULLY-NATIVE forced lifecycle
# (real CORE2 forcing + ice EVP + oce_fluxes) with which_ALE='zstar' instead of 'linfs'.
# This is the FIRST test of the REAL freshwater-flux path: use_virt_salt=.false. (is_nonlinfs=1)
# => no virtual salt; the water_flux enters the SSH/thickness (compute_hbar_ale + vert_vel_ale
# Wvel-=water_flux), the ice-thermo levitating split feeds water_flux (oce_fluxes), and the
# salinity surface BC uses real_salt_flux. Reduced column physics (PP, no GM/Redi, no sw_pene)
# by default to ISOLATE the zstar + freshwater changes; add FER_GM=1 REDI=1 MIX_KPP=1 SW_PENE=1
# for the full production stack.
#
#   tools/run_lifecycle_zstar_native_gate_core2.sh [run_dir] [nsteps] [whichEVP]
#   (env passthrough: FER_GM REDI MIX_KPP SW_PENE)
set -euo pipefail
F3=/home/a/a270088/fesom3
RUN="${1:-/scratch/a/a270088/lifecycle_zstar_native_core2}"
NSTEPS="${2:-3}"
WHICHEVP="${3:-0}"
WHICH_ALE=zstar bash "$F3/tools/run_lifecycle_fullynative_gate_core2.sh" "$RUN" "$NSTEPS" "$WHICHEVP"
