#!/usr/bin/env bash
# M7c TKE + sw_pene FULLY-NATIVE forced lifecycle byte-gate (CORE2 1-rank) — the reduced
# PRODUCTION linfs+TKE gate. Activates the surface WIND forcing (forc_tke_surf=|stress_node_surf|/
# density_0 ≠ 0 -> the Neumann (cd*forc_tke_surf**(3./2.))/dzt(1) term), with GM/Redi on (work_*_tke
# config). ⚠️ work_*_tke has use_sw_pene=.true., so SW_PENE=1 is set ALONGSIDE MIX_TKE=1 (else the
# shortwave heating term is silently absent and the tracer column diverges).
#
#   tools/run_lifecycle_tke_native_gate_core2.sh [run_dir] [nsteps] [whichEVP]
set -euo pipefail
F3=/home/a/a270088/fesom3
RUN="${1:-/scratch/a/a270088/lifecycle_tke_native_core2}"
NSTEPS="${2:-3}"
WHICHEVP="${3:-0}"

MIX_TKE=1 SW_PENE=1 FER_GM=1 REDI=1 \
    bash "$F3/tools/run_lifecycle_fullynative_gate_core2.sh" "$RUN" "$NSTEPS" "$WHICHEVP"
