#!/usr/bin/env bash
# M4e GM+Redi FULLY-NATIVE forced lifecycle byte-gate (CORE2 1-rank) — the production tracer
# physics integration. Runs the REAL FESOM2 forced lifecycle (use_ice + native CORE2 forcing +
# ice EVP + oce_fluxes) with the work_core Gent-McWilliams bolus advection AND Redi isopycnal
# diffusion KEPT ON (FER_GM=1 REDI=1), vs FESOM3's FULLY-NATIVE lifecycle with GM+Redi enabled
# (FESOM3_FER_GM=1 FESOM3_REDI=1). The whole air-sea forcing is computed natively (8 NCAR fields
# + NCAR bulk + runoff/Ssurf), the sea-ice + flux coupling is native, AND the GM+Redi tracer
# physics runs inside step_oce on the forced state. Compares the 13 NODE substeps x 5 probes over
# N steps (--ignore-substep=2). This is the M4 capstone: GM bolus + Redi diffusion + sea ice +
# native forcing all live in one byte-exact coupled run.
#
# The GM/Redi-OFF no-regression baseline is run_lifecycle_fullynative_gate_core2.sh; the UNFORCED
# GM-only / GM+Redi gates are run_lifecycle_gm_gate_core2.sh / run_lifecycle_redi_gate_core2.sh.
#
#   tools/run_lifecycle_gmredi_native_gate_core2.sh [run_dir] [nsteps] [whichEVP]
set -euo pipefail
F3=/home/a/a270088/fesom3
RUN="${1:-/scratch/a/a270088/lifecycle_gmredi_native_core2}"
NSTEPS="${2:-3}"
WHICHEVP="${3:-0}"

FER_GM=1 REDI=1 bash "$F3/tools/run_lifecycle_fullynative_gate_core2.sh" "$RUN" "$NSTEPS" "$WHICHEVP"
