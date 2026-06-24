#!/usr/bin/env bash
# M4f MULTI-RANK GM+Redi FULLY-NATIVE forced lifecycle byte-gate (CORE2 dist_<NP>) — the M4
# capstone at multi-rank. Runs the REAL FESOM2 forced lifecycle (use_ice + native CORE2 forcing +
# ice EVP + oce_fluxes) with work_core Gent-McWilliams bolus advection AND Redi isopycnal diffusion
# KEPT ON (FER_GM=1 REDI=1) at NP ranks, vs FESOM3's FULLY-NATIVE MULTI-RANK lifecycle with GM+Redi
# enabled (FESOM3_FER_GM=1 FESOM3_REDI=1). The GM+Redi tracer physics runs inside step_oce through
# the optional partit — every M4 producer/GM/Redi routine does owned-loop bounds + the FESOM2 halo
# exchanges. Compares the per-rank gid-keyed 13 NODE substeps x 5 probes over N steps.
#
# ⚠️ multi-rank levante needs the env.sh KNEM flag (L35) — sourced by the gate runner.
# The GM/Redi-OFF MR no-regression baseline is run_lifecycle_fullynative_gate_multirank.sh.
#
#   tools/run_lifecycle_gmredi_native_gate_multirank.sh [np] [nsteps] [whichEVP] [run_dir]
set -euo pipefail
F3=/home/a/a270088/fesom3
NP="${1:-2}"
NSTEPS="${2:-3}"
WHICHEVP="${3:-0}"
RUN="${4:-/scratch/a/a270088/lifecycle_gmredi_native_mr${NP}}"

FER_GM=1 REDI=1 bash "$F3/tools/run_lifecycle_fullynative_gate_multirank.sh" "$NP" "$NSTEPS" "$WHICHEVP" "$RUN"
