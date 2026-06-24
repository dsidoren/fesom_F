#!/usr/bin/env bash
# M5d MULTI-RANK KPP + sw_pene FULLY-NATIVE forced lifecycle byte-gate (CORE2 dist_<NP>) — the M5
# capstone at multi-rank. Runs the REAL FESOM2 forced lifecycle (use_ice + native CORE2 forcing +
# ice EVP + oce_fluxes) with the FULL work_core vertical physics KEPT ON (mix_scheme='KPP' +
# use_sw_pene + Fer_GM + Redi) at NP ranks, vs FESOM3's FULLY-NATIVE MULTI-RANK lifecycle with the
# same four ON (FESOM3_MIX_KPP/SW_PENE/FER_GM/REDI). The KPP module is optional-`partit` from the
# start (owned-loop bounds + the oracle's exchange_nod(blmc/diffK/ghats/viscA) + smooth_blmc), the
# sw_3d/ghats tracer terms loop owned nodes, and cal_shortwave_rad loops owned+halo via partit — so
# this is pure WIRING (the M4f lesson). Compares the per-rank gid-keyed 13 NODE substeps x 5 probes.
#
# ⚠️ multi-rank levante needs the env.sh KNEM flag (L35) — sourced by the gate runner.
# The 1-rank production gate is run_lifecycle_kpp_native_gate_core2.sh; the GM/Redi (PP, no-KPP) MR
# gate is run_lifecycle_gmredi_native_gate_multirank.sh.
#
#   tools/run_lifecycle_kpp_native_gate_multirank.sh [np] [nsteps] [whichEVP] [run_dir]
set -euo pipefail
F3=/home/a/a270088/fesom3
NP="${1:-2}"
NSTEPS="${2:-3}"
WHICHEVP="${3:-0}"
RUN="${4:-/scratch/a/a270088/lifecycle_kpp_native_mr${NP}}"

MIX_KPP=1 SW_PENE=1 FER_GM=1 REDI=1 \
    bash "$F3/tools/run_lifecycle_fullynative_gate_multirank.sh" "$NP" "$NSTEPS" "$WHICHEVP" "$RUN"
