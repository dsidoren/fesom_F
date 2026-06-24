#!/usr/bin/env bash
# M7d MULTI-RANK TKE + sw_pene FULLY-NATIVE forced lifecycle byte-gate (CORE2 dist_<NP>) — the
# linfs+TKE capstone at multi-rank. Runs the REAL FESOM2 forced lifecycle (use_ice + native CORE2
# forcing + ice EVP + oce_fluxes) with mix_scheme='cvmix_TKE' + use_sw_pene + Fer_GM + Redi at NP
# ranks, vs FESOM3's FULLY-NATIVE MULTI-RANK lifecycle with the same (FESOM3_MIX_TKE/SW_PENE/FER_GM/
# REDI). oce_mixing_tke.calc_cvmix_tke is optional-`partit` from the start (owned-loop over nNodO +
# exchange_nod(tke_Kv)/exchange_nod(tke_Av) BEFORE the element average; tke is NEVER exchanged) — so
# this is pure WIRING (the M4f/M5d lesson). Compares the per-rank gid-keyed 13 NODE substeps x 5 probes.
#
# ⚠️ multi-rank levante needs the env.sh KNEM flag — sourced by the gate runner.
#
#   tools/run_lifecycle_tke_native_gate_multirank.sh [np] [nsteps] [whichEVP] [run_dir]
set -euo pipefail
F3=/home/a/a270088/fesom3
NP="${1:-2}"
NSTEPS="${2:-3}"
WHICHEVP="${3:-0}"
RUN="${4:-/scratch/a/a270088/lifecycle_tke_native_mr${NP}}"

MIX_TKE=1 SW_PENE=1 FER_GM=1 REDI=1 \
    bash "$F3/tools/run_lifecycle_fullynative_gate_multirank.sh" "$NP" "$NSTEPS" "$WHICHEVP" "$RUN"
