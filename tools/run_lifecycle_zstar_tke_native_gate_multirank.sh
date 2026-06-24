#!/usr/bin/env bash
# M7e zstar+TKE PRODUCTION byte-gate (CORE2 dist_<NP>) — the FULL default config with TKE at
# MULTI-RANK: which_ALE='zstar' + mix_scheme='cvmix_TKE' + GM + Redi + sw_pene + native ice/forcing.
# Pure WIRING on top of the M6a-4 zstar MR + the M7d TKE MR block (calc_cvmix_tke is optional-partit;
# tke is partition-local). Flags only, no new kernel.
#
#   tools/run_lifecycle_zstar_tke_native_gate_multirank.sh [np] [nsteps] [whichEVP] [run_dir]
set -euo pipefail
F3=/home/a/a270088/fesom3
NP="${1:-2}"
NSTEPS="${2:-3}"
WHICHEVP="${3:-0}"
RUN="${4:-/scratch/a/a270088/lifecycle_zstar_tke_native_mr${NP}}"

WHICH_ALE=zstar MIX_TKE=1 SW_PENE=1 FER_GM=1 REDI=1 \
    bash "$F3/tools/run_lifecycle_fullynative_gate_multirank.sh" "$NP" "$NSTEPS" "$WHICHEVP" "$RUN"
