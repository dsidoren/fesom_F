#!/usr/bin/env bash
# M7e zstar+TKE PRODUCTION byte-gate (CORE2 1-rank) — the FULL default config with TKE:
# which_ALE='zstar' (M6 free surface) + mix_scheme='cvmix_TKE' + GM + Redi + sw_pene + native
# ice/forcing. Expected PURE CONFIG flip (zstar done M6; TKE is an Av/Kv producer that does not
# touch the ALE thickness machinery — for zstar hnode/Z_3d_n are time-varying, read live by
# calc_cvmix_tke). Flags only, no new kernel.
#
#   tools/run_lifecycle_zstar_tke_native_gate_core2.sh [run_dir] [nsteps] [whichEVP]
set -euo pipefail
F3=/home/a/a270088/fesom3
RUN="${1:-/scratch/a/a270088/lifecycle_zstar_tke_native_core2}"
NSTEPS="${2:-3}"
WHICHEVP="${3:-0}"

WHICH_ALE=zstar MIX_TKE=1 SW_PENE=1 FER_GM=1 REDI=1 \
    bash "$F3/tools/run_lifecycle_fullynative_gate_core2.sh" "$RUN" "$NSTEPS" "$WHICHEVP"
