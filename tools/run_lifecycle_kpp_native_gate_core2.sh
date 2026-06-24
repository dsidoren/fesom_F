#!/usr/bin/env bash
# M5c KPP + sw_pene FULLY-NATIVE forced lifecycle byte-gate (CORE2 1-rank) — THE production
# CORE2 tracer-physics gate. Runs the REAL FESOM2 forced lifecycle (use_ice + native CORE2
# forcing + ice EVP + oce_fluxes) with the FULL work_core vertical physics KEPT ON:
#   mix_scheme='KPP'  (MIX_KPP=1)  — the K-Profile-Parameterization boundary-layer scheme
#   use_sw_pene=.true.(SW_PENE=1)  — shortwave penetration (cal_shortwave_rad -> sw_3d tracer term)
#   Fer_GM=.true.     (FER_GM=1)   — Gent-McWilliams bolus advection (M4)
#   Redi=.true.       (REDI=1)     — Redi isopycnal diffusion (M4)
# vs FESOM3's FULLY-NATIVE lifecycle with the same four ON (FESOM3_MIX_KPP/SW_PENE/FER_GM/REDI).
# This is the production CORE2 column physics: KPP Kv/Av + mo_convect, GM bolus + Redi, AND the
# sw_3d 3D shortwave heating in the tracer TDMA, all on the natively-forced coupled state.
# Compares the 13 NODE substeps x 5 probes over N steps (--ignore-substep=2).
#
# The KPP-only (no GM/Redi/sw_pene) UNFORCED gate is run_lifecycle_kpp_gate_core2.sh; the GM+Redi
# (PP, no-KPP) native gate is run_lifecycle_gmredi_native_gate_core2.sh. The ghats nonlocal-flux
# variant (DEAD in production) is run_lifecycle_kppnonlcl_native_gate_core2.sh.
#
#   tools/run_lifecycle_kpp_native_gate_core2.sh [run_dir] [nsteps] [whichEVP]
set -euo pipefail
F3=/home/a/a270088/fesom3
RUN="${1:-/scratch/a/a270088/lifecycle_kpp_native_core2}"
NSTEPS="${2:-3}"
WHICHEVP="${3:-0}"

MIX_KPP=1 SW_PENE=1 FER_GM=1 REDI=1 \
    bash "$F3/tools/run_lifecycle_fullynative_gate_core2.sh" "$RUN" "$NSTEPS" "$WHICHEVP"
