#!/usr/bin/env bash
# M6a-4 MULTI-RANK forced/native zstar byte-gate (CORE2 dist_<NP>). The REAL FESOM2 forced
# lifecycle with which_ALE='zstar' at NP ranks vs FESOM3's FULLY-NATIVE MULTI-RANK lifecycle
# (fesom_lifecycle_native_mr, FESOM3_WHICH_ALE=zstar). The zstar kernels are optional-`partit`
# from the start: update_stiff_mat_ale / vert_vel_ale / update_thickness_ale loop owned, and the
# exchanges are in place (exchange_elem(helem), exchange_nod(Wvel/hnode_new/ssh_rhs/ssh_rhs_old));
# dhe is read at owned el(i) only (both triangles of an owned edge are owned — M2.12c-2 invariant).
# So this is pure WIRING (the M4f/M5d lesson). Reduced column physics (PP) by default to ISOLATE
# the zstar + freshwater MR wiring; add FER_GM=1 REDI=1 MIX_KPP=1 SW_PENE=1 for the full production.
#
# ⚠️ multi-rank levante needs the env.sh KNEM flag (L35) — sourced by the gate runner.
#
#   tools/run_lifecycle_zstar_native_gate_multirank.sh [np] [nsteps] [whichEVP] [run_dir]
#   (env passthrough: FER_GM REDI MIX_KPP SW_PENE)
set -euo pipefail
F3=/home/a/a270088/fesom3
NP="${1:-2}"
NSTEPS="${2:-3}"
WHICHEVP="${3:-0}"
RUN="${4:-/scratch/a/a270088/lifecycle_zstar_native_mr${NP}}"

WHICH_ALE=zstar bash "$F3/tools/run_lifecycle_fullynative_gate_multirank.sh" "$NP" "$NSTEPS" "$WHICHEVP" "$RUN"
