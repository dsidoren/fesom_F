#!/usr/bin/env bash
# M5c ghats-variant byte-gate (CORE2 1-rank): the production KPP+sw_pene+GM+Redi native lifecycle
# PLUS the KPP nonlocal counter-gradient flux turned ON (KPP_NONLCL=1 -> use_kpp_nonlclflx=.true.
# on BOTH sides). This term (FESOM2 oce_ale_tracer.F90:892-987, the ghats*blmc surface/bulk/bottom
# flux on T and S) is DEAD in the production work_core config (use_kpp_nonlclflx defaults .false.
# and is absent from the namelists), so the production gate (run_lifecycle_kpp_native_gate_core2.sh)
# never exercises it. This variant turns it on so the otherwise-untested ghats term is byte-verified
# rather than left as ungated code (the L39 ungated-landmine lesson).
#
# Compares the 13 NODE substeps x 5 probes over N steps (--ignore-substep=2).
#
#   tools/run_lifecycle_kppnonlcl_native_gate_core2.sh [run_dir] [nsteps] [whichEVP]
set -euo pipefail
F3=/home/a/a270088/fesom3
RUN="${1:-/scratch/a/a270088/lifecycle_kppnonlcl_native_core2}"
NSTEPS="${2:-3}"
WHICHEVP="${3:-0}"

MIX_KPP=1 SW_PENE=1 FER_GM=1 REDI=1 KPP_NONLCL=1 \
    bash "$F3/tools/run_lifecycle_fullynative_gate_core2.sh" "$RUN" "$NSTEPS" "$WHICHEVP"
