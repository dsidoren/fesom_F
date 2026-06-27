#!/usr/bin/env bash
# M8c JRA55-do PRODUCTION forcing byte-gate (CORE2 mesh dist_<NP>). The production target the user
# actually runs: JRA55-do-v1.4.0 atmosphere (3-hourly uas/vas/huss/tas/rsds/rlds/prra/prsn,
# gregorian calendar + include_fleapyear=.true.), zstar + cvmix_TKE + GM + Redi + sw_pene, native
# ice/forcing, start year 1958. This is the FIRST fesom3 byte-gate on JRA55 (M2.10->M8b used the
# CORE2 NCAR substitute); it exercises the never-before-gated gregorian julday/calendar_date paths
# + the tmid=0 mid-point time-axis shift, and (because JRA55 is 3-hourly) the forcing_sbc_do record
# crossings fire every 6 steps even within day 1.
#
# CHL_SWEENEY=0 (default, step-1 isolate) => chl='None'/const 0.1 on BOTH codes (matches the F3
# constant chl) so the gate isolates the JRA55 atmosphere read from the chl-Sweeney monthly read.
# CHL_SWEENEY=1 = the full production config (needs the F3 chl-Sweeney read — M8c step 2).
#
#   tools/run_lifecycle_jra55_gate_multirank.sh [np] [nsteps] [whichEVP] [run_dir]
#   START_CLOCK="0 d 1958" overrides the cold start; CHL_SWEENEY=1 for the full production gate.
set -euo pipefail
F3=/home/a/a270088/fesom3
JRA55POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
NP="${1:-2}"
NSTEPS="${2:-12}"
WHICHEVP="${3:-0}"
RUN="${4:-/scratch/a/a270088/lifecycle_jra55_mr${NP}}"
START_CLOCK="${START_CLOCK:-0 1 1958}"
CHL_SWEENEY="${CHL_SWEENEY:-0}"

FORCING_SET=JRA55 CHL_SWEENEY="$CHL_SWEENEY" \
FORCING_OVERRIDE="$JRA55POOL" START_CLOCK="$START_CLOCK" \
WHICH_ALE=zstar MIX_TKE=1 SW_PENE=1 FER_GM=1 REDI=1 \
    bash "$F3/tools/run_lifecycle_fullynative_gate_multirank.sh" "$NP" "$NSTEPS" "$WHICHEVP" "$RUN"
