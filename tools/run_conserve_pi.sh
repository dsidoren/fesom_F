#!/usr/bin/env bash
# CONSERVATION + NO-LEAKAGE GATE for the bottom-at-vertices change (docs/plans/
# 20260910-fesom3-bottom-at-vertices.md, Task 2). This is the numerical regression net
# that replaces the retired FESOM2 byte-gates: it is FESOM_F-internal, needs no oracle,
# and must stay green through every task of that plan.
#
# What it checks, via src/drivers/fesom_conserve.F90:
#   zstar + Redi, np 1 -> the Redi path. diff_ver_part_redi_expl is the ONLY caller of
#       the tr_xynodes node-average whose denominator the depth-independent area changed,
#       and it never runs unless Redi is on. Note that conservation alone does NOT validate
#       that denominator -- the Redi tendency is a telescoping flux divergence, so it
#       conserves whatever tr_xynodes contains. The denominator is pinned by test_bottom's
#       'Redi: node-averaged gradient is exact on the WET area' check instead.
#   zstar, np 1 and 2 -> HARD conservation gate. Total heat and salt content must not
#       drift beyond TOL. This is the sharp test of the bottom change: if the flux areas
#       and the cell volumes ever disagree -- the central risk when area() becomes
#       depth-independent and nlevels() flips to vertex-defined -- it shows up here as
#       drift long before it shows up anywhere else.
#   linfs, np 1 -> INVARIANTS ONLY (T10 velocity leakage, helem == mean(hnode),
#       zbar_e_bot, finiteness). Drift is reported but NOT gated, because the linear free
#       surface is not tracer-conserving by construction: hnode is frozen, so the surface
#       vertical advective flux -w*T*area at nzmin is a real source/sink with no thickness
#       change to balance it. Measured baseline on pi over 20 steps: heat -2.2e-04,
#       salt -1.2e-06, non-monotone (a free-surface adjustment transient).
#
# Measured zstar baseline on pi, 20 steps, before any bottom change:
#   np=1  heat -4.87e-15  salt -1.66e-14
#   np=2  heat  0.00e+00  salt -7.80e-15
# TOL=1e-12 leaves roughly two orders of headroom over that.
#
# Usage:  bash tools/run_conserve_pi.sh [NSTEPS]
set -euo pipefail

F3="${F3:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BIN="${BIN:-$F3/build_intel_dp/bin/fesom_conserve}"
PIMESH="${PIMESH:-/home/a/a270088/port2/fesom2/tests/data/MESHES/pi}"
ICFILE="${ICFILE:-/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc}"
NSTEPS="${1:-20}"
TOL="${TOL:-1e-12}"
MPIFLAGS="${MPIFLAGS:---mca pml ob1 --mca btl self,vader --oversubscribe}"

[ -x "$BIN" ] || { echo "run_conserve_pi: missing $BIN (build first)"; exit 1; }
[ -d "$PIMESH" ] || { echo "run_conserve_pi: missing pi mesh $PIMESH"; exit 1; }

export FESOM3_MESH_DIR="$PIMESH" FESOM3_IC_FILE="$ICFILE" FESOM3_NSTEPS="$NSTEPS"

fail=0

for np in 1 2; do
    echo "=== zstar, np=$np, $NSTEPS steps (conservation gate, tol=$TOL) ==="
    if FESOM3_WHICH_ALE=zstar FESOM3_CONSERVE_TOL="$TOL" \
         mpirun $MPIFLAGS -n "$np" "$BIN" 2>&1 | tail -6; then
        :
    else
        echo "run_conserve_pi: FAILED (zstar np=$np)"; fail=1
    fi
done

# Physics combinations. The area/areasvol consumers are spread across the vertical
# diffusion TDMA, the Redi isoneutral flux and the KPP non-local / shortwave terms, and
# each is only reachable with its own scheme switched on -- a gate that runs only the
# default configuration proves nothing about the others.
while read -r label vars; do
    [ -z "$label" ] && continue
    echo "=== zstar + $label, np=1, $NSTEPS steps (conservation gate, tol=$TOL) ==="
    if env $vars FESOM3_WHICH_ALE=zstar FESOM3_CONSERVE_TOL="$TOL" \
         mpirun $MPIFLAGS -n 1 "$BIN" 2>&1 | tail -7; then
        :
    else
        echo "run_conserve_pi: FAILED (zstar+$label np=1)"; fail=1
    fi
done <<'CFG'
Redi FESOM3_REDI=1
KPP FESOM3_MIX_KPP=1
TKE FESOM3_MIX_TKE=1
GM+Redi+TKE FESOM3_FER_GM=1 FESOM3_REDI=1 FESOM3_MIX_TKE=1
GM+KPP FESOM3_FER_GM=1 FESOM3_MIX_KPP=1
CFG

echo "=== linfs, np=1, $NSTEPS steps (invariants only; drift reported, not gated) ==="
if FESOM3_WHICH_ALE=linfs mpirun $MPIFLAGS -n 1 "$BIN" 2>&1 | tail -5; then
    :
else
    echo "run_conserve_pi: FAILED (linfs invariants np=1)"; fail=1
fi

if [ "$fail" -ne 0 ]; then
    echo "run_conserve_pi: GATE FAILED"
    exit 1
fi
echo "run_conserve_pi: GATE OK"
