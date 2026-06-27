#!/usr/bin/env bash
# M8a-1 unit-level GATE for mod_clock: advance the clock from a cold-start 0 1 1948 and verify
# timenew/daynew/yearnew/month/day_in_month against a hand-computed table (day/month/year
# rollovers) + clock_nsteps for s/d/m/y units. Pure calendar arithmetic, no mesh/forcing.
#
#   tools/run_clocktest.sh [np]
set -euo pipefail
F3=/home/a/a270088/fesom3
NP="${1:-1}"
RUN="${RUN:-/scratch/a/a270088/clocktest}"
mkdir -p "$RUN"
source "$F3/env.sh" intel >/dev/null 2>&1
FESOM3_RESTART_IN="$RUN/" mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n "$NP" \
    "$F3/build_intel_dp/bin/fesom_clocktest"
