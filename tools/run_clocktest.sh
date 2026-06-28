#!/usr/bin/env bash
# Unit-level GATE for mod_clock. Two parts, both in the fesom_clocktest driver:
#   1. (M8a-1) advance the clock from a cold-start 0 1 1948 and verify
#      timenew/daynew/yearnew/month/day_in_month against a hand-computed table (day/month/year
#      rollovers) + clock_nsteps for s/d/m/y units.
#   2. (Task 0.1) clock_finish -> clock_init round-trip: write the .clock and read it straight
#      back, asserting BOTH lines round-trip exactly, that equal lines => r_restart=.false. and
#      differing lines => r_restart=.true., and the year-rollover branch
#      (daynew==ndpyr .and. timenew==86400 -> 0.0 / 1 / yearold+1).
# Pure calendar arithmetic + .clock file I/O, no mesh/forcing.
#
#   tools/run_clocktest.sh [np]
#   F3=<repo-or-worktree> tools/run_clocktest.sh   # gate a worktree build (defaults to main repo)
set -euo pipefail
F3="${F3:-/home/a/a270088/fesom3}"
BUILD="${BUILD:-build_intel_dp}"
NP="${1:-1}"
RUN="${RUN:-/scratch/a/a270088/clocktest}"
mkdir -p "$RUN"
source "$F3/env.sh" intel >/dev/null 2>&1
FESOM3_RESTART_IN="$RUN/" mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n "$NP" \
    "$F3/$BUILD/bin/fesom_clocktest"
