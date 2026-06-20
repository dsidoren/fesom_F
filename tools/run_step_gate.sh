#!/usr/bin/env bash
# End-to-end M2.9b step-ASSEMBLY byte-gate. Drives the REAL FESOM2 oce_timestep_ale
# (compute_vel_nodes -> pressure_bv -> PP+mo_convect -> compute_vel_rhs -> viscosity ->
# impl_vert_visc_ale -> SSH CG -> update_vel/hbar/eta_n -> vert_vel_ale ->
# solve_tracers_ale -> update_thickness_ale) on prescribed pi 1-rank state, and the
# FESOM3 assembled step (mod_step_oce::step_oce) on the IDENTICAL state, then compares
# every per-substep NODE dump (density/pressure/bvfreq / Kv / ssh_rhs / d_eta / hbar /
# eta_n / hnode_new / w / T / S / hnode) for max|delta|=0 with tools/dump_diff.py.
# The SW_AB substep (id=2: sw_alpha/sw_beta, dead in M2 — KPP/GM/Redi only) is FESOM2-
# only and intentionally ignored.
#
# Unlike run_pressure_gate (kernels in isolation), this gates the ASSEMBLY: each kernel
# reads the PREVIOUS kernel's LIVE output (uvnode from compute_vel_nodes, Av/Kv from PP,
# d_eta from CG, ...), so a match proves the data flow + dispatch are byte-faithful.
#
# Prereqs (see docs/HANDOFF.md):
#   - FESOM2 rebuilt with src/fesom_step_dump.F90 wired into ocean_setup
#     (build/lib64/libfesom.so; cmake re-run to GLOB the new file).
#   - pi mesh has the hand-crafted dist_1/.
#   - FESOM3 CLEAN-built: build_intel_dp/bin/fesom_stepdump.
set -euo pipefail
F3=/home/a/a270088/fesom3
RUN="${1:-/scratch/a/a270088/stepdump_pi}"

echo "[1/3] FESOM2 oracle: REAL oce_timestep_ale (1-rank pi)"
bash "$F3/tools/run_stepdump_pi.sh" "$RUN" "$RUN/step_f2" >/dev/null

echo "[2/3] FESOM3: assembled step_oce (1-rank pi)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM_DUMP_FILE="$RUN/step_f3" FESOM_DUMP_MAXSTEPS=1
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_stepdump" >/dev/null 2>&1

echo "[3/3] compare (SW_AB substep 2 ignored — dead in M2)"
python3 "$F3/tools/dump_diff.py" "$RUN/step_f2" "$RUN/step_f3" --glob --ignore-substep=2
