#!/usr/bin/env bash
# End-to-end M2.1 EOS/pressure/N^2 + M2.2 PGF + M2.3/M2.4 vel_rhs + M2.4 biharmonic
# viscosity + M2.5 implicit vertical viscosity operator byte-gate: dump FESOM2's
# density_m_rho0/hpressure/bvfreq/pgf/coriolis/uvnode_rhs/uv_rhs*/visc_u_c/visc_v_c/
# uv_rhs_visc/uv_rhs_ivv (1-rank pi, the REAL pressure_bv +
# pressure_force_4_linfs_fullcell + compute_vel_rhs[momadv_opt=2] +
# visc_filt_bidiff[opt_visc=7] + impl_vert_visc_ale on prescribed T/S/UV/eta/w_e/w_i/Av/
# stress_surf), dump FESOM3's (1-rank pi, transcribed oce_pressure_bv -> oce_pgf ->
# oce_dyn_velrhs incl. momentum_adv_scalar -> oce_dyn_visc -> oce_dyn_ivertvisc on the
# identical prescribed inputs), and compare for max|delta|=0 on all 28 fields. The vel_rhs
# is the FULL UV_rhs (Coriolis + AB2 + PGF + SSH gradient + momentum advection); biharmonic
# viscosity then the implicit-vertical-viscosity TDMA are SEPARATE operators run after it.
# The M2.4 momadv intermediate uvnode_rhs and the viscosity 1st-Laplacian visc_u_c/visc_v_c
# are gated separately. (The driver prints a visc strength diagnostic: with UV bumped to
# 2.0/1.5 m/s, the flow-aware gamma0/gamma1/gamma2 selection fires on ~20/79/1% of edges;
# and an ivertvisc diagnostic: max|d(uv_rhs)| and the w_i>0/<0 split.)
#
# Prereqs (one-time, see docs/HANDOFF.md):
#   - FESOM2 rebuilt with src/fesom_pressure_dump.F90 wired into ocean_setup
#     (build/lib64/libfesom.so; fesom.x loads it at runtime).
#   - pi mesh has the hand-crafted dist_1/ (shared with the geometry/advection gates).
#   - FESOM3 built: build_intel_dp/bin/fesom_pressuredump.
set -euo pipefail
F3=/home/a/a270088/fesom3
RUN="${1:-/scratch/a/a270088/pressuredump_pi}"

echo "[1/3] FESOM2 oracle pressure dump (1-rank pi)"
bash "$F3/tools/run_pressuredump_pi.sh" "$RUN" "$RUN/pressure_f2.bin" >/dev/null

echo "[2/3] FESOM3 pressure dump (1-rank pi)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_PRESSURE_OUT="$RUN/pressure_f3.bin"
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_pressuredump" >/dev/null 2>&1

echo "[3/3] compare"
python3 "$F3/tools/pressure_diff.py" "$RUN/pressure_f2.bin" "$RUN/pressure_f3.bin"
