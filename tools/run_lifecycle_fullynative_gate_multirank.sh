#!/usr/bin/env bash
# M3f-4 MULTI-RANK FULLY-NATIVE forced lifecycle byte-gate (CORE2 dist_<NP>). The REAL forced
# FESOM2 lifecycle (use_ice=.true., real CORE2 forcing + ice EVP + oce_fluxes) at NP ranks vs
# FESOM3's FULLY NATIVE multi-rank lifecycle: FESOM3 reads NOTHING from the oracle atmosphere —
# it computes the WHOLE air-sea forcing itself each step over owned+halo (8 NCAR + NCAR bulk +
# stresses + runoff + Ssurf), then ocean2ice -> ice_timestep -> oce_fluxes -> step_oce, all
# through the optional partit. Both codes share dist_<NP>/myList, so each global probe id is
# owned by the same rank on both (the L8 same-partition rule); the per-rank gid-keyed dump_shim
# / mod_dump are matched by GLOBAL id (dump_diff.py --glob). 13 NODE substeps x 5 probes x N
# steps (SW_AB substep 2 ignored). NOTE: the oracle flux/atmflux dumps skip at npes/=1, so the
# node-substep gate is the validation (the 1-rank fully-native self-check proved the fluxes).
#
#   tools/run_lifecycle_fullynative_gate_multirank.sh [np] [nsteps] [whichEVP] [run_dir]
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
NP="${1:-2}"
NSTEPS="${2:-3}"
WHICHEVP="${3:-0}"
RUN="${4:-/scratch/a/a270088/lifecycle_fullynative_mr${NP}}"
# M4f: FER_GM=1 / REDI=1 turn ON the work_core GM bolus / Redi isopycnal diffusion in BOTH the
# oracle and the FESOM3 MR native lifecycle. Default 0 = the proven M3f-4 GM/Redi-off no-regr gate.
FER_GM="${FER_GM:-0}"
REDI="${REDI:-0}"
# M5d: MIX_KPP=1 swaps PP->KPP, SW_PENE=1 turns on shortwave penetration, KPP_NONLCL=1 turns on the
# ghats nonlocal flux — in BOTH the oracle and the FESOM3 MR native lifecycle. Defaults 0.
MIX_KPP="${MIX_KPP:-0}"
SW_PENE="${SW_PENE:-0}"
KPP_NONLCL="${KPP_NONLCL:-0}"
# M7d: MIX_TKE=1 runs mix_scheme='cvmix_TKE' (prognostic TKE) at MULTI-RANK in BOTH sides. Pair
# with SW_PENE=1 (work_*_tke use_sw_pene=.true.). Default 0.
MIX_TKE="${MIX_TKE:-0}"
# M6a-4: WHICH_ALE=zstar runs the production free surface + real freshwater flux at MULTI-RANK.
WHICH_ALE="${WHICH_ALE:-linfs}"
# M8b: FORCING_OVERRIDE is the atmosphere directory read by BOTH codes (oracle FORCING/CORE2 swap +
# FESOM3_FORCING_DIR). Default = the day-1 STUB so existing gates are byte-unchanged; the long gate
# points it at the full-year CORE2 pool. START_CLOCK seeds the cold start (oracle fesom.clock +
# FESOM3_START_CLOCK) so a non-day-1 start is byte-gated identically on both sides.
STUBDIR=/home/a/a270088/port2/fesom2/test/input/global
FORCING_OVERRIDE="${FORCING_OVERRIDE:-$STUBDIR}"
START_CLOCK="${START_CLOCK:-0 1 1948}"
# M8c: FORCING_SET=JRA55 runs the production JRA55-do atmosphere (gregorian + leap) on BOTH codes
# (oracle uses work_core's JRA55 namelist.forcing; FESOM3 sets FESOM3_FORCING=JRA55). CHL_SWEENEY=1
# = production monthly chl on BOTH; default 0 = chl='None'/const 0.1 (isolates the atmosphere read).
FORCING_SET="${FORCING_SET:-CORE2}"
CHL_SWEENEY="${CHL_SWEENEY:-0}"

echo "[1/3] FESOM2 oracle FORCED lifecycle ($NSTEPS steps, $NP-rank CORE2 dist_$NP, use_ice, whichEVP=$WHICHEVP; WHICH_ALE=$WHICH_ALE FER_GM=$FER_GM REDI=$REDI MIX_KPP=$MIX_KPP MIX_TKE=$MIX_TKE SW_PENE=$SW_PENE NONLCL=$KPP_NONLCL)"
WHICH_ALE="$WHICH_ALE" FER_GM="$FER_GM" REDI="$REDI" MIX_KPP="$MIX_KPP" MIX_TKE="$MIX_TKE" SW_PENE="$SW_PENE" KPP_NONLCL="$KPP_NONLCL" \
    FORCING_OVERRIDE="$FORCING_OVERRIDE" START_CLOCK="$START_CLOCK" FORCING_SET="$FORCING_SET" CHL_SWEENEY="$CHL_SWEENEY" \
    bash "$F3/tools/run_lifecycle_forced_core2.sh" "$RUN" "$RUN/lifef_f2" "$RUN/flux_f2" "$NSTEPS" "$RUN/atmflux_f2" "$WHICHEVP" "$NP" | tail -4

echo "[2/3] FESOM3 FULLY NATIVE MR lifecycle ($NSTEPS steps, $NP-rank dist_$NP, whichEVP=$WHICHEVP) — NO prescribed atmosphere"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE"
export FESOM3_FORCING_DIR="$FORCING_OVERRIDE"  # M8b: native NCAR read (stub day-1 or CORE2 pool)
export FESOM3_START_CLOCK="$START_CLOCK"        # M8b: cold start (matches the oracle fesom.clock)
export FESOM3_RUNOFF_FILE="$POOL/CORE2_runoff.nc"   # native runoff (M3f-3b)
export FESOM3_SSS_FILE="$POOL/PHC2_salx.nc"          # native SSS restoring (M3f-3b)
export FESOM3_WHICHEVP="$WHICHEVP"
if [ "$FORCING_SET" = JRA55 ]; then export FESOM3_FORCING=JRA55; fi  # M8c: JRA55-do production atmosphere
if [ "$CHL_SWEENEY" = 1 ]; then export FESOM3_CHL_SWEENEY=1; fi       # M8c: production Sweeney monthly chl read
if [ "$WHICH_ALE" != linfs ]; then export FESOM3_WHICH_ALE="$WHICH_ALE"; fi  # M6a-4: zstar MR
if [ "$FER_GM" = 1 ]; then export FESOM3_FER_GM=1; fi   # M4f: GM bolus in the MR native lifecycle
if [ "$REDI" = 1 ]; then export FESOM3_REDI=1; fi       # M4f: Redi isopycnal diffusion
if [ "$MIX_KPP" = 1 ]; then export FESOM3_MIX_KPP=1; fi # M5d: KPP vertical mixing
if [ "$MIX_TKE" = 1 ]; then export FESOM3_MIX_TKE=1; fi # M7d: TKE vertical mixing (cvmix_TKE)
if [ "$SW_PENE" = 1 ]; then export FESOM3_SW_PENE=1; fi # M5d: shortwave penetration (sw_3d term)
if [ "$KPP_NONLCL" = 1 ]; then export FESOM3_KPP_NONLCL=1; fi # M5d: ghats nonlocal flux (gate variant)
export FESOM_DUMP_FILE="$RUN/lifen_f3" FESOM_DUMP_MAXSTEPS="$NSTEPS" FESOM3_NSTEPS="$NSTEPS"
# M8b: login needs unlimited stack (L20); Levante compute nodes reject it -> batch uses 204800 KB.
if [ -n "${SLURM_JOB_ID:-}" ]; then ulimit -s 204800; else ulimit -s unlimited; fi
# M8b: srun inside a SLURM allocation (dist_864), else the proven login-node mpirun (byte-unchanged).
if [ -n "${SLURM_JOB_ID:-}" ]; then
    srun -l -n "$NP" "$F3/build_intel_dp/bin/fesom_lifecycle_native_mr" > "$RUN/run_f3.log" 2>&1 || \
        { echo "  FESOM3 run failed"; tail -40 "$RUN/run_f3.log"; exit 1; }
else
    mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n "$NP" \
        "$F3/build_intel_dp/bin/fesom_lifecycle_native_mr" > "$RUN/run_f3.log" 2>&1 || \
        { echo "  FESOM3 run failed"; tail -40 "$RUN/run_f3.log"; exit 1; }
fi
grep -E 'nod2D|IC\(|native runoff|ENABLED|step [0-9]|done' "$RUN/run_f3.log" || true

echo "[3/3] compare per-rank (gid-keyed; 13 NODE substeps x 5 probes x $NSTEPS steps; SW_AB substep 2 ignored)"
python3 "$F3/tools/dump_diff.py" "$RUN/lifef_f2" "$RUN/lifen_f3" --glob --ignore-substep=2
