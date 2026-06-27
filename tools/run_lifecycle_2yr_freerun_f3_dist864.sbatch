#!/bin/bash
# M8c Step 3 (physical validation path) — FREE-RUNNING 2-year FESOM3 stability run at dist_864.
# Per the day-107 finding: the 2-yr headline is NOT bit-identical to FESOM2 because of spurious
# DENORMAL m_snow (~1e-308 m, physically-zero snow) flipping the ice albedo branch — a latent
# FESOM2 quirk both codes generate, NOT a port bug (all snow/ice code + flags byte-identical,
# m_ice/a_ice exactly clean). So instead of max|delta|=0 we VALIDATE the 2-year FESOM3 run is
# PHYSICALLY stable & sane: full production config (JRA55-do 1958 + zstar+TKE+GM+Redi+sw_pene +
# Sweeney chl), 35040 steps (2 model years), NO oracle / NO dump. Diagnostics: per-step GLOBAL
# (MPI_MAX over all ranks) max|eta|/max|uv| stay bounded (no blow-up), the cryosphere grows &
# stays bounded (a_ice in [0,1], m_ice), warm/salty ceilings Tmax/Smax don't drift away, no
# NaN/Inf, and the rollover stack fires.
#
#   sbatch tools/run_lifecycle_2yr_freerun_f3_dist864.sbatch
#   NSTEPS=17520 sbatch tools/...     # 1 model year
#SBATCH --job-name=m8c_2yr_free
#SBATCH -p compute
#SBATCH --nodes=7
#SBATCH --ntasks=864
#SBATCH --time=02:00:00
#SBATCH -A ab0995
#SBATCH -o /scratch/a/a270088/m8c_2yr_free.%j.out
#SBATCH -e /scratch/a/a270088/m8c_2yr_free.%j.err
set -uo pipefail
F3=/home/a/a270088/fesom3
NSTEPS="${NSTEPS:-35040}"          # 2 model years (1958+1959, both non-leap) at 48 steps/day
WHICHEVP="${WHICHEVP:-0}"
START_CLOCK="${START_CLOCK:-0 1 1958}"
RUN="/scratch/a/a270088/m8c_2yr_free_evp${WHICHEVP}_n${NSTEPS}"
mkdir -p "$RUN"

COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0

source "$F3/env.sh" intel >/dev/null 2>&1          # KNEM fix (L35)
if [ -n "${SLURM_JOB_ID:-}" ]; then ulimit -s 204800; else ulimit -s unlimited; fi

# Full PRODUCTION config (matches the gate's F3 leg) — but NO FESOM_DUMP_FILE => dumps inactive,
# and no oracle leg. The fw_* diagnostic dumps in oce_fluxes are gated by the (inactive) dump, so
# they are no-ops here; physics is the production path.
export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE"
export FESOM3_FORCING_DIR="$POOL"  FESOM3_FORCING=JRA55
export FESOM3_START_CLOCK="$START_CLOCK"
export FESOM3_RUNOFF_FILE="$POOL/CORE2_runoff.nc"
export FESOM3_SSS_FILE="$POOL/PHC2_salx.nc"
export FESOM3_CHL_SWEENEY=1
export FESOM3_WHICH_ALE=zstar FESOM3_FER_GM=1 FESOM3_REDI=1 FESOM3_MIX_TKE=1 FESOM3_SW_PENE=1
export FESOM3_WHICHEVP="$WHICHEVP" FESOM3_NSTEPS="$NSTEPS"
export FESOM3_STEP_DIAG=1   # this run's PURPOSE is the per-step global MPI_MAX peak trajectory (gated OFF by default since M8-perf)

date
echo "################################################################"
echo "## M8c 2-YEAR FREE-RUN (F3 only): $NSTEPS steps, dist_864, whichEVP=$WHICHEVP"
echo "##   JRA55-do 1958 + zstar+TKE+GM+Redi+sw_pene + Sweeney chl ; START_CLOCK='$START_CLOCK'"
echo "################################################################"
rc=0
srun -l -n 864 "$F3/build_intel_dp/bin/fesom_lifecycle_native_mr" > "$RUN/run_f3.log" 2>&1 || rc=1

echo "## completion:"
grep -aE "done \($NSTEPS steps|done \([0-9]+ steps" "$RUN/run_f3.log" | tail -1 || echo "   !! did NOT reach 'done' — check $RUN/run_f3.log"
echo "## NaN/Inf check (should be EMPTY):"
grep -aiE "NaN|Infinity| Inf |\*\*\*\*" "$RUN/run_f3.log" | head -5 || true
grep -aciE "NaN|Infinity" "$RUN/run_f3.log" | awk '{print "   NaN/Inf line count = "$1}'
echo "## rollover events:"
grep -aE "YEAR ROLLOVER|roll_monthly_clim:" "$RUN/run_f3.log" | sort | uniq -c | head -30 || echo "   (none logged)"
echo "## peak |eta| / |uv| / |a_ice| trajectory (every ~month = 1488 steps):"
grep -aE "fesom_lifecycle_native_mr: step" "$RUN/run_f3.log" | awk 'NR==1 || NR%1488==0 || NR==35040' | tail -30 || true
date
echo "Exit rc=$rc"
exit $rc
