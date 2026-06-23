#!/usr/bin/env bash
# M3e oracle: run the REAL FESOM2 init (use_ice=.true.) on CORE2 1-rank up to ice_setup,
# then the ice_flux_dump_write shim (FESOM_OCEFLUX_DUMP) prescribes the analytic ocean +
# atmospheric forcing + atm-ocean momentum stress (stress_atmoce_x/y) + SSS-restoring
# climatology (Ssurf), runs the REAL ocean2ice + EVPdynamics (120 subcycles) -> uice/vice,
# the REAL ice FCT advection (ice_TG_rhs + ice_fct_solve), cut_off + thermodynamics, then the
# REAL oce_fluxes_mom + oce_fluxes (the air-sea coupling-out), dumps heat_flux/water_flux/
# virtual_salt/relax_salt + stress_surf and STOPS (before the time loop). Same forced init
# path as run_icethermodump_core2.sh (CORE2 NCAR forcing stubs + runoff/SSS, np=1 next_io_rank
# fix, reduced-M2 dynamics) — env-gated to the M3e flux dump.
#
# The flux-balancing config comes from the CORE2 namelists (the FESOM3 driver matches them):
#   - namelist.tra: surf_relax_S=1.929e-06, ref_sss_local=.true. (=> rsss=S_oc), ref_sss=34. (dead),
#     which_ALE forced to 'linfs' => use_virt_salt=.true.
#   - namelist.forcing: Ch_atm_ice=Ce_atm_ice=0.00175; l_snow=.true.
#   - density_0=1030 (o_PARAM), ocean_area (geometry). use_sw_pene=.false. (cal_shortwave_rad skipped).
# Calendar: CORE (noleap) forcing REQUIRES include_fleapyear=.false.
#
#   tools/run_icefluxdump_core2.sh [run_dir] [flux_dump_prefix] [whichEVP]
# whichEVP: 0 = standard EVP (default), 1 = modified EVP (mEVP).
set -euo pipefail
F2=/home/a/a270088/port2/fesom2
STUB="$F2/test/input/global"
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
RUN="${1:-/scratch/a/a270088/icefluxdump_core2}"
FLUXDUMP="${2:-$RUN/iceflux_f2.bin}"
WHICHEVP="${3:-0}"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1
rm -rf "$RUN"; mkdir -p "$RUN"
cp "$F2"/work_core/namelist.* "$RUN"/
cp "$F2"/work_core/namelist.forcing.CORE2 "$RUN"/namelist.forcing
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x
printf '0 1 1948\n0 1 1948\n' > "$RUN"/fesom.clock

python3 - "$RUN" "$STUB" "$POOL" "$WHICHEVP" <<'PY'
import re,sys
run,stub,pool,whichevp=sys.argv[1:5]
p=run+'/namelist.ice'; s=open(p).read()
s=re.sub(r"whichEVP\s*=\s*\d+",f"whichEVP       = {whichevp}",s,1)
open(p,'w').write(s)
p=run+'/namelist.config'; s=open(p).read()
s=re.sub(r"ResultPath\s*=\s*'[^']*'","ResultPath       = './'",s,1)
s=re.sub(r"which_ALE\s*=\s*'zlevel'","which_ALE          = 'linfs'",s,1)
s=re.sub(r"yearnew\s*=\s*1958","yearnew = 1948",s,1)
s=re.sub(r"run_length\s*=\s*\d+","run_length        = 1",s,1)
s=re.sub(r"include_fleapyear\s*=\s*\.true\.","include_fleapyear = .false.",s,1)   # CORE noleap
s=re.sub(r"use_sw_pene\s*=\s*\.true\.","use_sw_pene              = .false.",s,1)  # match ported step_oce
open(p,'w').write(s)
p=run+'/namelist.oce'; s=open(p).read()
s=re.sub(r"mix_scheme\s*=\s*'KPP'","mix_scheme         = 'PP'",s,1)
s=re.sub(r"Fer_GM\s*=\s*\.true\.","Fer_GM             = .false.",s,1)
s=re.sub(r"Redi\s*=\s*\.true\.","Redi               = .false.",s,1)
open(p,'w').write(s)
p=run+'/namelist.forcing'; s=open(p).read()
s=s.replace("FORCING/CORE2/", stub+"/")                                          # atm stubs absolute
s=re.sub(r"nm_runoff_file\s*=\s*'[^']*'",   f"nm_runoff_file ='{pool}/CORE2_runoff.nc'",s,1)
s=re.sub(r"nm_sss_data_file\s*=\s*'[^']*'", f"nm_sss_data_file ='{pool}/PHC2_salx.nc'",s,1)
open(p,'w').write(s)
PY

cd "$RUN"
export FESOM_OCEFLUX_DUMP="$FLUXDUMP"
ulimit -s unlimited
echo "run_icefluxdump_core2: 1 rank (FORCED use_ice, whichEVP=${WHICHEVP}, stops after oce_fluxes) -> ${FLUXDUMP}"
timeout 600 mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 ./fesom.x > "$RUN/run.log" 2>&1 || true
if [ -f "$FLUXDUMP" ]; then
    echo "run_icefluxdump_core2: done -> $(stat -c%s "$FLUXDUMP")B"
    grep -E 'ice_flux_dump_write|Ice is initialized|EVP scheme' "$RUN/run.log" | head
else
    echo "run_icefluxdump_core2: MISSING DUMP — see $RUN/run.log"; tail -30 "$RUN/run.log"; exit 1
fi
