#!/usr/bin/env bash
# M3d oracle: run the REAL FESOM2 init (use_ice=.true.) on CORE2 1-rank up to ice_setup,
# then the ice_thermo_dump_write shim (FESOM_THERMO_DUMP) prescribes the analytic ocean +
# ATMOSPHERIC forcing, runs the REAL ocean2ice + EVPdynamics (120 subcycles) -> uice/vice,
# the REAL ice FCT advection (ice_TG_rhs + ice_fct_solve), then the REAL cut_off +
# thermodynamics on the advected ice + prescribed atmosphere + do_ic3d surface T/S, dumps
# a_ice/m_ice/m_snow + flx_h/flx_fw + t_skin and STOPS (before the time loop). Same forced
# init path as run_icefctdump_core2.sh (CORE2 NCAR forcing stubs + runoff/SSS, np=1
# next_io_rank fix, reduced-M2 dynamics) — env-gated to the thermo dump.
#
# The thermo SCALAR config comes from the CORE2 namelists (the FESOM3 driver matches them):
#   - namelist.ice  &ice_therm: con/consn/Sice/h0/h0_s/hmin/Armin/emiss/albedos/albw=0.1 ...
#   - namelist.forcing  &forcing_exchange_coeff: Ch_atm_ice=Ce_atm_ice=0.00175; &nam_sbc l_snow=.true.
#   - namelist.tra: ref_sss_local=.true. (=> rsss=S_oc); which_ALE='linfs' => use_virt_salt=.true.
# Calendar: CORE (noleap) forcing REQUIRES include_fleapyear=.false. use_sw_pene=.false.
#
#   tools/run_icethermodump_core2.sh [run_dir] [thermo_dump_prefix] [whichEVP]
# whichEVP: 0 = standard EVP (default), 1 = modified EVP (mEVP).
set -euo pipefail
F2=/home/a/a270088/port2/fesom2
STUB="$F2/test/input/global"
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
RUN="${1:-/scratch/a/a270088/icethermodump_core2}"
THERMODUMP="${2:-$RUN/icethermo_f2.bin}"
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
export FESOM_THERMO_DUMP="$THERMODUMP"
ulimit -s unlimited
echo "run_icethermodump_core2: 1 rank (FORCED use_ice, whichEVP=${WHICHEVP}, stops after thermodynamics) -> ${THERMODUMP}"
timeout 600 mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 ./fesom.x > "$RUN/run.log" 2>&1 || true
if [ -f "$THERMODUMP" ]; then
    echo "run_icethermodump_core2: done -> $(stat -c%s "$THERMODUMP")B"
    grep -E 'ice_thermo_dump_write|Ice is initialized|EVP scheme' "$RUN/run.log" | head
else
    echo "run_icethermodump_core2: MISSING DUMP — see $RUN/run.log"; tail -30 "$RUN/run.log"; exit 1
fi
