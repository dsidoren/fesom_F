#!/usr/bin/env bash
# M3b oracle: run the REAL FESOM2 init (use_ice=.true.) on CORE2 1-rank up to ice_setup,
# then the evp_dump_write shim (FESOM_EVP_DUMP) prescribes the analytic ocean-forcing
# state (surface UV / hbar / stress_atmice), runs the REAL ocean2ice + EVPdynamics (120
# EVP subcycles) on the cold-start ice IC, dumps srfoce_u/v + uice/vice + sigma11/12/22
# and STOPS (before the time loop). This is the run_icedump_core2.sh path (CORE2 NCAR
# forcing stubs + runoff/SSS, np=1 next_io_rank fix, reduced-M2 dynamics) — env-gated to
# the EVP dump. cd_oce_ice/delta_min come from namelist.ice (doubles).
#
# Calendar: CORE (noleap) forcing REQUIRES include_fleapyear=.false. use_sw_pene=.false.
# The shim fires in init (after ice_setup), so run_length is irrelevant.
#
#   tools/run_evpdump_core2.sh [run_dir] [evp_dump_prefix] [whichEVP]
# whichEVP: 0 = standard EVP (default), 1 = modified EVP (mEVP).
set -euo pipefail
F2=/home/a/a270088/port2/fesom2
STUB="$F2/test/input/global"
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
RUN="${1:-/scratch/a/a270088/evpdump_core2}"
EVPDUMP="${2:-$RUN/evp_f2.bin}"
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
export FESOM_EVP_DUMP="$EVPDUMP"
ulimit -s unlimited
echo "run_evpdump_core2: 1 rank (FORCED use_ice, whichEVP=${WHICHEVP}, stops after EVPdynamics) -> ${EVPDUMP}"
timeout 600 mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 ./fesom.x > "$RUN/run.log" 2>&1 || true
if [ -f "$EVPDUMP" ]; then
    echo "run_evpdump_core2: done -> $(stat -c%s "$EVPDUMP")B"
    grep -E 'evp_dump_write|Ice is initialized|EVP scheme' "$RUN/run.log" | head
else
    echo "run_evpdump_core2: MISSING DUMP — see $RUN/run.log"; tail -30 "$RUN/run.log"; exit 1
fi
