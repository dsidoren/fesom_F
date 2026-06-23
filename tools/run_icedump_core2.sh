#!/usr/bin/env bash
# M3a oracle: run the REAL FESOM2 init (use_ice=.true.) on CORE2 1-rank up to ice_setup,
# then the fesom_ice_dump shim (FESOM_ICE_DUMP) dumps the cold-start ice IC
# (a_ice/m_ice/m_snow) + the FCT mass matrix and STOPS (before the time loop). This is
# the run_lifecycle_forced_core2.sh init path (CORE2 NCAR forcing stubs + runoff/SSS, the
# np=1 next_io_rank fix, reduced-M2 dynamics) — but env-gated to the ice-foundation dump.
#
# Calendar: CORE (noleap) forcing REQUIRES include_fleapyear=.false. use_sw_pene=.false.
# (matches the ported step_oce). The shim fires in init, so run_length is irrelevant.
#
#   tools/run_icedump_core2.sh [run_dir] [ice_dump_prefix]
set -euo pipefail
F2=/home/a/a270088/port2/fesom2
STUB="$F2/test/input/global"
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
RUN="${1:-/scratch/a/a270088/icedump_core2}"
ICEDUMP="${2:-$RUN/ice_f2.bin}"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1
rm -rf "$RUN"; mkdir -p "$RUN"
cp "$F2"/work_core/namelist.* "$RUN"/
cp "$F2"/work_core/namelist.forcing.CORE2 "$RUN"/namelist.forcing
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x
printf '0 1 1948\n0 1 1948\n' > "$RUN"/fesom.clock

python3 - "$RUN" "$STUB" "$POOL" <<'PY'
import re,sys
run,stub,pool=sys.argv[1:4]
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
export FESOM_ICE_DUMP="$ICEDUMP"
ulimit -s unlimited
echo "run_icedump_core2: 1 rank (FORCED use_ice, stops after ice_setup) -> ${ICEDUMP}"
timeout 600 mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 ./fesom.x > "$RUN/run.log" 2>&1 || true
if [ -f "$ICEDUMP" ]; then
    echo "run_icedump_core2: done -> $(stat -c%s "$ICEDUMP")B"
    grep -E 'fesom_ice_dump|Ice is initialized|EVP scheme' "$RUN/run.log" | head
else
    echo "run_icedump_core2: MISSING DUMP — see $RUN/run.log"; tail -30 "$RUN/run.log"; exit 1
fi
