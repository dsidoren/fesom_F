#!/usr/bin/env bash
# M2.11c-2 oracle: REAL FESOM2 FORCED multi-step lifecycle on CORE2 1-rank (use_ice=.true.,
# CORE2 NCAR forcing at the 1948 stubs + pool runoff/SSS, the ice EVP + oce_fluxes producing
# the air-sea fluxes). Emits BOTH the built-in per-substep dump_shim (FESOM_DUMP_FILE, the 13
# NODE substeps) AND the NEW per-step full-field flux dump (FESOM_FLUX_DUMP, fesom_flux_dump.F90:
# heat_flux/water_flux/virtual_salt/relax_salt + stress_surf, written BEFORE oce_timestep_ale) so
# FESOM3 can PRESCRIBE the M3 air-sea gap. Reduced-M2 dynamics (linfs/PP/no-GM/no-Redi/opt_visc=7).
#
# Calendar: CORE (noleap) forcing REQUIRES include_fleapyear=.false. (else FESOM2 stops with a
# calendar-consistency error). use_sw_pene=.false.: matches the ported step_oce (no sw_3d term);
# cal_shortwave_rad is then skipped so heat_flux is the raw obudget value. Needs the np=1
# next_io_rank fix + the output()/write_initial_conditions() npes==1 early-returns.
#
#   tools/run_lifecycle_forced_core2.sh [run_dir] [dump_prefix] [flux_prefix] [nsteps]
set -euo pipefail
F2=/home/a/a270088/port2/fesom2
STUB="$F2/test/input/global"
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
RUN="${1:-/scratch/a/a270088/lifecycle_forced_core2}"
DUMP="${2:-$RUN/lifef_f2}"
FLUX="${3:-$RUN/flux_f2}"
NSTEPS="${4:-3}"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1
rm -rf "$RUN"; mkdir -p "$RUN"
cp "$F2"/work_core/namelist.* "$RUN"/
cp "$F2"/work_core/namelist.forcing.CORE2 "$RUN"/namelist.forcing
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x
printf '0 1 1948\n0 1 1948\n' > "$RUN"/fesom.clock

python3 - "$RUN" "$STUB" "$POOL" "$NSTEPS" <<'PY'
import re,sys
run,stub,pool,nsteps=sys.argv[1:5]
p=run+'/namelist.config'; s=open(p).read()
s=re.sub(r"ResultPath\s*=\s*'[^']*'","ResultPath       = './'",s,1)
s=re.sub(r"which_ALE\s*=\s*'zlevel'","which_ALE          = 'linfs'",s,1)
s=re.sub(r"yearnew\s*=\s*1958","yearnew = 1948",s,1)
s=re.sub(r"run_length\s*=\s*\d+",f"run_length        = {nsteps}",s,1)
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
export FESOM_DUMP_FILE="$DUMP" FESOM_DUMP_MAXSTEPS="$NSTEPS" FESOM_FLUX_DUMP="$FLUX"
ulimit -s unlimited
echo "run_lifecycle_forced_core2: 1 rank, $NSTEPS steps (FORCED use_ice) -> ${DUMP}.* + ${FLUX}.*"
timeout 600 mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 ./fesom.x > "$RUN/run.log" 2>&1 || true
if ls "${DUMP}".* >/dev/null 2>&1 && ls "${FLUX}".* >/dev/null 2>&1; then
    echo "run_lifecycle_forced_core2: done -> dump $(stat -c%s "${DUMP}".00000)B, flux $(stat -c%s "${FLUX}".00000)B"
    grep -E 'FESOM Run|FDBG step|forcing init' "$RUN/run.log" | head
else
    echo "run_lifecycle_forced_core2: MISSING DUMP — see $RUN/run.log"; tail -30 "$RUN/run.log"; exit 1
fi
