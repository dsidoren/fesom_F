#!/usr/bin/env bash
# M2.11c-2 oracle: REAL FESOM2 FORCED multi-step lifecycle on CORE2 1-rank (use_ice=.true.,
# CORE2 NCAR forcing at the 1948 stubs + pool runoff/SSS, the ice EVP + oce_fluxes producing
# the air-sea fluxes). Emits BOTH the built-in per-substep dump_shim (FESOM_DUMP_FILE, the 13
# NODE substeps) AND the NEW per-step full-field flux dump (FESOM_FLUX_DUMP, fesom_flux_dump.F90:
# heat_flux/water_flux/virtual_salt/relax_salt + stress_surf, written BEFORE oce_timestep_ale) so
# FESOM3 can PRESCRIBE the M3 air-sea gap. Reduced-M2 dynamics (linfs/PP/no-GM/no-Redi/opt_visc=7).
#
# Calendar: CORE (noleap) forcing REQUIRES include_fleapyear=.false. (else FESOM2 stops with a
# calendar-consistency error). use_sw_pene: reduced to .false. by default (heat_flux = raw obudget,
# cal_shortwave_rad skipped); SW_PENE=1 keeps work_core .true. (M5c). Needs the np=1 next_io_rank
# fix + the output()/write_initial_conditions() npes==1 early-returns.
#
#   tools/run_lifecycle_forced_core2.sh [run_dir] [dump_prefix] [flux_prefix] [nsteps] [atmflux] [whichEVP] [np]
# [np] (default 1) selects CORE2 dist_<np> (FESOM2 auto-reads dist_<npes> from the mesh path);
# np>1 is the M3f-4 multi-rank oracle — the dump_shim is gid-keyed -> per-rank <DUMP>.<mype5>,
# and the 1-rank-only flux/atmflux dumps skip harmlessly (npes/=1), so the per-rank node-substep
# dump is the gate target. The 1-rank io workarounds are inactive at npes>1 (normal io_gather).
set -euo pipefail
F2=/home/a/a270088/port2/fesom2
STUB="$F2/test/input/global"
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
RUN="${1:-/scratch/a/a270088/lifecycle_forced_core2}"
DUMP="${2:-$RUN/lifef_f2}"
FLUX="${3:-$RUN/flux_f2}"
NSTEPS="${4:-3}"
ATMFLUX="${5:-$RUN/atmflux_f2}"   # M3f: per-step post-bulk atmospheric forcing (fesom_atmflux_dump)
WHICHEVP="${6:-0}"                # M3f: ice EVP solver (0=std EVP, 1=mEVP) — must match FESOM3
NP="${7:-1}"                      # M3f-4: number of ranks (CORE2 dist_<NP>)
# M4e: FER_GM=1 / REDI=1 keep the work_core Fer_GM/Redi=.true. (GM bolus / Redi isopycnal
# diffusion ON) instead of the reduced-M2 sed-off. Default 0 = the proven GM/Redi-off gate.
FER_GM="${FER_GM:-0}"
REDI="${REDI:-0}"
# M5c: MIX_KPP=1 keeps work_core mix_scheme='KPP' (the production boundary-layer scheme) instead
# of the reduced sed-to-PP; SW_PENE=1 keeps use_sw_pene=.true. (the shortwave-penetration tracer
# term); KPP_NONLCL=1 injects use_kpp_nonlclflx=.true. into &tracer_phys (the ghats nonlocal flux,
# DEAD in work_core by default). Defaults 0 reproduce the M3f/M4e reduced (PP, no-sw_pene) gate.
MIX_KPP="${MIX_KPP:-0}"
SW_PENE="${SW_PENE:-0}"
KPP_NONLCL="${KPP_NONLCL:-0}"
# M6a-3: WHICH_ALE selects the vertical coordinate (default linfs reduced gate). zstar keeps
# the production free surface (Shchepetkin PGF + stretch) + the REAL freshwater flux path.
WHICH_ALE="${WHICH_ALE:-linfs}"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1
rm -rf "$RUN"; mkdir -p "$RUN"
cp "$F2"/work_core/namelist.* "$RUN"/
cp "$F2"/work_core/namelist.forcing.CORE2 "$RUN"/namelist.forcing
sed -i "s/^whichEVP *=.*/whichEVP = ${WHICHEVP}/" "$RUN"/namelist.ice   # M3f: match FESOM3 EVP variant
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x
printf '0 1 1948\n0 1 1948\n' > "$RUN"/fesom.clock

python3 - "$RUN" "$STUB" "$POOL" "$NSTEPS" "$FER_GM" "$REDI" "$MIX_KPP" "$SW_PENE" "$KPP_NONLCL" "$WHICH_ALE" <<'PY'
import re,sys
run,stub,pool,nsteps,fer_gm,redi,mix_kpp,sw_pene,kpp_nonlcl,which_ale=sys.argv[1:11]
p=run+'/namelist.config'; s=open(p).read()
s=re.sub(r"ResultPath\s*=\s*'[^']*'","ResultPath       = './'",s,1)
s=re.sub(r"which_ALE\s*=\s*'zlevel'",f"which_ALE          = '{which_ale}'",s,1)
s=re.sub(r"yearnew\s*=\s*1958","yearnew = 1948",s,1)
s=re.sub(r"run_length\s*=\s*\d+",f"run_length        = {nsteps}",s,1)
s=re.sub(r"include_fleapyear\s*=\s*\.true\.","include_fleapyear = .false.",s,1)   # CORE noleap
# M5c: SW_PENE=1 keeps work_core use_sw_pene=.true. (cal_shortwave_rad + the sw_3d tracer term);
# else reduce it off (the M3f/M4e gate, no sw_3d consumer).
if sw_pene!='1':
    s=re.sub(r"use_sw_pene\s*=\s*\.true\.","use_sw_pene              = .false.",s,1)
open(p,'w').write(s)
p=run+'/namelist.oce'; s=open(p).read()
# M5c: MIX_KPP=1 keeps work_core mix_scheme='KPP'; else reduce to PP (the M2-M4 reduced core).
if mix_kpp!='1':
    s=re.sub(r"mix_scheme\s*=\s*'KPP'","mix_scheme         = 'PP'",s,1)
# M4e: keep work_core Fer_GM/Redi=.true. when FER_GM=1 / REDI=1; else reduce them off.
if fer_gm!='1':
    s=re.sub(r"Fer_GM\s*=\s*\.true\.","Fer_GM             = .false.",s,1)
if redi!='1':
    s=re.sub(r"Redi\s*=\s*\.true\.","Redi               = .false.",s,1)
open(p,'w').write(s)
# M5c: KPP_NONLCL=1 injects use_kpp_nonlclflx=.true. into &tracer_phys (absent in work_core ->
# default .false.; this turns ON the ghats nonlocal counter-gradient flux for the gate variant).
if kpp_nonlcl=='1':
    p=run+'/namelist.tra'; s=open(p).read()
    s=re.sub(r"(&tracer_phys\s*\n)", r"\1use_kpp_nonlclflx  = .true.\n", s, 1)
    open(p,'w').write(s)
p=run+'/namelist.forcing'; s=open(p).read()
s=s.replace("FORCING/CORE2/", stub+"/")                                          # atm stubs absolute
s=re.sub(r"nm_runoff_file\s*=\s*'[^']*'",   f"nm_runoff_file ='{pool}/CORE2_runoff.nc'",s,1)
s=re.sub(r"nm_sss_data_file\s*=\s*'[^']*'", f"nm_sss_data_file ='{pool}/PHC2_salx.nc'",s,1)
open(p,'w').write(s)
PY

cd "$RUN"
export FESOM_DUMP_FILE="$DUMP" FESOM_DUMP_MAXSTEPS="$NSTEPS" FESOM_FLUX_DUMP="$FLUX"
export FESOM_ATMFLUX_DUMP="$ATMFLUX"   # M3f native-flux gate (harmless extra dump for the M2.11c-2 gate)
ulimit -s unlimited
echo "run_lifecycle_forced_core2: $NP rank(s), $NSTEPS steps (FORCED use_ice) -> ${DUMP}.*"
timeout 600 mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n "$NP" ./fesom.x > "$RUN/run.log" 2>&1 || true
if ls "${DUMP}".* >/dev/null 2>&1; then
    echo "run_lifecycle_forced_core2: done -> dump $(stat -c%s "${DUMP}".00000)B ($(ls "${DUMP}".* | wc -l) rank file(s))"
    grep -E 'FESOM Run|FDBG step|forcing init' "$RUN/run.log" | head
else
    echo "run_lifecycle_forced_core2: MISSING DUMP — see $RUN/run.log"; tail -30 "$RUN/run.log"; exit 1
fi
