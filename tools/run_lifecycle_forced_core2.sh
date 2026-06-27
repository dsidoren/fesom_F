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
# M8b: FORCING_OVERRIDE replaces the namelist `FORCING/CORE2/` atm prefix (default = the day-1
# STUB, so existing gates are byte-unchanged). The long gate points it at the full-year CORE2
# pool so the run can cross 6-hourly wind / daily radiation record boundaries. START_CLOCK
# (default the proven cold start) seeds fesom.clock so a non-day-1 start can be byte-gated too.
FORCING_OVERRIDE="${FORCING_OVERRIDE:-$STUB}"
START_CLOCK="${START_CLOCK:-0 1 1948}"
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
# M7c: MIX_TKE=1 sets mix_scheme='cvmix_TKE' (prognostic TKE producer) + supplies namelist.cvmix
# &param_tke (tke_cd=3.75). Pair with SW_PENE=1 (work_*_tke has use_sw_pene=.true.). Default 0.
MIX_TKE="${MIX_TKE:-0}"
# M6a-3: WHICH_ALE selects the vertical coordinate (default linfs reduced gate). zstar keeps
# the production free surface (Shchepetkin PGF + stretch) + the REAL freshwater flux path.
WHICH_ALE="${WHICH_ALE:-linfs}"
# M8c: FORCING_SET=JRA55 runs the production JRA55-do atmosphere (gregorian + leap, work_core
# default namelist.forcing, start 1958) instead of the CORE2 NCAR substitute. CHL_SWEENEY=1 keeps
# chl='Sweeney' (production monthly climatology); default 0 => chl='None' const 0.1 (matches the F3
# constant chl — isolates the atmosphere read for the step-1 gate).
FORCING_SET="${FORCING_SET:-CORE2}"
CHL_SWEENEY="${CHL_SWEENEY:-0}"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1
rm -rf "$RUN"; mkdir -p "$RUN"
cp "$F2"/work_core/namelist.* "$RUN"/
if [ "$FORCING_SET" = JRA55 ]; then
    cp "$F2"/work_core/namelist.forcing "$RUN"/namelist.forcing       # JRA55-do default (uas/vas/huss/...)
else
    cp "$F2"/work_core/namelist.forcing.CORE2 "$RUN"/namelist.forcing # CORE2 NCAR substitute
fi
# M7c: cvmix_TKE needs the &param_tke namelist group (tke_cd=3.75, the namelist-over-codedefault).
[ "$MIX_TKE" = 1 ] && cp "$F2"/work_tke_dump/namelist.cvmix "$RUN"/
sed -i "s/^whichEVP *=.*/whichEVP = ${WHICHEVP}/" "$RUN"/namelist.ice   # M3f: match FESOM3 EVP variant
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x
printf '%s\n%s\n' "$START_CLOCK" "$START_CLOCK" > "$RUN"/fesom.clock   # M8b: cold start at START_CLOCK

python3 - "$RUN" "$FORCING_OVERRIDE" "$POOL" "$NSTEPS" "$FER_GM" "$REDI" "$MIX_KPP" "$SW_PENE" "$KPP_NONLCL" "$WHICH_ALE" "$MIX_TKE" "$FORCING_SET" "$CHL_SWEENEY" <<'PY'
import re,sys
run,stub,pool,nsteps,fer_gm,redi,mix_kpp,sw_pene,kpp_nonlcl,which_ale,mix_tke,forcing_set,chl_sweeney=sys.argv[1:14]
jra=(forcing_set=='JRA55')
p=run+'/namelist.config'; s=open(p).read()
s=re.sub(r"ResultPath\s*=\s*'[^']*'","ResultPath       = './'",s,1)
s=re.sub(r"which_ALE\s*=\s*'zlevel'",f"which_ALE          = '{which_ale}'",s,1)
# M8c: JRA55 keeps the work_core production start (1958) + leap calendar; CORE2 substitutes noleap/1948.
if not jra:
    s=re.sub(r"yearnew\s*=\s*1958","yearnew = 1948",s,1)
s=re.sub(r"run_length\s*=\s*\d+",f"run_length        = {nsteps}",s,1)
if not jra:
    s=re.sub(r"include_fleapyear\s*=\s*\.true\.","include_fleapyear = .false.",s,1)   # CORE noleap
# M5c: SW_PENE=1 keeps work_core use_sw_pene=.true. (cal_shortwave_rad + the sw_3d tracer term);
# else reduce it off (the M3f/M4e gate, no sw_3d consumer).
if sw_pene!='1':
    s=re.sub(r"use_sw_pene\s*=\s*\.true\.","use_sw_pene              = .false.",s,1)
open(p,'w').write(s)
p=run+'/namelist.oce'; s=open(p).read()
# M7c: MIX_TKE=1 -> mix_scheme='cvmix_TKE'. Else M5c: keep 'KPP' when MIX_KPP=1, else reduce to 'PP'.
if mix_tke=='1':
    s=re.sub(r"mix_scheme\s*=\s*'KPP'","mix_scheme         = 'cvmix_TKE'",s,1)
elif mix_kpp!='1':
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
if not jra:
    s=s.replace("FORCING/CORE2/", stub+"/")              # M8b: atm dir absolute (stub or CORE2 pool)
# M8c: JRA55 default is chl='Sweeney' (monthly clim); CHL_SWEENEY!=1 reduces it to 'None'/const 0.1
# (matches the F3 constant chl — isolates the atmosphere read for the step-1 gate).
if jra and chl_sweeney!='1':
    s=re.sub(r"chl_data_source\s*=\s*'Sweeney'","chl_data_source    ='None'",s,1)
s=re.sub(r"nm_runoff_file\s*=\s*'[^']*'",   f"nm_runoff_file ='{pool}/CORE2_runoff.nc'",s,1)
s=re.sub(r"nm_sss_data_file\s*=\s*'[^']*'", f"nm_sss_data_file ='{pool}/PHC2_salx.nc'",s,1)
open(p,'w').write(s)
PY

cd "$RUN"
export FESOM_DUMP_FILE="$DUMP" FESOM_DUMP_MAXSTEPS="$NSTEPS" FESOM_FLUX_DUMP="$FLUX"
export FESOM_ATMFLUX_DUMP="$ATMFLUX"   # M3f native-flux gate (harmless extra dump for the M2.11c-2 gate)
# M8b: login dump writer needs a big stack (L20 unlimited); Levante compute nodes reject unlimited
# so the dist_864 batch uses the job_2yr_864 value (204800 KB).
if [ -n "${SLURM_JOB_ID:-}" ]; then ulimit -s 204800; else ulimit -s unlimited; fi
echo "run_lifecycle_forced_core2: $NP rank(s), $NSTEPS steps (FORCED use_ice) -> ${DUMP}.*"
# M8b: inside a SLURM allocation (dist_864 acceptance gate) launch with srun; on the login node
# keep the proven mpirun line verbatim (byte-unchanged for the dist_2/dist_8 gates).
if [ -n "${SLURM_JOB_ID:-}" ]; then
    timeout 1200 srun -l -n "$NP" ./fesom.x > "$RUN/run.log" 2>&1 || true
else
    timeout 600 mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n "$NP" ./fesom.x > "$RUN/run.log" 2>&1 || true
fi
if ls "${DUMP}".* >/dev/null 2>&1; then
    echo "run_lifecycle_forced_core2: done -> dump $(stat -c%s "${DUMP}".00000)B ($(ls "${DUMP}".* | wc -l) rank file(s))"
    # M8b: `| head` closes the pipe after 10 lines; at long runs (48-step M8b gate) grep is still
    # writing -> SIGPIPE (141) -> with `set -e`+pipefail this aborted the whole gate before [2/3].
    # Latent since the legacy 3-5 step gates never produced >10 matches. `|| true` makes it benign.
    grep -E 'FESOM Run|FDBG step|forcing init' "$RUN/run.log" | head || true
else
    echo "run_lifecycle_forced_core2: MISSING DUMP — see $RUN/run.log"; tail -30 "$RUN/run.log"; exit 1
fi
