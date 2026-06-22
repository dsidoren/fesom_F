#!/usr/bin/env bash
# Production-validation MULTI-RANK FREE-RUNNING lifecycle byte-gate (CORE2 dist_<NP>,
# unforced reduced-M2, dt=1800). Runs the REAL FESOM2 multi-step lifecycle (the instrumented
# fesom.x at NP ranks -> its built-in gid-keyed per-substep dump_shim) and FESOM3's MR
# lifecycle driver (fesom_lifecycle_mr -> mod_step_oce::step_oce per step through the optional
# partit), then compares the 13 NODE substeps x 5 probes over N steps PER RANK (dump_diff.py
# --glob, matched by GLOBAL probe id; SW_AB substep 2 = sw_alpha/beta, dead in M2, ignored).
#
# This closes the one combination the M2.12 single-step prescribe-and-stop gates never
# exercised: multi-rank x FREE-RUNNING x many steps — the AB2 velocity history + the cross-
# step halo state (w_e/w_i/eta carried forward) evolving in place. Both codes share dist_<NP>/
# myList, so each global probe id is owned by the same rank on both (the L8 same-partition
# rule). The IC (do_ic3d phc3.0) runs at multi-rank through the lifted extrap_nod3D.
#
#   tools/run_lifecycle_gate_multirank.sh [np] [nsteps] [run_dir]
set -euo pipefail
F3=/home/a/a270088/fesom3
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
NP="${1:-2}"
NSTEPS="${2:-3}"
RUN="${3:-/scratch/a/a270088/lifecycle_mr${NP}}"

echo "[1/3] FESOM2 oracle lifecycle ($NSTEPS steps, $NP-rank CORE2 dist_$NP, unforced)"
bash "$F3/tools/run_lifecycle_core2.sh" "$RUN" "$RUN/life_f2" "$NSTEPS" "$NP" | tail -3

echo "[2/3] FESOM3 MR lifecycle ($NSTEPS steps, $NP-rank CORE2 dist_$NP)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE"
export FESOM_DUMP_FILE="$RUN/life_f3" FESOM_DUMP_MAXSTEPS="$NSTEPS" FESOM3_NSTEPS="$NSTEPS"
ulimit -s unlimited
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n "$NP" \
    "$F3/build_intel_dp/bin/fesom_lifecycle_mr" > "$RUN/run_f3.log" 2>&1 || \
    { echo "  FESOM3 run failed"; tail -40 "$RUN/run_f3.log"; exit 1; }
grep -E 'nod2D|IC\(|step |done' "$RUN/run_f3.log" || true

echo "[3/3] compare per-rank (gid-keyed; 13 NODE substeps x 5 probes x $NSTEPS steps; SW_AB substep 2 ignored)"
python3 "$F3/tools/dump_diff.py" "$RUN/life_f2" "$RUN/life_f3" --glob --ignore-substep=2
