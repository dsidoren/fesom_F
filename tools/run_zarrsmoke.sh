#!/usr/bin/env bash
# M9 Task 0.2 GATE: round-trip the hand-rolled Zarr v2 writer (mod_io_zarr). fesom_zarrsmoke writes
# toy 1-D f8 / 2-D f4 / 2-D i4 arrays (non-dividing chunks => partial-chunk padding; distinct per-cell
# values => catches a C-order transpose bug); tools/zarr_diff.py --roundtrip reads them back with
# zarr AND xarray and asserts max|Δ|=0 (f8/i4) / exact (f4). Self-contained: no mesh/MPI/oracle.
#
#   tools/run_zarrsmoke.sh
set -euo pipefail
F3=/home/a/a270088/fesom3
BUILD="${BUILD:-$F3/build_intel_dp}"
RUN="${RUN:-/scratch/a/a270088/zarrsmoke}"
PY=/work/ab0995/a270088/mambaforge/bin/python3
mkdir -p "$RUN"
source "$F3/env.sh" intel >/dev/null 2>&1

# Reconfigure (re-glob new src/io + src/drivers files) then incremental build. The explicit
# reconfigure is needed the first time a brand-new target/file appears (CONFIGURE_DEPENDS alone
# can't make a target that doesn't exist in the current cache yet).
cmake "$BUILD" >/dev/null
cmake --build "$BUILD" --target fesom_zarrsmoke -j 4

rm -rf "$RUN/zarrsmoke.zarr" "$RUN/zarrsmoke_lz4.zarr"
FESOM3_ZARR_OUT="$RUN" "$BUILD/bin/fesom_zarrsmoke"

"$PY" "$F3/tools/zarr_diff.py" --roundtrip "$RUN/zarrsmoke.zarr"
"$PY" "$F3/tools/zarr_diff.py" --lz4       "$RUN/zarrsmoke_lz4.zarr"
echo "run_zarrsmoke: GATE GREEN"
