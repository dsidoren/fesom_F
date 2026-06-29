#!/usr/bin/env bash
# Restart Stage 4 (Task 4.1) GATE: the IN-PROCESS write -> corrupt -> read round-trip that proves the
# READ path (zarr_read_chunk + decomp_gather + halo exchange) reproduces the WRITE bit-for-bit over the
# FULL local extent (owned + halo + eXDim) for the full oce+ice prognostic state incl. EVP sigma and the
# real(MP) mesh%hbar/hnode. fesom_restartroundtrip fills synthetic OWNED values, writes a checkpoint,
# halo-exchanges the originals (reference), CORRUPTS every live array (restart-owned region -> wild
# sentinel; un-owned tail -> 0), restart_reads, then compares each live array to its reference. The whole
# comparison is in-process Fortran: it prints per-field max|Δ| and 'ROUNDTRIP OK' / error-stops on any
# nonzero Δ. Runs AB_order=2 (26 fields) AND AB_order=3 (30 fields: +urhs_AB3/vrhs_AB3/temp_M2/salt_M2),
# at np 1 AND 2 (np=2 = real cross-rank gather + halo exchange — a wrong variant/gather shows up as a
# nonzero Δ in the halo/eXDim region).
#
#   tools/run_restartroundtrip.sh           # default: np 1 2
#   tools/run_restartroundtrip.sh 1         # just 1-rank
set -euo pipefail
F3="${F3:-/home/a/a270088/fesom3}"          # override to run from a git worktree checkout
BUILD="${BUILD:-$F3/build_intel_dp}"
PIMESH=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi
RUN="${RUN:-/scratch/a/a270088/restartroundtrip}"
NPS="${*:-1 2}"
source "$F3/env.sh" intel >/dev/null 2>&1

cmake "$BUILD" >/dev/null
cmake --build "$BUILD" --target fesom_restartroundtrip -j 4

export FESOM3_MESH_DIR="$PIMESH"
export FESOM3_CHUNK_HORIZ="${FESOM3_CHUNK_HORIZ:-1000}"   # small => multiple chunks + writer subset
ulimit -s unlimited
MPIRUN=(mpirun --mca pml ob1 --mca btl self,vader --mca btl_vader_single_copy_mechanism none --oversubscribe)

rc=0
for AB in 2 3; do
    [ "$AB" = 2 ] && NEXP=26 || NEXP=30
    export FESOM3_AB_ORDER="$AB"
    for NP in $NPS; do
        DIR="$RUN/ab${AB}_np$NP"
        rm -rf "$DIR"; mkdir -p "$DIR"
        export FESOM3_RESTARTRT_DIR="$DIR"
        echo "=== fesom_restartroundtrip AB_order=$AB np=$NP (expect $NEXP fields, max|Δ|=0) ==="
        if "${MPIRUN[@]}" -n "$NP" "$BUILD/bin/fesom_restartroundtrip" 2>&1 \
                | grep -E "fields=|global max|ROUNDTRIP"; then :; else rc=1; fi
    done
done

[ "$rc" = 0 ] && echo "run_restartroundtrip: GATE GREEN (np: $NPS)" || { echo "run_restartroundtrip: GATE FAILED"; exit 1; }
