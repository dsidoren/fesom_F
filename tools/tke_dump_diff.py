#!/usr/bin/env python3
"""tke_dump_diff.py — compare two sets of TKE text dumps (tke_dump_mod format).

The oracle (gen_modules_cvmix_tke.F90 tke_dump_mod) and FESOM3 (fesom_tkereplay /
fesom_tkedump) both emit per-step per-gid rows:

    # step=1 tag=vshear2 rank=0 N=<nnod> ncomp=<nl>
    <gid> <v1> <v2> ... <vN>     (each value es24.16 — exact double round-trip)

For each (step, tag) it parses A/tke_dump_s<step>_<tag>_rank<R>.txt and the B counterpart
into {gid: tuple(floats)}, computes max|Δ| over all gid×component, and reports. Fortran↔
Fortran target is max|Δ|=0 EXACTLY (no tolerance floor). Exit 0 iff every tag matches at
or below --threshold (default 0.0); else exit 1 with the worst (step, tag, gid, comp).

    tke_dump_diff.py <dir_a> <dir_b> [--steps 1,2,3] [--tags tke,tkeav,...]
                     [--ranks 0] [--threshold 0.0]
"""
import argparse
import os
import sys

# M7a-1 gate output tags (column outputs + the Part-localizing intermediates).
DEFAULT_TAGS = "tke,tkeav,tkekv,lmix,pr,tbpr,tspr,tdif,tdis,twin,tiwf,tbck,ttot"


def parse_file(path):
    """Return {gid: tuple(float, ...)} from a tke_dump text file (header skipped)."""
    out = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            gid = int(parts[0])
            out[gid] = tuple(float(x) for x in parts[1:])
    return out


def compare_pair(pa, pb):
    """max|Δ|, worst (gid, comp, a, b), n_mismatch, n_common over two parsed dicts."""
    a = parse_file(pa)
    b = parse_file(pb)
    common = a.keys() & b.keys()
    worst = 0.0
    worst_loc = None
    n_mismatch = 0
    for gid in common:
        va, vb = a[gid], b[gid]
        for c, (x, y) in enumerate(zip(va, vb)):
            d = abs(x - y)
            if d > 0.0:
                n_mismatch += 1
            if d > worst:
                worst = d
                worst_loc = (gid, c, x, y)
    only_a = a.keys() - b.keys()
    only_b = b.keys() - a.keys()
    return worst, worst_loc, n_mismatch, len(common), only_a, only_b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir_a")
    ap.add_argument("dir_b")
    ap.add_argument("--steps", default="1,2,3")
    ap.add_argument("--tags", default=DEFAULT_TAGS)
    ap.add_argument("--ranks", default="0")
    ap.add_argument("--threshold", type=float, default=0.0)
    args = ap.parse_args()

    steps = [int(s) for s in args.steps.split(",") if s != ""]
    tags = [t for t in args.tags.split(",") if t != ""]
    ranks = [int(r) for r in args.ranks.split(",") if r != ""]

    overall_worst = 0.0
    overall_loc = None
    failed = []
    missing = []
    npairs = 0

    for step in steps:
        for tag in tags:
            for rank in ranks:
                fn = f"tke_dump_s{step}_{tag}_rank{rank}.txt"
                pa = os.path.join(args.dir_a, fn)
                pb = os.path.join(args.dir_b, fn)
                if not (os.path.exists(pa) and os.path.exists(pb)):
                    missing.append(fn)
                    continue
                npairs += 1
                worst, loc, nmis, ncom, oa, ob = compare_pair(pa, pb)
                status = "OK  " if worst <= args.threshold else "FAIL"
                extra = ""
                if loc is not None and worst > args.threshold:
                    gid, c, x, y = loc
                    extra = f"  worst gid={gid} comp={c} a={x!r} b={y!r}"
                if oa or ob:
                    extra += f"  [onlyA={len(oa)} onlyB={len(ob)}]"
                print(f"  [{status}] s{step} {tag:<10} max|Δ|={worst:.3e} "
                      f"nmis={nmis}/{ncom*max(1,1)}{extra}")
                if worst > overall_worst:
                    overall_worst = worst
                    overall_loc = (step, tag, loc)
                if worst > args.threshold:
                    failed.append((step, tag, worst))

    print("-" * 72)
    if missing:
        print(f"tke_dump_diff: {len(missing)} pair(s) MISSING (e.g. {missing[0]})")
    if npairs == 0:
        print("tke_dump_diff: NO pairs compared — check dirs/steps/tags")
        return 1
    if failed:
        s, t, loc = overall_loc
        print(f"tke_dump_diff: DIVERGENCE — {len(failed)} tag(s) > {args.threshold:g}; "
              f"worst max|Δ|={overall_worst:.6e} at step {s} tag {t}")
        return 1
    print(f"tke_dump_diff: MATCH — {npairs} pairs, max|Δ|={overall_worst:.3e} "
          f"(threshold {args.threshold:g})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
