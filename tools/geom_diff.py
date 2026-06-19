#!/usr/bin/env python3
"""Compare two mesh-geometry dumps (FESOM2 oracle vs FESOM3) field-by-field for the
M1 geometry byte-gate (closes deferred M0.7).

Both files are written by the same self-describing format (FESOM2
src/fesom_geom_dump.F90 / FESOM3 src/infra/mod_geom_dump.F90):

  char8  "FGEOMDMP" | int32 nod2D,elem2D,edge2D,nl
  records: char24 name | int32 dtype(0=real64,1=int32) | int32 d1,d2 | d1*d2 values

Usage:  geom_diff.py <fesom2.bin> <fesom3.bin> [--signal 0.0]

Exit status 0 iff every common field matches within the signal threshold
(default 0 == byte-identical) and shapes/dims agree.
"""
import sys
import struct
import numpy as np


def read_dump(path):
    with open(path, "rb") as f:
        buf = f.read()
    if buf[:8] != b"FGEOMDMP":
        raise ValueError(f"{path}: bad magic {buf[:8]!r}")
    off = 8
    nod2D, elem2D, edge2D, nl = struct.unpack_from("<4i", buf, off)
    off += 16
    dims = dict(nod2D=nod2D, elem2D=elem2D, edge2D=edge2D, nl=nl)
    fields = {}
    n = len(buf)
    while off < n:
        name = buf[off:off + 24].decode("ascii", "replace").strip()
        off += 24
        dtype, d1, d2 = struct.unpack_from("<3i", buf, off)
        off += 12
        count = d1 * d2
        if dtype == 0:
            data = np.frombuffer(buf, dtype="<f8", count=count, offset=off)
            off += 8 * count
        elif dtype == 1:
            data = np.frombuffer(buf, dtype="<i4", count=count, offset=off)
            off += 4 * count
        else:
            raise ValueError(f"{path}: unknown dtype {dtype} for {name}")
        fields[name] = (dtype, d1, d2, data)
    return dims, fields


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    signal = 0.0
    for a in sys.argv[1:]:
        if a.startswith("--signal"):
            signal = float(a.split("=", 1)[1]) if "=" in a else 0.0
    if len(args) < 2:
        print(__doc__)
        sys.exit(2)
    p2, p3 = args[0], args[1]
    d2, f2 = read_dump(p2)
    d3, f3 = read_dump(p3)

    print(f"FESOM2: {p2}")
    print(f"FESOM3: {p3}")
    print(f"dims  FESOM2={d2}")
    print(f"dims  FESOM3={d3}")
    ok = True
    if d2 != d3:
        print("FAIL: mesh dims differ")
        ok = False

    names = list(f2.keys())
    extra3 = [k for k in f3 if k not in f2]
    print(f"\n{'field':20s} {'shape':>14s} {'dtype':>6s} {'max|delta|':>16s}  status")
    print("-" * 70)
    for name in names:
        if name not in f3:
            print(f"{name:20s} {'-':>14s} {'-':>6s} {'MISSING in F3':>16s}  FAIL")
            ok = False
            continue
        t2, a1, a2, v2 = f2[name]
        t3, b1, b2, v3 = f3[name]
        shape = f"{a1}x{a2}"
        if (t2, a1, a2) != (t3, b1, b2):
            print(f"{name:20s} {shape:>14s} {t2:>6d} shape/dtype F3={b1}x{b2}  FAIL")
            ok = False
            continue
        if v2.shape != v3.shape:
            print(f"{name:20s} {shape:>14s} {t2:>6d} {'len mismatch':>16s}  FAIL")
            ok = False
            continue
        if t2 == 0:
            d = np.abs(v2 - v3)
            mx = float(d.max()) if d.size else 0.0
            status = "PASS" if mx <= signal else "FAIL"
        else:
            nmis = int((v2 != v3).sum())
            mx = float(nmis)
            status = "PASS" if nmis == 0 else "FAIL"
        if status == "FAIL":
            ok = False
        print(f"{name:20s} {shape:>14s} {t2:>6d} {mx:>16.6e}  {status}")

    if extra3:
        print(f"\n(only in FESOM3, ignored: {extra3})")

    print("\n" + ("GEOMETRY BYTE-GATE: PASS (max|delta|=0)" if ok else
                  "GEOMETRY BYTE-GATE: FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
