#!/usr/bin/env python3
"""Extract every fw_/substep record for one gid across a step window from F2+F3 dumps.
Usage: extract_gid.py RUN_DIR GID [minstep] [maxstep]
"""
import sys, struct, glob
HDR = struct.Struct("<4i24s")
run, gid = sys.argv[1], int(sys.argv[2])
lo = int(sys.argv[3]) if len(sys.argv) > 3 else 0
hi = int(sys.argv[4]) if len(sys.argv) > 4 else 10**9

def grab(prefix):
    out = {}
    for path in glob.glob(prefix + ".*"):
        with open(path, "rb") as f:
            data = f.read()
        off, n = 0, len(data)
        while off < n:
            step, sub, g, nlev, raw = HDR.unpack_from(data, off); off += HDR.size
            vals = struct.unpack_from("<%dd" % nlev, data, off); off += 8 * nlev
            if g == gid and lo <= step <= hi:
                name = raw.split(b"\x00")[0].decode("ascii", "replace").strip()
                out[(step, sub, name)] = vals
    return out

a = grab(run + "/lifef_f2")
b = grab(run + "/lifen_f3")
keys = sorted(set(a) | set(b))
print("gid %d  window [%d,%d]   (F2=oracle, F3)" % (gid, lo, hi))
print("%-5s %-4s %-12s %-26s %-26s %-10s" % ("step", "sub", "field", "F2", "F3", "|d|"))
for k in keys:
    va = a.get(k); vb = b.get(k)
    if va is None or vb is None:
        print("  %-5d %-4d %-12s  ONLY IN %s" % (k[0], k[1], k[2], "F2" if vb is None else "F3"))
        continue
    d = max((abs(x-y) for x, y in zip(va, vb)), default=0.0) if len(va) == len(vb) else float("inf")
    star = "  <<<" if d > 0 else ""
    print("%-5d %-4d %-12s %-26.18e %-26.18e %-10.3e%s" % (k[0], k[1], k[2], va[0], vb[0], d, star))
