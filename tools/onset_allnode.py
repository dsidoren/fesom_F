#!/usr/bin/env python3
"""onset_allnode.py — localize the FIRST byte-divergence in windowed DUMP_ALL dumps.

The 5-probe gate can only see divergence once it reaches a probe; when all 5 probes
light up simultaneously at SSH_RHS that is the fingerprint of a GLOBAL scalar (the
freshwater `net` from integrate_nod_2D) spilling into water_flux at every node — it
does NOT localize the origin. This tool reads the all-owned-node DUMP_ALL window and
answers the decisive question: at the first diverging (step, substep), does ONE node
diverge (=> local origin, inspect that node's field) or ~ALL nodes (=> global spill,
instrument the freshwater components in oce_fluxes)?

Memory-frugal: loads the oracle side into a dict, then STREAMS the F3 side (no 2nd
full dict), so the resident set is ~one window's worth, not two.

Usage:  onset_allnode.py RUN_DIR [minstep] [maxstep]
        reads RUN_DIR/lifef_f2.* (oracle) and RUN_DIR/lifen_f3.* (F3)
"""
import sys, os, struct, glob
from collections import defaultdict

HDR = struct.Struct("<4i24s")        # step, substep, gid, nlevels, name[24]
SUB = {0:"INIT",1:"PRESS_BV",2:"SW_AB",3:"PGF",4:"MIXING",5:"VEL_RHS",6:"VISC_FILT",
       7:"IMPL_VISC",8:"SSH_RHS",9:"SSH_SOLVE",10:"UPD_VEL",11:"HBAR",12:"ETA_N",
       13:"ALE",14:"GM_BOLUS",15:"TRACERS",16:"THICKNESS"}

run = sys.argv[1]
lo  = int(sys.argv[2]) if len(sys.argv) > 2 else 0
hi  = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9

def load_dict(prefix):
    out = {}
    for path in glob.glob(prefix + ".*"):
        with open(path, "rb") as f:
            data = f.read()
        off, n = 0, len(data)
        while off < n:
            step, substep, gid, nlev, raw = HDR.unpack_from(data, off); off += HDR.size
            vals = struct.unpack_from("<%dd" % nlev, data, off); off += 8 * nlev
            if lo <= step <= hi:
                name = raw.split(b"\x00")[0].decode("ascii", "replace").strip()
                out[(step, substep, gid, name)] = vals
    return out

# Underflow floor: a level-pair where BOTH magnitudes are below this is denormal/underflow
# noise (e.g. a snow thickness of ~1e-308 m = physically zero) that gets absorbed the moment it
# meets any normal-magnitude quantity; it cannot drive a real divergence, so it is ignored.
NEG = 1e-100

def maxd(x, y):
    if len(x) != len(y):
        return float("inf")
    m = 0.0
    for p, q in zip(x, y):
        if max(abs(p), abs(q)) < NEG:   # both denormal/underflow -> not a physical divergence
            continue
        d = abs(p - q)
        if d > m:
            m = d
    return m

print("loading oracle (lifef_f2) ...", flush=True)
oracle = load_dict(run + "/lifef_f2")
print("  oracle records in window [%d,%d]: %d" % (lo, hi, len(oracle)), flush=True)

# Stream F3, compare against oracle, collect divergences only.
div = []                       # (step, substep, gid, name, d, f2_0, f3_0)
f3_count = 0
missing_in_oracle = 0
print("streaming F3 (lifen_f3) ...", flush=True)
for path in glob.glob(run + "/lifen_f3.*"):
    with open(path, "rb") as f:
        data = f.read()
    off, n = 0, len(data)
    while off < n:
        step, substep, gid, nlev, raw = HDR.unpack_from(data, off); off += HDR.size
        vals = struct.unpack_from("<%dd" % nlev, data, off); off += 8 * nlev
        if not (lo <= step <= hi):
            continue
        name = raw.split(b"\x00")[0].decode("ascii", "replace").strip()
        f3_count += 1
        k = (step, substep, gid, name)
        ov = oracle.get(k)
        if ov is None:
            missing_in_oracle += 1
            continue
        d = maxd(ov, vals)
        if d > 0.0:
            div.append((step, substep, gid, name, d, ov[0], vals[0]))
print("  f3 records: %d   (keys missing in oracle: %d)" % (f3_count, missing_in_oracle), flush=True)

if not div:
    print("\nNO divergence anywhere in window [%d,%d] — onset is AFTER step %d." % (lo, hi, hi))
    sys.exit(0)

div.sort(key=lambda r: (r[0], r[1], r[2]))   # step, substep, gid (execution order)
steps_div = sorted(set(d[0] for d in div))
onset = steps_div[0]
print("\n==================== DIVERGENCE LOCALIZATION ====================")
print("steps in window with ANY diverging node:", steps_div)
print("ONSET step = %d" % onset)

prev = onset - 1
if prev < lo:
    print("WARNING: onset is the window's first step (%d) — cannot confirm the prior step is\n"
          "         clean at all nodes. Widen the window (lower FESOM_DUMP_MINSTEP) to be sure." % onset)
else:
    clean_prev = not any(d[0] == prev for d in div)
    print("step %d (one before onset) clean at ALL owned nodes? %s"
          % (prev, "YES" if clean_prev else "NO — divergence pre-dates the onset!"))

# Per-(substep,field) diverging-node count at the onset step, with lowest-gid example.
cnt = defaultdict(int); ex = {}
for step, substep, gid, name, d, f2, f3 in div:
    if step == onset:
        cnt[(substep, name)] += 1
        if (substep, name) not in ex or gid < ex[(substep, name)][0]:
            ex[(substep, name)] = (gid, d, f2, f3)
first_ss = min(ss for (ss, nm) in cnt)
print("\nper-(substep,field) diverging-node COUNT at onset step %d:" % onset)
for (ss, nm) in sorted(cnt):
    gid, d, f2, f3 = ex[(ss, nm)]
    flag = "   <== FIRST substep" if ss == first_ss else ""
    print("  sub %2d %-10s %-18s ndiv=%7d  lowgid=%7d |d|=%.4e F2=%.12e F3=%.12e%s"
          % (ss, SUB.get(ss, "?"), nm, cnt[(ss, nm)], gid, d, f2, f3, flag))

# Verdict: how many nodes diverge at the FIRST substep of the onset step?
first_counts = [cnt[(ss, nm)] for (ss, nm) in cnt if ss == first_ss]
mincount = min(first_counts)
print("\n---------------------------------------------------------------")
print("FIRST diverging substep at onset: %d (%s); min node-count among its fields = %d"
      % (first_ss, SUB.get(first_ss, "?"), mincount))
if mincount <= 8:
    fields = [nm for (ss, nm) in cnt if ss == first_ss and cnt[(ss, nm)] == mincount]
    print("VERDICT: *** LOCAL ORIGIN *** — only %d node(s) in field(s) %s." % (mincount, fields))
    print("         Inspect that node's FP-form for this field (likely a wet/dry, cavity,")
    print("         boundary, or ice-edge node with a divergent branch).")
    # dump the exact origin record(s)
    for (ss, nm) in sorted(cnt):
        if ss == first_ss and cnt[(ss, nm)] == mincount:
            for step, substep, gid, name, d, f2, f3 in div:
                if step == onset and substep == ss and name == nm:
                    print("    ORIGIN: step %d sub %d %-16s gid %7d |d|=%.4e F2=%.15e F3=%.15e"
                          % (step, substep, SUB.get(substep), gid, d, f2, f3))
else:
    print("VERDICT: *** GLOBAL SPILL *** — %d nodes diverge at the first substep." % mincount)
    print("         The origin is an upstream GLOBAL scalar (freshwater `net` via")
    print("         integrate_nod_2D in oce_fluxes), not a single dumped node field.")
    print("         Next: instrument the per-node freshwater-flux integrand (flux/relax_salt/")
    print("         water_flux) in oce_fluxes on BOTH codes and re-run this window.")
print("================================================================")
