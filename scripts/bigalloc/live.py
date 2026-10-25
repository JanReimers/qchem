"""The LIVE SET AT PEAK of a bigalloc log, grouped by the first qchem source frames (addr2line -f -C -i, inlines expanded).
usage: live.py TRACE [DEPTH=4] [TOP=25]"""
import re, sys, subprocess, collections, functools

@functools.lru_cache(maxsize=None)
def resolve(lib, off):
    r = subprocess.run(["addr2line", "-f", "-C", "-i", "-e", lib, off], capture_output=True, text=True).stdout.split("\n")
    out = []
    for fn, loc in zip(r[0::2], r[1::2]):
        loc = re.sub(r" \(discriminator.*", "", loc)
        if "/src/" not in loc or "qchem" not in loc: continue
        out.append("[%s] %s" % (loc.split("/src/", 1)[1], re.sub(r"\(.*", "", fn)[-60:]))
    return out

@functools.lru_cache(maxsize=None)
def symtab(lib):
    d = {}
    for l in subprocess.run(["nm", "--defined-only", lib], capture_output=True, text=True).stdout.split("\n"):
        x = l.split()
        if len(x) == 3: d[x[2]] = int(x[0], 16)
    return d

NOISE = re.compile(r"Math/Blaze\.C|Imp/Factory\.C|HamiltonianTerm\.C")
def frames(bt, depth):
    out = []
    for f in bt:
        m = re.match(r"(\S+?)\(([^+)]*)\+(0x[0-9a-f]+)\)", f)
        if not m or "qchem" not in m[1]: continue
        lib, sym, off = m[1], m[2], int(m[3], 16)
        if sym:
            base = symtab(lib).get(sym)
            if base is None: continue
            off += base
        for x in resolve(lib, hex(max(off - 1, 0))):
            if NOISE.search(x): continue
            if not out or out[-1] != x: out.append(x)
        if len(out) >= depth: break
    return "\n          <- ".join(out[:depth]) or bt[0][:110]

def peak_set(path):
    live, cur, peak, peakset = {}, 0, 0, {}
    lines = open(path, errors="replace").read().split("\n")
    i = 0
    while i < len(lines):
        l = lines[i]
        if l.startswith("A "):
            _, p, n = l.split(); n = int(n); bt = []
            i += 1
            while i < len(lines) and lines[i] != "E": bt.append(lines[i]); i += 1
            live[p] = (n, tuple(bt)); cur += n
            if cur > peak: peak, peakset = cur, dict(live)
        elif l.startswith("F "):
            p = l.split()[1]
            if p in live: cur -= live.pop(p)[0]
        i += 1
    return peak, peakset

if __name__ == "__main__":
    depth = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    top = int(sys.argv[3]) if len(sys.argv) > 3 else 25
    peak, ps = peak_set(sys.argv[1])
    g = collections.defaultdict(lambda: [0, 0])
    for n, bt in ps.values():
        k = frames(bt, depth); g[k][0] += n; g[k][1] += 1
    print("PEAK live (logged >= threshold): %.1f MB in %d blocks" % (peak / 2**20, len(ps)))
    for k, (n, c) in sorted(g.items(), key=lambda x: -x[1][0])[:top]:
        print("%8.1f MB %4d  %s" % (n / 2**20, c, k))
