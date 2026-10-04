#!/usr/bin/env python3
"""log2deck -- reconstruct a DRAFT input deck from the banner lines of an old qchem log (D-ENV-RERUN, 2026-10-04).

The campaign logs in ~/Code/qchem6-runs were produced by `gpwprobe` with MNO_*/NIO_*/GATE1_* environment knobs that were not saved
(only two `.cmd` files exist).  What the run PRINTED about itself -- system line, grids, symmetry, +U manifolds, the [ortho] pivot line, the
scf mixer line per stage, the vet-stage trim, restart/save paths, the chi actions -- is enough for a deck that `rundeck --check` can vet.
It is a DRAFT: every value the banner does not carry is listed in the deck's `_doc` and must be judged by the person who ran it
(the convergence threshold above all; 1e-6 is assumed).  The basis span is chosen by trying the candidates through `rundeck --check` and
keeping the one whose function count equals the log's own `nFunctions`.

    scripts/log2deck.py <log> [--out decks/campaigns/<Material>/<name>.json] [--rundeck build/Release/CLIapps/rundeck]
"""
import json, os, re, subprocess, sys

def parse(path):
    L = open(path, errors='replace').read().split('\n')
    d = {"_notes": []}
    name = os.path.basename(path)
    for l in L:
        m = re.match(r"\[(.+?) run\] system: (\d+) atoms, (\d+) valence e, multiplicity (\d+)", l)
        if m: d['label'], d['mult'] = m.group(1), int(m.group(4)); break
    if 'label' not in d: return None
    lab = d['label']
    first = lab.split()[0]
    d['structure'] = (first + "_AFM2") if "AFM-II" in lab else first
    d['imposed'] = any("symmetry: IMPOSED" in l for l in L)
    d['grey'] = any("grey" in l and "symmetry: IMPOSED" in l for l in L)
    for l in L:
        m = re.search(r"xcMesh=(\w+) \(nR=(\d+) L=(\d+)\)", l)
        if m: d['xc'] = (m.group(1), int(m.group(2)), int(m.group(3))); break
    for l in L:
        if "+U:" in l and "manifolds:" in l:
            man = []
            for mm in re.finditer(r"\(site (\d+), l=(\d+), U=([\d.eE+-]+) eV, radial: ([^)]*(?:\([^)]*\))?[^)]*)\)", l):
                site, ll, U, rad = int(mm.group(1)), int(mm.group(2)), float(mm.group(3)), mm.group(4)
                e = {"site": site, "l": ll, "U_eV": U}
                if "ONE contracted" in rad: e["atomicRadial"] = True
                if "ORTHO-atomic" in rad: e["orthoAtomic"] = True
                man.append(e)
            d['hubbard'] = man; break
    for l in L:
        m = re.match(r"\[ortho\] (CholeskyPivoted|Cholesky): tol=([\d.eE+-]+)", l)
        if m: d['ortho'], d['orthoTol'] = m.group(1), float(m.group(2)); break
    stages = []
    for l in L:
        m = re.search(r"scf\] mixer: Kerker\(G0=([\d.]+)\) \+ Pulay history\(depth (\d+), start (\d+)\) alpha=([\d.]+);.*accel: (\w+);\s+kT=([\d.eE+-]+) MOM=(on|off) NMaxIter=(\d+)", l)
        if m:
            s = dict(G0=float(m.group(1)), depth=int(m.group(2)), start=int(m.group(3)), alpha=float(m.group(4)), acc=m.group(5),
                     kT=float(m.group(6)), mom=m.group(7) == "on", nmax=int(m.group(8)))
            if not stages or stages[-1] != s: stages.append(s)
    d['stages'] = stages
    for l in L:
        m = re.match(r"\[basis trim\] vet stage pass \d+: .*orthoTol=([\d.eE+-]+)", l)
        if m: d['vet'] = True; d['vetTol'] = float(m.group(1)); break
    for l in L:
        m = re.match(r"\[basis trim\] STATED", l)
        if m: d['_notes'].append("a STATED shell trim was used (basis.trim) -- the banner line does not list Z:l:alpha in a parseable form; add it by hand")
    for l in L:
        m = re.search(r"EXACT RESUME from (\S+)", l)
        if m: d['restart'] = m.group(1); break
    for l in L:
        m = re.search(r"state saved to (\S+)", l)
        if m: d['save'] = m.group(1); break
    for l in L:
        m = re.search(r"nFunctions\s+(\d+)", l)
        if m: d['nfun'] = int(m.group(1)); break
    post = []
    for l in L:
        m = re.search(r"\[chi0\] independent-particle channel response, q-mesh (\d+)x", l)
        if m and not any(p.get('independentResponse') for p in post): post.append({"independentResponse": {"nq": int(m.group(1))}})
        if "self-consistent chi (q=0) from the last iterate" in l and not any('hubbardLinearResponse' in p for p in post):
            post.append({"hubbardLinearResponse": {}})
    d['post'] = post
    km = re.search(r"k(\d)(\d)(\d)", name)
    nk = None
    for l in L:
        m = re.match(r"\[IBZ\] (\d+) k-points", l)
        if m: nk = int(m.group(1)); break
    n = round(nk ** (1 / 3)) if nk else 1
    if nk and n ** 3 == nk: d['kmesh'] = [n, n, n]
    elif km: d['kmesh'] = [int(km.group(i)) for i in (1, 2, 3)]
    else:
        d['kmesh'] = [1, 1, 1]
        nq = [list(p.values())[0].get('nq', 1) for p in post if 'independentResponse' in p]
        if nq and nq[0] > 1: d['kmesh'] = [nq[0]] * 3; d['_notes'].append("k-mesh taken from the chi0 q-mesh (the q-mesh must divide it)")
    if nk is None and not km: d['_notes'].append("k-mesh not in the banner ([IBZ] line absent): Gamma assumed unless a q-mesh implies more")
    if nk: d['_notes'].append("k-mesh from the log's [IBZ] line: %d k-points" % nk)
    return d

def deck(d, src, basis_data, spherical):
    sol = {"multiplicity": d['mult'], "seed": "IonicSAD", "imposeSymmetry": d['imposed']}
    if d['grey']: sol["greyImposition"] = True
    if 'ortho' in d: sol["ortho"], sol["orthoTol"] = d['ortho'], d['orthoTol']
    if d.get('xc'):
        k, nr, L = d['xc']
        if k == "Uniform": sol["xcMesh"] = {"cellKind": "Uniform"}
        elif (nr, L) != (40, 29): sol["xcMesh"] = {"cellKind": "Becke", "nRadial": nr, "angularDegree": L, **({"angular": "GaussLegendre"} if L < 29 else {})}
    if d.get('hubbard'): sol["hubbard"] = d['hubbard']
    if d['post']: sol["forceComplex"] = True
    st = d['stages'] or [dict(G0=1.0, depth=8, start=5, alpha=0.45, acc="Null", kT=5e-3, mom=False, nmax=200)]
    stages = []
    for i, s in enumerate(st):
        scf = {"NMaxIter": s['nmax'], "minDeltaRho": 1e-6, "minDeltaE": 1e30, "minDeltaFD": 1e30, "minVirial": 1e30, "minFD": 1e30,
               "startingRelaxRo": s['alpha'], "mergeTol": 1e-4, "pulayDepth": s['depth'], "pulayStart": s['start'], "kerkerG0": s['G0'],
               "useMOM": s['mom'], "smearingkT": s['kT']}
        if len(st) > 1 and i + 1 < len(st): scf["stopOnAccelExhausted"] = True
        stages.append({"accelerator": s['acc'], "scf": scf})
    basis = {"data": basis_data, "spherical": spherical}
    if d.get('vet'): basis["vet"] = True; sol["orthoTol"] = d['vetTol']
    out = {"_doc": ["DRAFT reconstructed from the banner of " + src + " (scripts/log2deck.py, D-ENV-RERUN 2026-10-04): the original `gpwprobe` environment command was not saved.",
                    "NOT in the banner, ASSUMED: convergence threshold minDeltaRho=1e-6 on the default measure, mergeTol 1e-4, the minDelta{E,FD}/minVirial/minFD gates off (1e30).",
                    "Judge these against the run's author before trusting a re-run's numbers.  `rundeck --check` has vetted that it parses, resolves and builds its lattice and basis."] + d['_notes'],
           "structure": d['structure'], "kmesh": d['kmesh'], "basis": basis, "solid": sol}
    if len(stages) == 1:
        out["scf"] = stages[0]["scf"]; sol["accelerator"] = stages[0]["accelerator"]
    else:
        out["schedule"] = stages
    if d['post']: out["postSCF"] = d['post']
    if d.get('save'): out["state"] = {"save": "auto"}
    if d.get('restart'): out.setdefault("state", {})["restartFrom"] = d['restart']
    return out

def nfun(rundeck, path):
    r = subprocess.run([rundeck, path, "--check"], capture_output=True, text=True)
    if r.returncode != 0: return None, (r.stderr.strip() or r.stdout.strip())[:300]
    m = re.search(r"(\d+) functions", r.stdout)
    return (int(m.group(1)) if m else None), r.stdout

if __name__ == "__main__":
    a = sys.argv[1:]
    log = a[0]
    out = a[a.index("--out") + 1] if "--out" in a else None
    rd = a[a.index("--rundeck") + 1] if "--rundeck" in a else "build/Release/CLIapps/rundeck"
    d = parse(log)
    if d is None: print("SKIP  %s: no qchem run banner" % log); sys.exit(0)
    cands = [("VALENCE_LOWQ_VA", True), ("VALENCE_LOWQ_SR", False), ("VALENCE_LOWQ_VA", False), ("VALENCE_LOWQ_SPH", True), ("VALENCE_LOWQ_VB", True)]
    if not d.get('hubbard'): cands = [cands[1], cands[0]] + cands[2:]
    tmp = (out or "/tmp/log2deck.json") + ".tmp"
    os.makedirs(os.path.dirname(os.path.abspath(tmp)), exist_ok=True)
    chosen, note, err = None, "", ""
    for bd, sph in cands:
        dk = deck(d, os.path.basename(log), bd, sph); json.dump(dk, open(tmp, "w"), indent=2)
        n, msg = nfun(rd, tmp)
        if n is None: err = msg; continue
        if d.get('vet') or 'nfun' not in d:
            chosen = (bd, sph, n); note = "function count not comparable (%s)" % ("vet-stage trim" if d.get('vet') else "no nFunctions in log"); break
        if n == d['nfun']: chosen = (bd, sph, n); note = "function count %d == the log's" % n; break
        err = "built %d functions, log says %d" % (n, d['nfun'])
    if not chosen: print("FAIL  %s: %s" % (log, err)); sys.exit(1)
    bd, sph, n = chosen
    dk = deck(d, os.path.basename(log), bd, sph); dk["_doc"].append("basis span chosen: %s%s; %s" % (bd, " spherical" if sph else "", note))
    json.dump(dk, open(out or tmp, "w"), indent=2); 
    if out and os.path.exists(tmp): os.remove(tmp)
    print("OK    %s -> %s  [%s%s, %s]" % (os.path.basename(log), out, bd, " sph" if sph else "", note))
