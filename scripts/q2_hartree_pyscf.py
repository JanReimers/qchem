#!/usr/bin/env python3
# Q2 (lane S 2026-11-25, design note lane-s-design-ij-delta.md): PySCF's EXACT J as the reference energy and the COST REFERENCE for the
# molecular Hartree-route table.  LDA (Slater + VWN5, our "LDA"), the SAME basis files our facade reads (Cartesian), direct SCF:
#   - density_fit OFF (exact 4-centre J through the direct-SCF integral screen), direct_scf_tol STATED and printed;
#   - incremental Fock ON (PySCF's default for direct SCF: J of ΔD, rebuilt in full every `rebuild_nsteps`) -- said so in the output;
#   - XC on a converged grid (200, 590), so the energy is the reference our arms are judged against (|E - ref|).
# Timing split: SETUP = mol build + grid build + 1e integrals (everything before the first Fock); PER-ITERATION J = the wall time of each
# get_j call inside the SCF (wrapped), and the iteration count.  Peak RSS from getrusage.  Threads: OMP_NUM_THREADS as the caller sets it.
#
#     OMP_NUM_THREADS=4 ~/Code/pyscf-env/bin/python scripts/q2_hartree_pyscf.py H2O [CH4 SF6 TeF6] [--pp] [--tol 1e-13]
import os, re, sys, time, resource, json
from pyscf import gto, dft, lib

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'BasisSetData')
L = {'S': 0, 'P': 1, 'D': 2, 'F': 3}
def g94(name, el):
    """One element's block of a Gaussian94 .bsd file as a PySCF basis (as pyscf_molecule_oracle.py, which runs its cases on import)."""
    txt = open(os.path.join(DATA, name)).read()
    m = re.search(r'\n %s +0\n(.*?)\n \*\*\*\*' % el.upper(), txt, re.S)
    if not m: raise SystemExit(f'{name} has no block for {el}')
    lines, out, i = m.group(1).split('\n'), [], 0
    while i < len(lines):
        t = lines[i].split()
        if not t: i += 1; continue
        n = int(t[1])
        prims = [list(map(float, lines[i+1+k].replace('D', 'E').split())) for k in range(n)]
        if t[0] == 'SP':
            out.append([0] + [[p[0], p[1]] for p in prims]); out.append([1] + [[p[0], p[2]] for p in prims])
        else:
            out.append([L[t[0]]] + [[p[0], p[1]] for p in prims])
        i += n + 1
    return out

MOLS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src', 'Structure', 'Data', 'molecules.json')

def atoms(name):
    e = json.load(open(MOLS))[name]
    return [[a['el'], tuple(a['xyz'])] for a in e['atoms']]

def run(name, basis, pseudo, tol):
    t0 = time.perf_counter()
    at = atoms(name)
    els = sorted({a[0] for a in at})
    m = gto.M(atom=at, unit='Bohr', basis={e: g94(basis, e) for e in els}, cart=True,
              pseudo=({e: pseudo for e in els} if pseudo else None), verbose=0, max_memory=8000)
    mf = dft.RKS(m); mf.xc = 'lda,vwn'; mf.grids.atom_grid = (200, 590)
    mf.direct_scf = True; mf.direct_scf_tol = tol; mf.conv_tol = 1e-10; mf.max_cycle = 200
    mf.grids.build(); h1 = mf.get_hcore(); s1 = mf.get_ovlp()
    mf.get_hcore = lambda *a, **k: h1; mf.get_ovlp = lambda *a, **k: s1
    tSetup = time.perf_counter() - t0
    jt = []
    getj = mf.get_j
    def timed_get_j(*a, **k):
        t = time.perf_counter(); r = getj(*a, **k); jt.append(time.perf_counter() - t); return r
    mf.get_j = timed_get_j
    # RKS for a pure functional builds J through get_jk/get_j inside get_veff; time the whole veff too (J + XC on the 200/590 grid)
    vt = []
    getv = mf.get_veff
    def timed_veff(*a, **k):
        t = time.perf_counter(); r = getv(*a, **k); vt.append(time.perf_counter() - t); return r
    mf.get_veff = timed_veff
    t1 = time.perf_counter(); E = mf.kernel(); tScf = time.perf_counter() - t1
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    J = mf.get_j(dm=mf.make_rdm1()); EJ = 0.5 * float((mf.make_rdm1() * J).sum())
    print(f'{name:5s} {basis:22s} pp={pseudo or "AE":12s} E={E:.10f} conv={mf.converged} iters={mf.cycles} EJ={EJ:.10f}  '
          f'setup={tSetup:.2f}s scf={tScf:.2f}s J/call={sum(jt)/max(1,len(jt)):.3f}s x{len(jt)} veff/call={sum(vt)/max(1,len(vt)):.3f}s x{len(vt)} '
          f'peakRSS={rss:.0f}MB direct_scf_tol={tol:g} incremental_fock=on(rebuild_nsteps={getattr(mf,"rebuild_nsteps","?")}) '
          f'threads={lib.num_threads()}', flush=True)

if __name__ == '__main__':
    args = sys.argv[1:]
    tol = 1e-13
    if '--tol' in args: i = args.index('--tol'); tol = float(args[i+1]); del args[i:i+2]
    pp = '--pp' in args; args = [a for a in args if a != '--pp']
    for name in args or ['H2O']:
        if pp: run(name, 'valence_lowq_vb.bsd', 'gth-pade', tol)
        else:  run(name, 'dzvp.bsd', None, tol)
