#!/usr/bin/env python3
# The independent route for the MOLECULAR DFT anchors (pin 10): PySCF with EXACT J and XC (no fit) on the SAME basis
# file our facade reads (Cartesian Gaussians), on a 200/590 grid.  Written 2026-10-31 (lane B) to re-judge the anchors
# moved by the molecular default XC mesh (doc/Records/OpenWork_History5.md, "default XC mesh").
#
# Run it with the PySCF venv:
#     ~/Code/pyscf-env/bin/python scripts/pyscf_molecule_oracle.py
#
# Prints, per case: the exact energy and the A1_coul Coulomb-fitted one (PySCF's density fitting with our own
# fit-basis file).  HF rows check the import: they must reproduce our HF anchors to 1e-8.
import os, re
from pyscf import gto, scf, dft

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'BasisSetData')
L = {'S': 0, 'P': 1, 'D': 2, 'F': 3}

def g94(name, el):
    """One element's block of a Gaussian94 .bsd file as a PySCF basis (contracted shells; SP split into S and P)."""
    txt = open(os.path.join(DATA, name)).read()
    m = re.search(r'\n %s +0\n(.*?)\n \*\*\*\*' % el.upper(), txt, re.S)
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

def mol(atoms, basis, spin=0, pseudo=None, cart=True):
    els = sorted({a[0] for a in atoms})
    return gto.M(atom=atoms, unit='Bohr', basis={e: g94(basis, e) for e in els}, spin=spin, cart=cart,
                 pseudo=pseudo, verbose=0)

def ks(m, xc, fit=False):
    mf = (dft.UKS(m) if m.spin else dft.RKS(m))
    if fit: mf = mf.density_fit(auxbasis={e: g94('A1_coul.bsd', e) for e in {a[0] for a in m._atom}})
    mf.xc = xc; mf.grids.atom_grid = (200, 590); mf.conv_tol = 1e-11
    return mf.kernel()

xa = lambda a: '%.10f*slater,' % (a * 1.5)   # our Xalpha(alpha) = alpha/(2/3) x Dirac exchange

WATER = [['O', (0, 0, 0)], ['H', (0, 1.431, 1.107)], ['H', (0, -1.431, 1.107)]]   # molecules.json H2O
N2    = [['N', (-1.03, 0, 0)], ['N', (1.04, 0, 0)]]                                 # M_DFT MakeN2
O2    = [['O', (0, 0, 1.141)], ['O', (0, 0, -1.141)]]                               # M_DFT MakeO2
F2    = [['F', (0, 0, -1.3343)], ['F', (0, 0, 1.3343)]]                             # ~/Code/cp2k-runs/F2

w, o2 = mol(WATER, 'dzvp.bsd'), mol(O2, 'dzvp.bsd', spin=2)
print('H2O  HF            %.8f   (ours M_Calculation.WaterEnergy -76.022903)' % scf.RHF(w).run(conv_tol=1e-11).e_tot)
print('O2   UHF triplet   %.8f   (ours M_DFT.OxygenTripletHF -149.6295470276)' % scf.UHF(o2).run(conv_tol=1e-11).e_tot)
for name, m, xc in [('H2O  LDA (vwn5)   ', w, 'lda,vwn'), ('H2O  Xa 0.74      ', w, xa(0.74)),
                    ('H2O  Xa 0.74 sph  ', mol(WATER, 'dzvp.bsd', cart=False), xa(0.74)),
                    ('N2   Xa 0.75197   ', mol(N2, 'dzvp.bsd'), xa(0.75197)),
                    ('O2   LDA triplet  ', o2, 'lda,vwn'), ('O2   Xa 0.7 trip  ', o2, xa(0.7)),
                    ('F2   LDA GTH q7   ', mol(F2, 'valence_lowq_vb.bsd', pseudo={'F': 'gth-pade-q7'}), 'lda,vwn')]:
    print('%s exact %.8f   A1_coul J-fit %.8f' % (name, ks(m, xc), ks(m, xc, fit=True)))
