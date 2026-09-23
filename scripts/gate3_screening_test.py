#!/usr/bin/env python3
# scripts/gate3_screening_test.py -- doc/HubbardUPlan.md GATE 3, run ON PAPER before any kernel is written.
#
# THE QUESTION: does ONE screening length move EVERY measured manifold from its BARE on-site F0 onto its
# INDEPENDENT oracle?  If yes, the residual between our ACBN0 and linear response is plausibly dielectric
# screening and route (b) is worth building.  If the manifolds need systematically different lengths, the
# residual is NOT bulk screening and (b) is refuted.
#
# Run it with the PySCF venv (installed 2026-09-23, PySCF 2.14.0 on Python 3.14):
#     ~/Code/pyscf-env/bin/python scripts/gate3_screening_test.py
#
# WHY PySCF: it carries range-separated (omega) four-index integrals, so the screened on-site average can be
# evaluated WITHOUT building anything in our tree first.  The contractions below are OUR OWN -- read off the
# `[+U radial]` banner of a gpwprobe run, so this asks about our manifold, not a model one.
#
# ⚠ Targets are heterogeneous on purpose and must stay labelled: hp.x is a MATCHED-PP oracle (gth2upf ran our
# own PP and projector); the O 2p figure is a literature cRPA BOUND, not a value.  Turning that bound into a
# value is what ABINIT `ucrpa` is for -- and it is the single biggest improvement available to this test.
# GATE 3's REFUTATION TEST, on paper: does ONE screening length move BOTH manifolds onto their
# independent oracles?  Contractions are OUR OWN, read off the gpwprobe [+U radial] banner.
import numpy as np
from pyscf import gto
HA = 27.211386245988
MAN = {
 'Ni 3d (l=2)': (2, [0.180,0.39485815,0.86618308,1.90010803,4.16818406,9.14356349,20.05783621,44.0],
                    [-0.114,-0.154,-0.292,-0.308,-0.387,-0.047,0.009,-0.002], 5.27, 'hp.x (matched PP)'),
 'O  2p (l=1)': (1, [0.465,1.200,3.098,8.000],
                    [0.705,0.092,0.335,0.045],                                  4.0,  'cRPA bound >~4'),
}
def F0(l, exps, cs, omega, elem):
    mol = gto.M(atom=f'{elem} 0 0 0', verbose=0, spin=None,
                basis={elem: [[l] + [[e,c] for e,c in zip(exps,cs)]]})
    n = mol.nao
    with mol.with_range_coulomb(omega):
        eri = mol.intor('int2e').reshape(n,n,n,n)
    return np.mean([eri[m,m,p,p] for m in range(n) for p in range(n)])
print(f"{'manifold':14s} {'bare F0':>9s} {'target':>8s} {'needed ratio':>13s} {'omega needed':>13s}")
for name,(l,e,c,target,src) in MAN.items():
    elem = 'Ni' if 'Ni' in name else 'O'
    bare = F0(l,e,c,0.0,elem)*HA
    want = target/bare
    lo,hi = 0.01, 5.0                      # bisect on omega for the erfc-screened average
    for _ in range(60):
        mid=0.5*(lo+hi)
        r = F0(l,e,c,-mid,elem)*HA/bare
        if r > want: lo=mid
        else:        hi=mid
    print(f"{name:14s} {bare:8.2f}eV {target:7.2f}eV {want:13.3f} {0.5*(lo+hi):10.3f} a.u.  "
          f"(1/omega = {1/(0.5*(lo+hi)):.2f} bohr)   [{src}]")
