#!/usr/bin/env python3
# scripts/a1_bare_f0_check.py -- doc/HubbardUPlan.md A1 (RESOLVED 2026-09-23): the ACBN0 paper's eq-10
# "bare Ubar" and PySCF's plain shell-averaged F0 are DIFFERENT quantities by construction, not two
# measurements of the same integral.  Eq 10's numerator is an unrestricted double sum over the density
# matrix; its pair-count denominator (eq 10c) excludes same-orbital/same-spin self-pairs.  That asymmetry
# inflates Ubar above the plain average by a factor set by the spin-resolved occupation alone -- nothing to
# do with the integrals themselves.
#
# Run with the PySCF venv: ~/Code/pyscf-env/bin/python scripts/a1_bare_f0_check.py
#
# Verified against a live run: `gpwprobe nio` (NIO_ACBN0=1 GPW_SPHERICAL=1 NIO_U_RADIAL=ortho, U=0)
# reported bare N_up/N_dn = 4.9669/4.9667 (site 0) and 2.9597/2.8819 (site 1), and bare Ubar = 27.018 /
# 27.131 eV.  This script's uniform-within-spin approximation predicts 26.99 eV for both -- 0.1-0.5 % off,
# the residual being the real (crystal-field-split, non-uniform) per-orbital occupation the approximation
# ignores.
import numpy as np
from pyscf import gto
HA = 27.211386245988

l, exps, cs = 2, [0.180,0.39485815,0.86618308,1.90010803,4.16818406,9.14356349,20.05783621,44.0], \
                  [-0.114,-0.154,-0.292,-0.308,-0.387,-0.047,0.009,-0.002]     # Ni 3d, gpwprobe's [+U radial] banner
mol = gto.M(atom='Ni 0 0 0', verbose=0, spin=None, basis={'Ni': [[l] + [[e,c] for e,c in zip(exps,cs)]]})
n = mol.nao
eri = mol.intor('int2e').reshape(n,n,n,n)

F0_plain = np.mean([eri[m,m,p,p] for m in range(n) for p in range(n)]) * HA
print(f"PySCF plain shell average F0 (unweighted, matches src/BasisSet/Gaussian/tests/M_BareCoulomb.C's "
      f"convention): {F0_plain:.4f} eV")

def acbn0_Ubar(Na_per, Nb_per, eri):
    """eq 10a/10c: unrestricted numerator, pair-count denominator excluding same-spin self-pairs.
    Assumes P diagonal in the manifold's own (crystal-field-adapted) basis, occupation uniform within
    each spin -- exact only for a spherically symmetric fill; the real thing differs at the percent level."""
    m = len(Na_per)
    Nt = np.array(Na_per) + np.array(Nb_per)
    numU = sum(Nt[a]*Nt[c]*eri[a,a,c,c] for a in range(m) for c in range(m))
    Na, Nb = sum(Na_per), sum(Nb_per)
    Na2, Nb2 = sum(x*x for x in Na_per), sum(x*x for x in Nb_per)
    denU = (Na+Nb)**2 - Na2 - Nb2
    return numU/denU

print()
print(f"{'Na(up)':>7} {'Nb(dn)':>7} {'Ubar(eV)':>10} {'ratio to F0_plain':>18}")
for Na, Nb in [(5,0), (5,3), (4.980,3.368), (5,5)]:
    U = acbn0_Ubar([Na/5]*5, [Nb/5]*5, eri) * HA
    print(f"{Na:7.3f} {Nb:7.3f} {U:10.3f} {U/F0_plain:18.4f}")

print()
print("Live gpwprobe nio run (U=0, AFM-II collapsed -- trap 3, a diagnostic only):")
for label, Na, Nb, actual in [("site 0", 4.9669, 4.9667, 27.0180), ("site 1", 2.9597, 2.8819, 27.1308)]:
    predicted = acbn0_Ubar([Na/5]*5, [Nb/5]*5, eri) * HA
    print(f"  {label}: Na={Na} Nb={Nb} -> predicted {predicted:.3f} eV vs actual banner {actual:.3f} eV"
          f"  ({100*abs(predicted-actual)/actual:.2f} % apart)")
