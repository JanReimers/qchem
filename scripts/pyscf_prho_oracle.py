#!/usr/bin/env python3
# The independent route for the A1 P-rho XC cell (VxcFit::Auxiliary at XCρ::Fit, V-XC-ONE-TERM step 3, lane S 2026-11-10; pin 10).
# At PySCF's converged density (exact XC, the A1_coul J-fit -- the route scripts/pyscf_molecule_oracle.py prints), project rho onto
# the A1_exch Vxc fit basis in the OVERLAP metric (c = S^-1 <a|rho>, the analytic 3-centre overlap) and compare E_xc[P rho] with
# E_xc[rho] (LDA = Slater + VWN5) on PySCF's own grid.  That first-order shift is what our self-consistent P-rho energy should sit
# near; self-consistency lowers it further (variational), by second order.
# Measured 2026-11-10:  water dzvp  dE_xc -18.298 mHa (ours, SCF: -19.57);  F2 GTH q7  -4.782 mHa (ours, SCF: -4.91).
#
#     ~/Code/pyscf-env/bin/python -I scripts/pyscf_prho_oracle.py BasisSetData
import os, re, sys
import numpy as np
from pyscf import gto, dft
DATA = sys.argv[1]
L = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4}
def g94(name, el):
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
def run(atoms, basis, pseudo=None):
    els = sorted({a[0] for a in atoms})
    m = gto.M(atom=atoms, unit='Bohr', basis={e: g94(basis, e) for e in els}, cart=True, pseudo=pseudo, verbose=0)
    aux = gto.M(atom=atoms, unit='Bohr', basis={e: g94('A1_exch.bsd', e) for e in els}, cart=True, verbose=0)
    mf = dft.RKS(m).density_fit(auxbasis={e: g94('A1_coul.bsd', e) for e in els})
    mf.xc = 'lda,vwn'; mf.grids.atom_grid = (200, 590); mf.conv_tol = 1e-11
    E = mf.kernel(); D = mf.make_rdm1()
    # p_a = sum_ij D_ij <ij|a> (3-centre overlap), c = S^-1 p
    pm = gto.mole.conc_mol(m, aux)
    nao = m.nao
    int3 = pm.intor('int3c1e_cart', shls_slice=(0, m.nbas, 0, m.nbas, m.nbas, pm.nbas))
    p = np.einsum('ij,ija->a', D, int3)
    c = np.linalg.solve(aux.intor('int1e_ovlp_cart'), p)
    g = mf.grids; ao = m.eval_gto('GTOval_cart', g.coords); chi = aux.eval_gto('GTOval_cart', g.coords)
    rho = np.einsum('gi,ij,gj->g', ao, D, ao); prho = chi @ c
    exc = lambda r: np.dot(g.weights, dft.libxc.eval_xc('lda,vwn', r)[0] * r)
    print('E_scf(exact XC)=%.8f  N[rho]=%.6f N[Prho]=%.6f  min Prho=%.3e  E_xc[rho]=%.8f E_xc[Prho]=%.8f  dE_xc=%+.6f Ha'
          % (E, np.dot(g.weights, rho), np.dot(g.weights, prho), prho.min(), exc(rho), exc(np.maximum(prho, 0)), exc(np.maximum(prho, 0)) - exc(rho)))
print('water'); run([['O', (0, 0, 0)], ['H', (0, 1.431, 1.107)], ['H', (0, -1.431, 1.107)]], 'dzvp.bsd')
print('F2 GTH'); run([['F', (0, 0, -1.3343)], ['F', (0, 0, 1.3343)]], 'valence_lowq_vb.bsd', pseudo={'F': 'gth-pade-q7'})
