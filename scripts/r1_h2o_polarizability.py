#!/usr/bin/env python3
# scripts/r1_h2o_polarizability.py -- the ORACLE for doc/LinearResponsePlan.md stage R1 (molecular CPHF/CPKS).
#
# Static dipole polarisability of H2O in OUR geometry (IntegrationTests/M_Calculation.C MakeWater, bohr) and
# OUR basis (dzvp, CARTESIAN -- the facade's default), two independent routes that must agree:
#   * analytic CPHF / CPKS (pyscf.scf.cphf, the Z-vector engine), and
#   * finite field on the SCF (central difference, F = 1e-4 au).
#
# Run with the PySCF venv:   ~/Code/pyscf-env/bin/python scripts/r1_h2o_polarizability.py
#
# BANKED 2026-09-28 (PySCF 2.14.0):
#   HF  E = -76.02290319593      (== our M_Calculation.WaterEnergy anchor -76.022903: the basis MATCHES)
#   HF  alpha_xx, yy, zz = 3.19770284  7.11975051  5.54545175   (CPHF; finite field agrees to ~1e-4, its step noise)
#   LDA E = -75.87729931 (PySCF lda,vwn5, grid level 5)
#   LDA alpha_xx, yy, zz = 3.44581010  7.33075620  5.91599197   (CPKS; finite field agrees to ~1e-4, its step noise)
# ⚠ The LDA number is a LOOSE oracle for us: our molecular LDA (fitted Coulomb + fitted XC on a coarse
#   default mesh) sits at -75.93246, 55 mHa from PySCF's, so expect agreement at the percent level only.
#   The TIGHT LDA gate is internal: the FD kernel (src/Response/tests) through the same solver.
#   The HF number is tight: same basis, same energy to 1e-9.
#
# THE FACTOR-2 TRAP (cost one wrong number here): for a closed shell the first-order density built from the
# ov amplitudes carries the OCCUPANCY 2 -- dD = C_v (2U) C_o^T + h.c. -- or the kernel is half-counted.
from pyscf import gto, scf, dft
from pyscf.scf import cphf
import numpy as np

mol=gto.M(atom=[('O',(0,0,0)),('H',(0,1.431,1.107)),('H',(0,-1.431,1.107))],
          unit='Bohr', basis='dzvp', cart=True, verbose=0)

def alpha_cphf(mf):
    e, C, occ = mf.mo_energy, mf.mo_coeff, mf.mo_occ
    o = occ>0; Co, Cv = C[:,o], C[:,~o]
    with mol.with_common_orig((0,0,0)):
        r = mol.intor('int1e_r', comp=3)
    h1 = np.einsum('xpq,pa,qi->xai', r, Cv, Co)
    vind = mf.gen_response(singlet=None, hermi=1)
    def fvind(x):
        x = x.reshape(-1, Cv.shape[1], Co.shape[1])
        dm = np.einsum('xai,pa,qi->xpq', 2*x, Cv, Co)      # occupancy 2 (the trap above)
        dm = dm + dm.transpose(0,2,1)
        return np.einsum('xpq,pa,qi->xai', vind(dm), Cv, Co).ravel()
    U = cphf.solve(fvind, e, occ, h1, None, max_cycle=100, tol=1e-12)[0]
    return -4*np.einsum('xai,yai->xy', h1, U)

def alpha_ff(make, F=1e-4):
    with mol.with_common_orig((0,0,0)):
        r = mol.intor('int1e_r', comp=3)
    a = np.zeros((3,3))
    for j in range(3):
        mu = []
        for s in (+1,-1):
            m = make(); h0 = m.get_hcore()
            m.get_hcore = lambda *args, s=s, j=j, h0=h0: h0 + s*F*r[j]   # H = h + F.r (electron charge -1)
            m.conv_tol = 1e-12; m.kernel()
            mu.append(-np.einsum('xij,ji->x', r, m.make_rdm1()))
        a[:,j] = (mu[0]-mu[1])/(2*F)                            # mu(F) = mu0 + alpha F
    return a

def rhf():
    m = scf.RHF(mol); return m
def lda():
    m = dft.RKS(mol); m.xc = 'lda,vwn'; m.grids.level = 5; return m

np.set_printoptions(precision=8, suppress=True)
for name, make in (('HF', rhf), ('LDA', lda)):
    mf = make(); mf.conv_tol = 1e-12; E = mf.kernel()
    print(f"{name}  E = {E:.11f}")
    print(f"{name}  alpha CPHF/CPKS\n", alpha_cphf(mf))
    print(f"{name}  alpha finite field\n", alpha_ff(make))
