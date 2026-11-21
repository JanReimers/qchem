#!/usr/bin/env python3
"""gga_pyscf_oracle -- the PySCF oracle for the molecular GGA gates (IntegrationTests/M_GGA.C, GGA stage 2, lane S 2026-11-21).

THE SAME MODEL as our MoleculeCalculation on the delta route:
  * basis DZVP (DGauss), CARTESIAN d shells (ours are Cartesian by default), read from OUR BasisSetData/dzvp.bsd;
  * the Coulomb FIT: DGauss A1 cfit, read from OUR BasisSetData/A1_coul.bsd, and DUNLAP'S CHARGE-CONSTRAINED fit in the Coulomb
    metric, as ours (Fitting::ConstrainedFF): c = G^-1 (p - lam n) with n.c = N, the ROBUST energy (rho|rho~) - 1/2 (rho~|rho~), and
    J = (mu nu|P) c + lam S, its exact derivative.  PySCF's own density fitting is UNCONSTRAINED: it maximizes the robust energy, so
    it sits ABOVE ours by 2-46 uHa on these systems (He the most: its A1 fit is four s functions) -- a different fit, not an error;
  * libxc LDA = LDA_X + LDA_C_VWN (VWN5) and PBE = GGA_X_PBE + GGA_C_PBE, unrestricted for an open shell (ours: spin-native);
  * the XC grid converged: two fine PySCF grid levels are printed, and ours is converged on its own Becke ladder (M_GGA.C).
Geometries in BOHR, copied from IntegrationTests (H2O from src/Structure/Data/molecules.json; N2 and O2 from M_DFT.C).

Run it with the PySCF venv:
    ~/Code/pyscf-env/bin/python scripts/gga_pyscf_oracle.py
"""
import os
import sys
import numpy as np
from pyscf import gto, dft, lib
from pyscf import df as pdf

# OUR OWN basis files, parsed by PySCF: the one way to be sure both codes read the same exponents and contractions.
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "BasisSetData")


L_OF = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4}


def read_g94(path, el):
    """One element's block of a Gaussian94 basis file, as a PySCF basis [[l, [exp, coef], ...], ...].  An SP shell becomes an s and a p
    shell on the same exponents.  Our files carry a free-text header PySCF's own parser rejects, hence this reader."""
    lines = open(path).read().splitlines()
    i = next(k for k, ln in enumerate(lines) if ln.split()[:2] == [el.upper(), "0"])
    out, k = [], i + 1
    while lines[k].strip() != "****":
        tok = lines[k].split()
        kind, n = tok[0].upper(), int(tok[1])
        rows = [[float(x.replace("D", "E")) for x in lines[k + 1 + j].split()] for j in range(n)]
        if kind == "SP":
            out.append([0] + [[r[0], r[1]] for r in rows])
            out.append([1] + [[r[0], r[2]] for r in rows])
        else:
            out.append([L_OF[kind]] + [[r[0], r[1]] for r in rows])
        k += 1 + n
    return out


def basis_of(atoms, fname):
    """{element: PySCF basis} from our Gaussian94 file, for every element in the atom string."""
    els = sorted({a.split()[0] for a in atoms.split(";")})
    return {el: read_g94(os.path.join(DATA, fname), el) for el in els}

MOLECULES = {
    "He":  ("He 0 0 0", 0),
    "Ne":  ("Ne 0 0 0", 0),
    "H2O": ("O 0 0 0; H 0 1.431 1.107; H 0 -1.431 1.107", 0),
    "N2":  ("N -1.03 0 0; N 1.04 0 0", 0),
    "O2":  ("O 0 0 1.141; O 0 0 -1.141", 2),     # triplet
    "N":   ("N 0 0 0", 3),                       # quartet
}
FUNCTIONALS = {"LDA": "lda_x,lda_c_vwn", "PBE": "gga_x_pbe,gga_c_pbe"}
# The other two stage-1 GGAs, on water only (the same libxc ids as Model::PBEsol / Model::BLYP: 116+133, 106+131).
EXTRA = {"PBEsol": "gga_x_pbe_sol,gga_c_pbe_sol", "BLYP": "gga_x_b88,gga_c_lyp"}
GRIDS = [(150, 590), (250, 974)]   # (radial, Lebedev points) per atom: two rungs, to show PySCF's own grid residual


def dunlap_veff(mf, auxmol):
    """mf.get_veff with Dunlap's charge-constrained Coulomb fit (see the module note), XC from PySCF's own numint on mf.grids."""
    mol = mf.mol
    int3c = pdf.incore.aux_e2(mol, auxmol, intor="int3c2e", aosym="s1")          # (mu nu|P)
    G = auxmol.intor("int2c2e")                                                   # (P|Q)
    S = mol.intor("int1e_ovlp")
    mf.grids.build()
    chi = auxmol.eval_gto("GTOval_cart" if auxmol.cart else "GTOval_sph", mf.grids.coords)
    n = chi.T @ mf.grids.weights                                                  # n_P = int chi_P (quadrature; fine grid)
    Ginv_n = np.linalg.solve(G, n)
    unrestricted = isinstance(mf, dft.uks.UKS)

    def get_veff(mol_=None, dm=None, dm_last=0, vhf_last=0, hermi=1):
        if dm is None: dm = mf.make_rdm1()
        dmt = dm[0] + dm[1] if unrestricted else dm
        p = np.einsum("ijP,ij->P", int3c, dmt)
        Ginv_p = np.linalg.solve(G, p)
        N = np.einsum("ij,ij", dmt, S)
        lam = (n @ Ginv_p - N) / (n @ Ginv_n)
        c = Ginv_p - lam * Ginv_n
        vj = np.einsum("ijP,P->ij", int3c, c) + lam * S
        ecoul = p @ c - 0.5 * c @ G @ c
        if unrestricted:
            _, exc, vxc = mf._numint.nr_uks(mol, mf.grids, mf.xc, dm)
            veff = vxc + vj[None, :, :]
        else:
            _, exc, vxc = mf._numint.nr_rks(mol, mf.grids, mf.xc, dm)
            veff = vxc + vj
        return lib.tag_array(veff, ecoul=ecoul, exc=exc, vj=None, vk=None)
    return get_veff


def energy(atoms, spin, xc, grid):
    mol = gto.M(atom=atoms, unit="Bohr", basis=basis_of(atoms, "dzvp.bsd"), cart=True, spin=spin, verbose=0)
    auxmol = gto.M(atom=atoms, unit="Bohr", basis=basis_of(atoms, "A1_coul.bsd"), cart=True, spin=spin, verbose=0)
    mf = dft.UKS(mol) if spin else dft.RKS(mol)
    mf.xc = xc
    mf.grids.atom_grid = grid
    mf.grids.prune = None
    mf.get_veff = dunlap_veff(mf, auxmol)
    mf.conv_tol = 1e-11
    mf.max_cycle = 200
    e = mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"PySCF did not converge: {atoms} {xc} {grid}")
    return e


def main():
    names = sys.argv[1:] or list(MOLECULES)
    for name in names:
        atoms, spin = MOLECULES[name]
        funcs = dict(FUNCTIONALS, **(EXTRA if name == "H2O" else {}))
        for fname, xc in funcs.items():
            es = [energy(atoms, spin, xc, g) for g in GRIDS]
            print(f"{name:4s} {fname}: " + "  ".join(f"{g[0]}/{g[1]} {e:.10f}" for g, e in zip(GRIDS, es))
                  + f"   (grid residual {es[1]-es[0]:+.2e})", flush=True)


if __name__ == "__main__":
    main()
