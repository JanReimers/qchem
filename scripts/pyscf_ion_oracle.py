#!/usr/bin/env python3
# pyscf_ion_oracle -- the SECOND-CODE oracle for gen-ae-seeds' ions (OpenWork row "All-electron ION seeds", controller ruling 2026-11-15).
# A second-code oracle for a spherically averaged, high-spin LSDA atom or ion in a STATED configuration: PySCF UKS (LDA = Slater + VWN5),
# a large uncontracted basis, with the occupations SET per (l, spin) block -- each l block filled from its lowest MO up with the
# configuration's per-spin counts, an open shell spread evenly over its 2l+1 orbitals (the spherical ensemble: the radial atom's state).
# Validated 2026-11-15 against NIST's LSD tables: Fe [Ar] 3d^6 4s^2 0.59 mHa and Fe+ [Ar] 3d^6 4s^1 0.62 mHa ABOVE NIST (the unc-ano-rcc
# basis gap), so it judges an ion NIST does not tabulate (Mn2+, Fe3+, ...) to well inside gen-ae-seeds' 1e-5 relative tolerance.
# Prints one JSON line {"E", "converged", "basis"}.  Run with the PySCF venv:
#   ~/Code/pyscf-env/bin/python scripts/pyscf_ion_oracle.py SYMBOL CHARGE "CONFIG" [basis]   e.g.  Fe 1 "[Ar] 3d^6 4s^1"
import sys, re, numpy as np
from pyscf import gto, dft, lib

NG = {'He':'1s2','Ne':'1s2 2s2 2p6','Ar':'1s2 2s2 2p6 3s2 3p6','Kr':'1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6',
      'Xe':'1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6','Rn':'1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6'}
L = 'spdf'
def shells(cfg):
    out = []
    m = re.match(r'\s*\[(\w+)\](.*)', cfg)
    if m: cfg = NG[m[1]] + ' ' + m[2]
    for t in cfg.replace('^','').split():
        mm = re.match(r'(\d)([spdf])(\d+)$', t); out.append((int(mm[1]), L.index(mm[2]), int(mm[3])))
    return out

def spin_counts(sh):
    """per l: (up, down) electron counts, Hund high spin within each open subshell (majority spin up)."""
    up = [0]*4; dn = [0]*4
    for n,l,k in sh:
        g = 2*l+1
        u = min(k, g); up[l] += u; dn[l] += k-u
    return up, dn

def run(sym, charge, cfg, basis='unc-ano-rcc'):
    sh = shells(cfg); up, dn = spin_counts(sh)
    nel = sum(k for _,_,k in sh)
    mol = gto.M(atom=f'{sym} 0 0 0', basis=basis, charge=charge, spin=sum(up)-sum(dn), verbose=0)
    assert mol.nelectron == nel, (mol.nelectron, nel)
    lab = np.array([L.index(re.search(r'\d+([spdfghi])', s.split()[2]).group(1)) if re.search(r'\d+([spdf])', s.split()[2]) else 9
                    for s in mol.ao_labels()])
    mf = dft.UKS(mol); mf.xc = 'lda,vwn'; mf.grids.level = 7; mf.conv_tol = 1e-10; mf.max_cycle = 300
    S = mol.intor('int1e_ovlp')
    def get_occ(mo_energy=None, mo_coeff=None):
        occ = []
        for s,(cnt) in enumerate((up, dn)):
            e, C = mo_energy[s], mo_coeff[s]
            w = np.array([[ (C[:,i] @ S[:, lab==l]) @ C[lab==l, i] for l in range(4)] for i in range(C.shape[1])])
            lmo = w.argmax(axis=1)
            o = np.zeros(len(e))
            for l in range(4):
                idx = [i for i in np.argsort(e) if lmo[i]==l]
                k = cnt[l]; g = 2*l+1; j = 0
                while k > 0:
                    blk = idx[j:j+g]; take = min(k, g)
                    for i in blk: o[i] = take/g
                    k -= take; j += g
            occ.append(o)
        return np.array(occ)
    mf.get_occ = get_occ
    E = mf.kernel()
    return E, mf.converged

if __name__ == '__main__':
    sym, q, cfg = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    basis = sys.argv[4] if len(sys.argv) > 4 else 'unc-ano-rcc'
    E, ok = run(sym, q, cfg, basis)
    import json
    print(json.dumps({"E": E, "converged": bool(ok), "basis": basis}))
