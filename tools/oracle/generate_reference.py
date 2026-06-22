#!/usr/bin/env python3
"""Generate reference primitive-Gaussian integrals (the oracle) for the PG1 M&D unit tests.

Independent engine: gbasis (theochem) -- pure-Python analytic integrals, normalized so each Cartesian
component has unit self-overlap (verified: S(c,c)=1).  PG1 computes UNnormalized integrals and applies
the same per-component normalization separately, so the C++ test compares
    GetNormalization(pa) * GetNormalization(pb) * Integrate(raw)   ==   gbasis[pa,pb].

Every matrix element is labelled by the explicit Cartesian powers (lx,ly,lz) of each function (from
gbasis's angmom_components_cart), so the C++ side matches by (n,l,m) and ordering never matters.

Output: reference_integrals.json (full double precision).  Capped at L<=3 (s,p,d,f); g (L=4) later.
"""
import json, itertools, numpy as np
from gbasis.contractions import GeneralizedContractionShell
from gbasis.integrals.overlap import overlap_integral
from gbasis.integrals.kinetic_energy import kinetic_energy_integral
from gbasis.integrals.point_charge import point_charge_integral
from gbasis.integrals.electron_repulsion import electron_repulsion_integral

LMAX = 3                                   # s,p,d,f  (PG1 Hermite machinery is capped at f for now)

# Fixed point charges for the nuclear-attraction sweep (shared with the C++ Structure).
NUCLEI = [{"Z": 6.0, "center": [0.0, 0.0, 0.4]}, {"Z": 2.0, "center": [-0.6, 0.5, -0.3]}]

def shell(alpha, center, L):
    return GeneralizedContractionShell(L, np.array(center, float), np.array([1.0]),
                                       np.array([float(alpha)]), 'cartesian')

def powers(L):
    return [tuple(int(x) for x in p) for p in shell(1.0,[0,0,0],L).angmom_components_cart]

# A deliberately varied set of primitive shells to thrash: exponents spanning ~1.5 decades, and three
# centre arrangements (coincident, axis-displaced, general).
EXPS    = [0.30, 1.0, 3.7]
CENTERS = [[0.0,0.0,0.0], [0.0,0.0,1.3], [0.7,-1.1,0.4]]

def shell_json(alpha, center, L):
    return {"alpha": alpha, "center": list(center), "L": L}

def two_centre(name, integral_fn):
    """integral_fn([sh_a, sh_b]) -> normalized matrix; emit per-element records labelled by powers."""
    out = []
    for La, Lb in itertools.product(range(LMAX+1), repeat=2):
        pa, pb = powers(La), powers(Lb)
        for ea, eb in itertools.product(EXPS, repeat=2):
            for ca, cb in itertools.product(range(len(CENTERS)), repeat=2):
                if ca == cb and ea == eb and La == Lb:
                    pass  # keep self/diagonal cases too -- they exercise normalization hardest
                A, B = CENTERS[ca], CENTERS[cb]
                # overlap_integral([a,b]) is the COMBINED-basis matrix [[Saa,Sab],[Sba,Sbb]];
                # we want the a-b cross block Sab = M[0:na, na:na+nb].
                M = integral_fn([shell(ea, A, La), shell(eb, B, Lb)])
                na = len(pa)
                elems = []
                for i, fa in enumerate(pa):
                    for j, fb in enumerate(pb):
                        elems.append({"a": list(fa), "b": list(fb), "value": float(M[i, na + j])})
                out.append({"a": shell_json(ea, A, La), "b": shell_json(eb, B, Lb), "elements": elems})
    print(f"  {name}: {len(out)} shell-pairs, "
          f"{sum(len(r['elements']) for r in out)} elements")
    return out

# 4-centre ERI: a bounded sweep (O(N^4) elements), distinct fixed exponents/centres per slot, L up to
# LERI.  electron_repulsion_integral([a,b,c,d]) returns the combined-basis tensor; we want the
# a,b,c,d cross block E[i, na+j, na+nb+k, na+nb+nc+l].
LERI = 2                                    # s,p,d for the ERI sweep (dddd already exercises a lot)
ERI_EXP    = [1.1, 0.9, 1.3, 0.7]           # one exponent per slot a,b,c,d
ERI_CENTER = [[0,0,0],[0,0,1.1],[0.6,0,0.3],[0.2,-0.8,0.5]]

def eri():
    out = []
    for La,Lb,Lc,Ld in itertools.product(range(LERI+1), repeat=4):
        ps = [powers(La),powers(Lb),powers(Lc),powers(Ld)]
        sh = [shell(ERI_EXP[s], ERI_CENTER[s], [La,Lb,Lc,Ld][s]) for s in range(4)]
        E  = electron_repulsion_integral(sh, notation='chemist')  # (ab|cd): a,b on elec 1 -- matches PG1
        n  = [len(p) for p in ps]
        off= [0, n[0], n[0]+n[1], n[0]+n[1]+n[2]]
        elems=[]
        for i,fa in enumerate(ps[0]):
         for j,fb in enumerate(ps[1]):
          for k,fc in enumerate(ps[2]):
           for l,fd in enumerate(ps[3]):
            v = float(E[off[0]+i, off[1]+j, off[2]+k, off[3]+l])
            elems.append({"a":list(fa),"b":list(fb),"c":list(fc),"d":list(fd),"value":v})
        out.append({"a":shell_json(ERI_EXP[0],ERI_CENTER[0],La),"b":shell_json(ERI_EXP[1],ERI_CENTER[1],Lb),
                    "c":shell_json(ERI_EXP[2],ERI_CENTER[2],Lc),"d":shell_json(ERI_EXP[3],ERI_CENTER[3],Ld),
                    "elements":elems})
    print(f"  eri: {len(out)} quads, {sum(len(r['elements']) for r in out)} elements")
    return out

# --- 3-centre (DFT): pyscf int3c1e (<ab|c> overlap) and int3c2e ((ab|c) Coulomb). -----------------
# gbasis has no 3-centre, so use pyscf/libcint and renormalize per-component via pyscf's OWN overlap
# diagonal (verified to reproduce the per-component unit-self-overlap convention to ~1e-16), so the
# libcint Cartesian normalization never leaks in.
from pyscf import gto
C3_EXP    = [1.1, 0.9, 1.3]
C3_CENTER = [[0.0,0.0,0.0],[0.0,0.0,1.1],[0.6,0.0,0.3]]

def _mol3(La,Lb,Lc):
    m = gto.Mole()
    m.atom  = [["H1",tuple(C3_CENTER[0])],["H2",tuple(C3_CENTER[1])],["H3",tuple(C3_CENTER[2])]]
    m.basis = {"H1":[[La,[C3_EXP[0],1.0]]],"H2":[[Lb,[C3_EXP[1],1.0]]],"H3":[[Lc,[C3_EXP[2],1.0]]]}
    m.cart=True; m.unit="Bohr"; m.spin=1; m.build(verbose=0)  # 3 H = 3 electrons (odd)
    return m

def three_centre(name, intor):
    out=[]
    for La,Lb,Lc in itertools.product(range(LMAX+1), repeat=3):
        m=_mol3(La,Lb,Lc)
        dnorm=np.sqrt(np.diag(m.intor("int1e_ovlp_cart")))     # per-component self-overlap
        M=m.intor(intor)/dnorm[:,None,None]/dnorm[None,:,None]/dnorm[None,None,:]
        pa,pb,pc=powers(La),powers(Lb),powers(Lc)
        na,nb=len(pa),len(pb)
        elems=[]
        for i,fa in enumerate(pa):
         for j,fb in enumerate(pb):
          for k,fc in enumerate(pc):
            elems.append({"a":list(fa),"b":list(fb),"c":list(fc),
                          "value":float(M[i, na+j, na+nb+k])})
        out.append({"a":shell_json(C3_EXP[0],C3_CENTER[0],La),"b":shell_json(C3_EXP[1],C3_CENTER[1],Lb),
                    "c":shell_json(C3_EXP[2],C3_CENTER[2],Lc),"elements":elems})
    print(f"  {name}: {len(out)} triples, {sum(len(r['elements']) for r in out)} elements")
    return out

def main():
    data = {"meta": {"engine": "gbasis", "LMAX": LMAX,
                     "normalization": "per-cartesian-component unit self-overlap"}}
    print("Generating reference integrals (gbasis):")
    data["overlap"] = two_centre("overlap", overlap_integral)
    data["kinetic"] = two_centre("kinetic", kinetic_energy_integral)
    ncoords  = np.array([n["center"] for n in NUCLEI], float)
    ncharges = np.array([n["Z"]      for n in NUCLEI], float)
    data["nuclei"]  = NUCLEI
    # point_charge_integral returns per-charge contributions (nbf,nbf,npts), already attractive
    # (negative for +Z), matching PG1's Nuclear convention; sum over the charges.
    data["nuclear"] = two_centre("nuclear", lambda b: point_charge_integral(b, ncoords, ncharges).sum(axis=-1))
    data["eri"]     = eri()
    data["overlap3c"]   = three_centre("overlap3c",   "int3c1e_cart")  # <ab|c>
    data["repulsion3c"] = three_centre("repulsion3c", "int3c2e_cart")  # (ab|c)
    with open("tools/oracle/reference_integrals.json", "w") as f:
        json.dump(data, f)
    print("wrote tools/oracle/reference_integrals.json")

if __name__ == "__main__":
    main()
