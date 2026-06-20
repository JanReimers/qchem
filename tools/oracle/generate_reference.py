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

LMAX = 3                                   # s,p,d,f  (PG1 Hermite machinery is capped at f for now)

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

def main():
    data = {"meta": {"engine": "gbasis", "LMAX": LMAX,
                     "normalization": "per-cartesian-component unit self-overlap"}}
    print("Generating reference integrals (gbasis):")
    data["overlap"] = two_centre("overlap", overlap_integral)
    with open("tools/oracle/reference_integrals.json", "w") as f:
        json.dump(data, f)
    print("wrote tools/oracle/reference_integrals.json")

if __name__ == "__main__":
    main()
