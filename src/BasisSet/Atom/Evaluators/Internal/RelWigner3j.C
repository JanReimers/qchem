// File: BasisSet/Atom/Evaluators/Internal/RelWigner3j.C  Wigner 3j symbols with half-integer arguments.
module;
export module qchem.BasisSet.Atom.Evaluators.Internal.RelWigner3j;

namespace qchem {

//
//  Wigner 3j symbols with half-integer j arguments, for relativistic (SphericalSpinor/RKB) angular
//  integrals.  A thin, stateless wrapper over the home-grown half-integer Racah core (Wigner::wigner3j);
//  no precomputed table, no submodule, no startup pass.  Two symbol types:
//
//  Parity symbol:   (ja  k  jb)   where ja=j(κa), jb=j(κb)
//                   (½   0  -½)
//
//  m-dependent:     (ja   k   jb)
//                   (-mja  q  mjb)   where q = mja-mjb
//
export class RelWigner3j
{
public:
    static const RelWigner3j w3j; //Static instance with short name for convenience.
    double operator()(int κa, int κb, int k) const;
    double operator()(int κa, int κb, int k, double mja, double mjb) const;
private:
    RelWigner3j()=default;
    static const int LMax = 4;
    static const int KMax = LMax+1;       //max |κ|
    static const int KkMax= 2*LMax+1;     //max k: ja+jb <= 2*(LMax+0.5) = 2*LMax+1
};

} // namespace qchem