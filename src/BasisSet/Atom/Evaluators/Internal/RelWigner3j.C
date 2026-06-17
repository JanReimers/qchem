// File: BasisSet/Atom/Evaluators/Internal/RelWigner3j.C  Wigner 3j symbols with half-integer arguments.
module;
export module qchem.BasisSet.Atom.Evaluators.Internal.RelWigner3j;

//
//  Pre-computed Wigner 3j symbols with half-integer j arguments, for relativistic
//  (SphericalSpinor/RKB) angular integrals.  Two symbol types are stored:
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
    static RelWigner3j w3j; //Static instance with short name for convenience.
    double operator()(int κa, int κb, int k) const;
    double operator()(int κa, int κb, int k, double mja, double mjb) const;
private:
    RelWigner3j();
    static const int LMax  = 4;
    static const int KMax  = LMax+1;       //max |κ|
    static const int KkMax = 2*LMax+1;     //max k: ja+jb <= 2*(LMax+0.5) = 2*LMax+1
    static const int MJMax = 2*KMax-1;     //max value of 2*mj = 2*(LMax+0.5) = 2*LMax+1

    int κ_idx (int κ)       const {return κ + KMax;}
    int mj_idx(double mj)   const {return (int)(2*mj) + MJMax;}

    // Data[κa_idx][κb_idx][k]
    double Data  [2*KMax+1][2*KMax+1][KkMax+1];
    // Data_m[κa_idx][κb_idx][k][mja_idx][mjb_idx]
    double Data_m[2*KMax+1][2*KMax+1][KkMax+1][2*MJMax+1][2*MJMax+1];
};
