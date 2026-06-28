// File: BasisSet/Atom/Evaluators/Internal/Wigner3j.C  Compile-time lookup table of Wigner 3j symbols.
//
// The 3j symbols for the (integer-l) angular ERI integrals are built at COMPILE time by the constexpr
// Racah formula -- no runtime "Initializing..." pass, no mutable global table, and no dependency on the
// wignerSymbols submodule (it remains only for the half-integer relativistic path + as a test oracle).
//
// The trick that keeps it constexpr WITHOUT a constexpr std::sqrt (not available before C++26): a 3j is
// phase * sqrt(rational) * sum(rational), so its SQUARE and its SIGN are both pure rationals/integers --
// computable at compile time from the constexpr factorial table.  We therefore store {(3j)^2, sign} and
// reconstruct the signed value as sign*sqrt(sq) at lookup.  The angular integrals consume mostly (3j)^2
// (the exchange terms), so they read sq() with NO sqrt at all; only the direct term (a product of four
// DISTINCT 3j's, sign-sensitive) takes the runtime sqrt via operator().
module;
#include <array>
export module qchem.BasisSet.Atom.Evaluators.Internal.Wigner3j;
import qchem.Factorials;   // qchem::Fact (constexpr factorial table)

namespace w3j_detail
{
    constexpr int phase(int n) {return (n&1)?-1:1;}          // (-1)^n (parity-correct for n<0 too)
    constexpr int iabs(int x) {return x<0?-x:x;}
    constexpr int imin(int a,int b) {return a<b?a:b;}
    constexpr int imax(int a,int b) {return a>b?a:b;}

    //! A 3j symbol stored as its square (rational) + its sign: value = sgn*sqrt(sq).
    struct W { double sq; signed char sgn; };

    //! Wigner 3j (integer arguments) via the Racah single-sum formula, returning {square, sign} so the
    //! whole computation stays rational (no sqrt).  Zero (sgn=0) outside the selection rules.
    constexpr W wig3j(int j1,int j2,int j3,int m1,int m2,int m3)
    {
        if (m1+m2+m3!=0)                                   return {0.0,0};
        if (j3<iabs(j1-j2) || j3>j1+j2)                    return {0.0,0};
        if (iabs(m1)>j1 || iabs(m2)>j2 || iabs(m3)>j3)     return {0.0,0};

        const auto& F=qchem::Fact;
        double delta=F[j1+j2-j3]*F[j1-j2+j3]*F[-j1+j2+j3]/F[j1+j2+j3+1];     // triangle coefficient
        double prod =F[j1+m1]*F[j1-m1]*F[j2+m2]*F[j2-m2]*F[j3+m3]*F[j3-m3];

        int tmin=imax(0,imax(j2-j3-m1,j1-j3+m2));
        int tmax=imin(j1+j2-j3,imin(j1-m1,j2+m2));
        double S=0.0;                                                        // Racah sum (rational)
        for (int t=tmin;t<=tmax;t++)
            S+=phase(t)/(F[t]*F[j1+j2-j3-t]*F[j1-m1-t]*F[j2+m2-t]*F[j3-j2+m1+t]*F[j3-j1-m2+t]);

        int sgn=phase(j1-j2-m3)*(S>0?1:(S<0?-1:0));                          // sqrt(delta*prod) is >=0
        return {delta*prod*S*S, static_cast<signed char>(sgn)};
    }

    constexpr int LMax=4;
    constexpr int NL=LMax+1, NK=2*LMax+1, NM=2*LMax+1;
    constexpr int Size=NL*NL*NK*NM*NM;
    constexpr int Index(int la,int lb,int k,int ma,int mb)
        {return (((la*NL+lb)*NK+k)*NM+(ma+LMax))*NM+(mb+LMax);}

    // The whole table, built once at compile time.  (The all-m=0 case used by the 3-argument lookups is
    // just ma=mb=0, so a single table serves both -- no separate m=0 table needed.)
    constexpr std::array<W,Size> MakeTable()
    {
        std::array<W,Size> a{};
        for (int la=0;la<=LMax;la++)
            for (int lb=0;lb<=LMax;lb++)
                for (int k=0;k<=2*LMax;k++)
                    for (int ma=-la;ma<=la;ma++)
                        for (int mb=-lb;mb<=lb;mb++)
                            a[Index(la,lb,k,ma,mb)]=wig3j(la,lb,k,ma,mb,-ma-mb);
        return a;
    }
    inline constexpr std::array<W,Size> Table=MakeTable();
}

//
//  Serve up Wigner 3j symbols for the angular ERI integrations:  (l_a l_b k; m_a m_b -m_a-m_b).
//
export class Wigner3j
{
public:
    static const Wigner3j w3j; //Static instance with short name for convenience.

    //! The signed 3j value (runtime sqrt).  Used by the sign-sensitive direct angular term.
    double operator()(int la, int lb, int k) const {return (*this)(la,lb,k,0,0);}
    double operator()(int la, int lb, int k, int ma, int mb) const;
    //! The 3j SQUARE -- no sqrt.  Used by the exchange angular terms (which square the symbol anyway).
    double sq(int la, int lb, int k) const {return sq(la,lb,k,0,0);}
    double sq(int la, int lb, int k, int ma, int mb) const;
private:
    Wigner3j()=default;
    static constexpr int LMax=w3j_detail::LMax;
};
