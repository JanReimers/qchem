// File: BasisSet/Atom/Evaluators/Internal/Wigner3j.C  Compile-time Wigner 3j symbols (home-grown).
//
// The Wigner 3j / Clebsch-Gordan symbols for the angular ERI integrals are computed from the Racah single-
// sum formula -- no wignerSymbols submodule, no runtime "Initializing..." pass, no mutable global table.
//
// The core works in DOUBLED-INTEGER arguments (tj = 2j, tm = 2m), so half-integer j (the relativistic
// SphericalSpinor path) is exact integer arithmetic too -- one core serves both the integer-l and the
// half-integer-j paths.  And the trick that keeps it constexpr WITHOUT a constexpr std::sqrt (not available
// before C++26): a 3j is phase*sqrt(rational)*sum(rational), so its SQUARE and SIGN are pure
// rationals/integers, computable at compile time from the constexpr factorial table.  We store {(3j)^2,
// sign} and reconstruct the signed value as sign*sqrt(sq) at lookup (the exchange angular terms consume
// (3j)^2 directly, with no sqrt at all).
module;
#include <array>
#include <cmath>      // std::sqrt, std::lround (runtime, for the signed-value helpers)
export module qchem.BasisSet.Atom.Evaluators.Internal.Wigner3j;
import qchem.Factorials;   // qchem::Fact (constexpr factorial table)

export namespace Wigner
{
    constexpr int phase(int n) {return (n&1)?-1:1;}          // (-1)^n (parity-correct for n<0 too)
    constexpr int iabs(int x) {return x<0?-x:x;}
    constexpr int imin(int a,int b) {return a<b?a:b;}
    constexpr int imax(int a,int b) {return a>b?a:b;}

    //! A 3j symbol stored as its square (rational) + its sign: value = sgn*sqrt(sq).
    struct W { double sq; signed char sgn; };

    //! Wigner 3j via the Racah single-sum formula in DOUBLED-INTEGER arguments (tj=2j, tm=2m), returning
    //! {square, sign} so the whole computation stays rational (no sqrt).  Zero (sgn=0) outside the
    //! selection rules.  Valid for integer or half-integer j (any tj/tm with the physical parities).
    constexpr W wig3j2(int tj1,int tj2,int tj3,int tm1,int tm2,int tm3)
    {
        if (tm1+tm2+tm3!=0)                                       return {0.0,0};
        if (tj3<iabs(tj1-tj2) || tj3>tj1+tj2)                     return {0.0,0};
        if (((tj1+tj2+tj3)&1)!=0)                                 return {0.0,0};   // j-sum must be integer
        if (iabs(tm1)>tj1 || iabs(tm2)>tj2 || iabs(tm3)>tj3)      return {0.0,0};
        if (((tj1+tm1)&1) || ((tj2+tm2)&1) || ((tj3+tm3)&1))      return {0.0,0};   // m,j same half-parity

        // Factorial arguments are all integers (the doubled sums are even -> exact /2).
        int a=(tj1+tj2-tj3)/2, b=(tj1-tj2+tj3)/2, c=(-tj1+tj2+tj3)/2, n=(tj1+tj2+tj3)/2+1;
        int j1mm1=(tj1-tm1)/2, j2pm2=(tj2+tm2)/2;
        const auto& F=qchem::Fact;
        double delta=F[a]*F[b]*F[c]/F[n];
        double prod =F[(tj1+tm1)/2]*F[j1mm1]*F[j2pm2]*F[(tj2-tm2)/2]*F[(tj3+tm3)/2]*F[(tj3-tm3)/2];

        int tmin=imax(0,imax((tj2-tj3-tm1)/2,(tj1-tj3+tm2)/2));
        int tmax=imin(a,imin(j1mm1,j2pm2));
        double S=0.0;                                                          // Racah sum (rational)
        for (int t=tmin;t<=tmax;t++)
            S+=phase(t)/(F[t]*F[a-t]*F[j1mm1-t]*F[j2pm2-t]*F[(tj3-tj2+tm1)/2+t]*F[(tj3-tj1-tm2)/2+t]);

        int sgn=phase((tj1-tj2-tm3)/2)*(S>0?1:(S<0?-1:0));                     // sqrt(delta*prod) >= 0
        return {delta*prod*S*S, static_cast<signed char>(sgn)};
    }

    //! Doubled-integer argument from a (half-)integer double: 2*x rounded to the nearest int (exact).
    inline int d2(double x) {return static_cast<int>(std::lround(2.0*x));}

    //! Signed 3j value for (possibly half-integer) double arguments (runtime sqrt).
    inline double wigner3j(double j1,double j2,double j3,double m1,double m2,double m3)
    {
        W w=wig3j2(d2(j1),d2(j2),d2(j3),d2(m1),d2(m2),d2(m3));
        return w.sgn*std::sqrt(w.sq);
    }

    //! Clebsch-Gordan <j1 m1 j2 m2 | J M> = (-1)^(j1-j2+M) sqrt(2J+1) (j1 j2 J; m1 m2 -M).  The (2J+1) is
    //! folded UNDER the single sqrt (sqrt((2J+1)*sq)) so common cases stay exact -- e.g. the s_1/2 CG is
    //! exactly 1 (2*1/2=1), not 1+-few ULP from sqrt(2)*sqrt(1/2).
    inline double clebschGordan(double j1,double j2,double J,double m1,double m2,double M)
    {
        W w=wig3j2(d2(j1),d2(j2),d2(J),d2(m1),d2(m2),d2(-M));
        return phase((int)std::lround(j1-j2+M))*w.sgn*std::sqrt((2*J+1)*w.sq);
    }
}

//
//  Integer-l 3j lookup table for the angular ERI integrations:  (l_a l_b k; m_a m_b -m_a-m_b).
//
namespace w3j_detail
{
    using Wigner::W;
    constexpr W wig3j(int j1,int j2,int j3,int m1,int m2,int m3)   // integer-l convenience (doubled inside)
        {return Wigner::wig3j2(2*j1,2*j2,2*j3,2*m1,2*m2,2*m3);}

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
//  Serve up Wigner 3j symbols for the (integer-l) angular ERI integrations.
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
