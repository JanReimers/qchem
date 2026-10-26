// File: BasisSet/Imp/AtomCentredSet.C  The atom-centred set's value type: radial normalisation.
module;
#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>
module qchem.BasisSet.AtomCentredSet;
import qchem.Math.Angular;   // Math::SphericalShell -- the set's harmonics (BareCoulomb's Gaunt factors)

namespace qchem::BasisSet
{

// int_0^inf r^{2l+2+2(n+n')} e^{-(a+a') r^2} dr = Gamma(p)/(2 (a+a')^p),  p = l + n + n' + 3/2.
double GaussianAtomCentredSet::RadialNorm2(int l, const std::vector<Math::Gaussian>& radial)
{
    double s=0.0;
    for (const auto& a : radial)
        for (const auto& b : radial)
        {
            const double p=l+a.n+b.n+1.5;
            s+=a.c*b.c*std::tgamma(p)/(2.0*std::pow(a.alpha+b.alpha,p));
        }
    return s;
}

size_t GaussianAtomCentredSet::Add(const rvec3_t& centre, int l, std::vector<Math::Gaussian> radial)
{
    if (l<0) throw std::invalid_argument("GaussianAtomCentredSet::Add: l must be >= 0");
    if (radial.empty()) throw std::invalid_argument("GaussianAtomCentredSet::Add: empty radial");
    const double n2=RadialNorm2(l,radial);
    if (!(n2>0.0)) throw std::invalid_argument("GaussianAtomCentredSet::Add: radial has zero norm");
    const double s=1.0/std::sqrt(n2);
    for (auto& t : radial) t.c*=s;
    itsM.push_back({centre,l,std::move(radial)});
    return itsM.size()-1;
}

//------------------------------------------------------------------------------ the set's own bare integrals
namespace
{
// int over the unit sphere of x^a y^b z^c: 0 unless all even, else 2 G((a+1)/2) G((b+1)/2) G((c+1)/2) / G((a+b+c+3)/2).
double SphereMoment(int a, int b, int c)
{
    if (a%2 || b%2 || c%2) return 0.0;
    return 2.0*std::tgamma(0.5*(a+1))*std::tgamma(0.5*(b+1))*std::tgamma(0.5*(c+1))/std::tgamma(0.5*(a+b+c+3));
}
using Harmonic = std::vector<Math::CartTerm>;
// The set's angular convention: SphericalShell(l)[l+m] at a POSITIVE scale to unit norm on the sphere.
std::vector<Harmonic> UnitHarmonics(int l)
{
    std::vector<Harmonic> h=Math::SphericalShell(l);
    for (Harmonic& y : h)
    {
        double n2=0.0;
        for (const auto& a : y) for (const auto& b : y) n2+=a.c*b.c*SphereMoment(a.p.n+b.p.n, a.p.l+b.p.l, a.p.m+b.p.m);
        for (auto& t : y) t.c/=std::sqrt(n2);
    }
    return h;
}
// The real Gaunt factor int Y_a Y_b Y_c dOmega (each a polynomial restricted to the unit sphere).
double Gaunt(const Harmonic& A, const Harmonic& B, const Harmonic& C)
{
    double g=0.0;
    for (const auto& a : A) for (const auto& b : B) for (const auto& c : C)
        g+=a.c*b.c*c.c*SphereMoment(a.p.n+b.p.n+c.p.n, a.p.l+b.p.l+c.p.l, a.p.m+b.p.m+c.p.m);
    return g;
}
// int_0^inf dr r^{Pout-k-1} e^{-bout r^2} int_0^r ds s^{Pin+k} e^{-bin s^2}: the s<r half of a Slater integral.  With
// s = r t the r-integral is a Gamma function and the t-integral, after u = bin t^2/(bout+bin t^2), the incomplete beta
// B_x(a,b), x = bin/(bout+bin), a = (Pin+k+1)/2, b = (Pout-k)/2 -- an INTEGER b >= 1 here (Pout and k share parity,
// k <= Pout-2), so B_x(a,b) = sum_j C(b-1,j) (-1)^j x^{a+j}/(a+j) exactly.
double SlaterHalf(int Pout, double bout, int Pin, double bin, int k)
{
    assert((Pout-k)%2==0 && Pout-k>=2 && "SlaterHalf: the outer density's power and k must share parity with b >= 1");
    const double c=0.5*(Pout+Pin+1), a=0.5*(Pin+k+1);
    const int b=(Pout-k)/2;
    const double x=bin/(bout+bin);
    double B=0.0, binom=1.0;
    for (int j=0;j<b;j++)
    {
        B+=(j%2 ? -binom : binom)*std::pow(x,a+j)/(a+j);
        binom*=double(b-1-j)/double(j+1);
    }
    return 0.25*std::tgamma(c)*std::pow(bout,a-c)*std::pow(bin,-a)*B;
}
// R^k between two radial DENSITIES r^{P1} e^{-b1 r^2} and r^{P2} e^{-b2 r^2} (the r^2 Jacobian included in P): the
// double integral with kernel r_<^k / r_>^{k+1}, both halves.
double SlaterPrimitive(int P1, double b1, int P2, double b2, int k)
{
    return SlaterHalf(P1,b1,P2,b2,k)+SlaterHalf(P2,b2,P1,b1,k);
}
}

ERI4Block GaussianAtomCentredSet::BareCoulomb(const std::vector<size_t>& cols) const
{
    const size_t m=cols.size();
    struct Fn { size_t M; int l; size_t mi; };
    std::vector<Fn> fn(m);
    for (size_t i=0;i<m;i++)
    {
        size_t M=0;
        while (M<Count() && cols[i]>=Offset(M)+Size(M)) M++;
        if (M==Count()) throw std::out_of_range("GaussianAtomCentredSet::BareCoulomb: column outside the set");
        fn[i]={M, L(M), cols[i]-Offset(M)};
    }
    for (size_t i=1;i<m;i++)
    {
        const rvec3_t d=Centre(fn[i].M)-Centre(fn[0].M);
        if (d.x*d.x+d.y*d.y+d.z*d.z>1e-20)
            throw std::invalid_argument("GaussianAtomCentredSet::BareCoulomb: the columns span more than one centre -- only one-centre integrals are built");
    }
    int lmax=0; for (const Fn& f : fn) lmax=std::max(lmax,f.l);
    std::vector<std::vector<Harmonic>> Y(2*lmax+1);              // unit harmonics for every l and k that can occur
    for (int l=0;l<=2*lmax;l++) Y[l]=UnitHarmonics(l);
    // The radial pair density R_a R_b r^2 as Gaussian terms {coefficient, power of r, exponent}.
    struct Term { double c; int P; double beta; };
    auto pairDensity=[&](const Fn& a, const Fn& b)
    {
        std::vector<Term> t;
        for (const auto& ga : itsM[a.M].radial) for (const auto& gb : itsM[b.M].radial)
            t.push_back({ga.c*gb.c, a.l+b.l+2+2*(ga.n+gb.n), ga.alpha+gb.alpha});
        return t;
    };
    std::vector<std::vector<std::vector<Term>>> rho(m, std::vector<std::vector<Term>>(m));
    for (size_t a=0;a<m;a++) for (size_t b=a;b<m;b++) rho[a][b]=pairDensity(fn[a],fn[b]);
    // G[a][b][k][q] = int Y_a Y_b Y_kq, zero by parity/triangle where it must be.
    auto gaunt=[&](size_t a, size_t b, int k, size_t q)
    { return Gaunt(Y[fn[a].l][fn[a].mi], Y[fn[b].l][fn[b].mi], Y[k][q]); };
    ERI4Block eri(m);
    for (size_t a=0;a<m;a++) for (size_t b=a;b<m;b++)
        for (size_t c=0;c<m;c++) for (size_t d=c;d<m;d++)
        {
            if (a*m+b>c*m+d) continue;                           // (ab|cd)=(cd|ab): each pair-of-pairs once
            const int lab=fn[a].l+fn[b].l, lcd=fn[c].l+fn[d].l;
            double v=0.0;
            for (int k=0;k<=std::min(lab,lcd);k++)
            {
                if ((lab+k)%2 || (lcd+k)%2) continue;            // the Gaunt factors vanish by parity
                double ang=0.0;
                for (size_t q=0;q<size_t(2*k+1);q++) ang+=gaunt(a,b,k,q)*gaunt(c,d,k,q);
                if (ang==0.0) continue;
                double Rk=0.0;
                for (const Term& s : rho[a][b]) for (const Term& t : rho[c][d]) Rk+=s.c*t.c*SlaterPrimitive(s.P,s.beta,t.P,t.beta,k);
                v+=4.0*M_PI/(2*k+1)*Rk*ang;
            }
            for (size_t s=0;s<2;s++) for (size_t t=0;t<2;t++)       // the 8-fold symmetry of real functions
            {
                const size_t i=s?b:a, j=s?a:b, p=t?d:c, r=t?c:d;
                eri(i,j,p,r)=v; eri(p,r,i,j)=v;
            }
        }
    return eri;
}

} //namespace
