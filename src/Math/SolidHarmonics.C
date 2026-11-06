// File: Math/SolidHarmonics.C  Real spherical harmonics and their Cartesian-Gaussian forms.
//
// A separable or local field on a Gaussian basis is integrated analytically by writing it as Cartesian
// Gaussians, so the same polynomials serve every structure: the periodic KB and +U projections
// (GPW_Evaluator) and the molecular PP integrals (Gaussian::EOrbital_PP_IBS).  Moved here from
// GPW_Evaluator.C so the two cannot drift.
module;
#include <cmath>
#include <map>
#include <stdexcept>
#include <vector>
export module qchem.Math.SolidHarmonics;
import qchem.Math;            // Pi
import qchem.Math.Angular;    // SphericalShell, CartTerm, Monomial
import qchem.Types;           // rvec3_t

export namespace qchem::Math
{

//! \brief A primitive Cartesian-Gaussian scalar function,
//! \f$g(r)=\sum_t c_t\,(r-C)^{m_t}\,e^{-\alpha|r-C|^2}\f$ (one shared exponent, a finite Cartesian-
//! monomial polynomial \f$(r-C)^{m}=x^{m_x}y^{m_y}z^{m_z}\f$ about the centre \f$C\f$) -- the family of
//! scalar functions a Gaussian basis can integrate ANALYTICALLY.  Pure function-language: the caller may
//! think of \f$g\f$ as a projector, a moment, anything.
struct CartGaussian
{
    rvec3_t center;               //!< \f$C\f$
    double  alpha;                //!< the shared exponent \f$\alpha\f$
    std::vector<CartTerm> terms;  //!< \f$\{(m_t,c_t)\}\f$ -- the polynomial about \f$C\f$
};

//! \brief Real (tesseral) spherical harmonic \f$Y_{lm}(\hat r)\f$, unit-normalised on the sphere, from the
//! UNIT vector \a x,y,z.  Only \f$\sum_m|Y_{lm}\rangle\langle Y_{lm}|\f$ enters a KB term, so any orthonormal
//! real degree-\a l set is equivalent; these are the standard Cartesian forms.  GTH/HGH projectors reach
//! \f$l=3\f$ at most, and a larger \a l throws.
double RealYlm(int l, int m, double x, double y, double z)
{
    using std::sqrt;
    const double pi=Pi;
    switch (l)
    {
    case 0: return 0.5/sqrt(pi);                                            // Y_00
    case 1:
        switch (m)
        {
        case -1: return sqrt(0.75/pi)*y;
        case  0: return sqrt(0.75/pi)*z;
        case  1: return sqrt(0.75/pi)*x;
        }
        break;
    case 2:
        switch (m)
        {
        case -2: return sqrt(15.0/(4*pi))*x*y;
        case -1: return sqrt(15.0/(4*pi))*y*z;
        case  0: return sqrt( 5.0/(16*pi))*(3*z*z-1.0);                     // x^2+y^2+z^2=1
        case  1: return sqrt(15.0/(4*pi))*x*z;
        case  2: return sqrt(15.0/(16*pi))*(x*x-y*y);
        }
        break;
    case 3:
        switch (m)
        {
        case -3: return sqrt( 35.0/(32*pi))*y*(3*x*x-y*y);
        case -2: return sqrt(105.0/(4 *pi))*x*y*z;
        case -1: return sqrt( 21.0/(32*pi))*y*(5*z*z-1.0);
        case  0: return sqrt(  7.0/(16*pi))*z*(5*z*z-3.0);
        case  1: return sqrt( 21.0/(32*pi))*x*(5*z*z-1.0);
        case  2: return sqrt(105.0/(16*pi))*z*(x*x-y*y);
        case  3: return sqrt( 35.0/(32*pi))*x*(x*x-3*y*y);
        }
        break;
    }
    throw std::runtime_error("RealYlm: unsupported (l,m)");
}

//! \brief The Cartesian-monomial expansion of \f$r^lY_{lm}(\hat r)\f$, a homogeneous degree-\a l polynomial.
//! \c SphericalShell supplies the RAW monomial shape; the absolute scale is fixed NUMERICALLY against
//! \c RealYlm (evaluate both at a generic unit direction and take the ratio, exact to machine precision
//! since they are the same polynomial, and checked at a second direction).  POSITIVE scale: the harmonic is
//! \c SphericalShell's raw polynomial times a positive unit-sphere norm, the convention an atom-centred set
//! declares and a \c SphericalShellRep shell describes; so this is \c RealYlm up to a sign per \a m, which
//! the KB sum \f$|b|^2\f$ over \a m cannot see.
std::vector<CartTerm> YlmCartesian(int l, int m)
{
    std::vector<CartTerm> terms=SphericalShell(l)[size_t(l+m)];
    auto poly=[&terms](double x, double y, double z)
    {
        double s=0.0;
        for (const auto& t : terms)
            s += t.c * std::pow(x,t.p.n) * std::pow(y,t.p.l) * std::pow(z,t.p.m);
        return s;
    };
    auto fix=[&](double x, double y, double z)                     // RealYlm(unit)*r^l / raw(x,y,z)
    {
        double r=std::sqrt(x*x+y*y+z*z);
        return RealYlm(l,m, x/r,y/r,z/r) * std::pow(r,l) / poly(x,y,z);
    };
    const double N=fix(0.3,0.5,0.7);                               // generic direction (no raw zeros, l<=3)
    if (std::fabs(fix(0.9,-0.2,0.4)-N)>1e-10*std::fabs(N)+1e-14)
        throw std::runtime_error("YlmCartesian: raw solid harmonic is not proportional to RealYlm (convention drift)");
    for (auto& t : terms) t.c*=std::fabs(N);
    return terms;
}

//! Multiply a Cartesian polynomial by \f$(x^2+y^2+z^2)^n\f$ (the \f$r^{2n}\f$ factor of a higher radial
//! term), combining duplicate monomials.
std::vector<CartTerm> MultiplyR2(std::vector<CartTerm> terms, int n)
{
    for (int k=0;k<n;k++)
    {
        std::map<Monomial,double> acc;
        for (const auto& t : terms)
            for (int ax=0; ax<3; ax++)
            {
                Monomial p=t.p; p[ax]+=2;
                acc[p]+=t.c;
            }
        terms.clear();
        for (const auto& [p,c] : acc) terms.push_back({p,c});
    }
    return terms;
}

} //namespace
