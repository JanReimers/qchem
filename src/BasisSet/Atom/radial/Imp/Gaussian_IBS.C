// File: BasisSet/Atom/radial/Imp/Gaussian_IBS.C
module;
#include <valarray>
#include <cmath>
#include <cassert>
module BasisSet.Atom.Gaussian_IBS;
import qchem.BasisSet.Atom.Internal.radial.GaussianRk;
import qchem.BasisSet.Atom.Internal.radial.GaussianIntegrals;
import Common.Constants;
import Common.IntPow;


inline double Overlap(double ea , double eb,size_t l_total)
{
    return Gaussian::Integral(ea+eb,l_total); //Already has 4*Pi and r^2 from dr.
}

inline double Grad2(double ea , double eb,size_t la, size_t lb)
{
    assert(la==lb);
    double t=ea+eb;
    size_t l1=la+1;
    return  l1*l1     * Gaussian::Integral(t,2*la-2)
            -2*l1 * t * Gaussian::Integral(t,2*la  )
            +4*ea*eb  * Gaussian::Integral(t,2*la+2);
}

inline double Inv_r1(double ea , double eb,size_t l_total)
{
    return Gaussian::Integral(ea+eb,l_total-1); //Already has 4*Pi
}
inline double Inv_r2(double ea , double eb,size_t l_total)
{
    return Gaussian::Integral(ea+eb,l_total-2); //Already has 4*Pi
}

inline double Repulsion(double eab, double ec,size_t la,size_t lc)
{    
    Gaussian::RkEngine cd(eab,ec,std::max(la,lc));
    return FourPi2*cd.Coulomb_R0(la,lc);
}

inline double Charge(double ea, size_t l)
{
    return ::Gaussian::Integral(ea,l);
}

Gaussian_IBS::ds_t Gaussian_IBS::norms() const
{
    ds_t ret(size());
    for (size_t i=0;i<size();i++) ret[i]=1.0/sqrt(::Overlap(es[i],es[i],2*l)); 
    return ret;
}

Gaussian_IBS::omls_t Gaussian_IBS::Overlap() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Overlap(es[i-1],es[j-1],2*l)*ns[i-1]*ns[j-1];

    return S;
}

Gaussian_IBS::omls_t Gaussian_IBS::Grad2() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Grad2(es[i-1],es[j-1],l,l)*ns[i-1]*ns[j-1];

    return S;
}

Gaussian_IBS::omls_t Gaussian_IBS::Inv_r1() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Inv_r1(es[i-1],es[j-1],2*l)*ns[i-1]*ns[j-1];

    return S;
}

Gaussian_IBS::omls_t Gaussian_IBS::Inv_r2() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Inv_r2(es[i-1],es[j-1],2*l)*ns[i-1]*ns[j-1];

    return S;
}

Gaussian_IBS::omls_t Gaussian_IBS::Repulsion() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Repulsion(es[i-1],es[j-1],l,l)*ns[i-1]*ns[j-1];

    return S;
}

Gaussian_IBS::omlv_t Gaussian_IBS::Charge() const
{
    omlv_t V(size());
    for (auto i:V.indices())
            V(i)= ::Charge(es[i-1],l)*ns[i-1];

    return V;
}

template <class T> Vector<T> convert(const std::valarray<T>& v) 
{
    Vector<T> ret(v.size());
    size_t i=0;
    for (auto vi:v) ret(i)=vi;
    return ret;
}

template <class v> v gaussian(double r,size_t l,const v& e, const v& n)
{
    return n*uintpow(r,l)*exp(-e*r*r);
}
template <class v> v grad_gaussian(double r,size_t l,const v& e, const v& n)
{
    double lr= r==0 ? 0 : l/r;
    return (lr-2*r*e)*gaussian(r,l,e,n);
}

Gaussian_IBS::Vec    Gaussian_IBS::operator() (const RVec3& r) const
{
    return convert(gaussian(norm(r),l,es,ns));
}

Gaussian_IBS::Vec3Vec Gaussian_IBS::Gradient(const RVec3& r) const
{
    ds_t grad=grad_gaussian(norm(r),l,es,ns);
    RVec3 rhat=r/norm(r);
    Vec3Vec ret(size());
    size_t i=0;
    for (auto& g:grad) ret(++i)=g*rhat;
    return ret;
}