// File: BasisSet/Atom/radial/Imp/Slater_IBS.C
module;
#include <valarray>
#include <cmath>
#include <cassert>
module BasisSet.Atom.Slater_IBS;
import qchem.BasisSet.Atom.Internal.radial.Slater.Rk;
import qchem.BasisSet.Atom.Internal.radial.Slater.Integrals;
import Common.Constants;
import Common.IntPow;


inline double Overlap(double ea , double eb,size_t l_total)
{
    return Slater::Integral(ea+eb,l_total); //Already has 4*Pi and r^2 from dr.
}

inline double Grad2(double ea , double eb,size_t la, size_t lb)
{
    assert(la==lb);
    double ab=ea+eb;
    int l=la; //Safer to do formulas with int.
    // int ll=l*(l+1);
    double Term1=(l+1)*(l+1)*Slater::Integral(ab,2*l-2); //SlaterIntegral already has 4*Pi
    double Term2=-(l+1)*ab* Slater::Integral(ab,2*l-1);
    double Term3=ea*eb*Slater::Integral(ab,2*l);
    return Term1+Term2+Term3;
}

inline double Inv_r1(double ea , double eb,size_t l_total)
{
    return Slater::Integral(ea+eb,l_total-1); //Already has 4*Pi
}
inline double Inv_r2(double ea , double eb,size_t l_total)
{
    return Slater::Integral(ea+eb,l_total-2); //Already has 4*Pi
}

inline double Repulsion(double eab, double ec,size_t la,size_t lc)
{    
    Slater::RkEngine cd(eab,ec,std::max(la,lc));
    return FourPi2*cd.Coulomb_R0(la,lc);
}

inline double Charge(double ea, size_t l)
{
    return ::Slater::Integral(ea,l);
}

//---------------------------------------------------------------------------
//
//  Start member functions.
//
void Slater_IBS::Register(ExponentGrouper& grouper)
{
    for (auto e:es) es_indices.push_back(grouper.Insert(e,l));
}


Slater_IBS::ds_t Slater_IBS::norms() const
{
    ds_t ret(size());
    for (size_t i=0;i<size();i++) ret[i]=1.0/sqrt(::Overlap(es[i],es[i],2*l)); 
    return ret;
}

Slater_IBS::omls_t Slater_IBS::Overlap() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Overlap(es[i-1],es[j-1],2*l)*ns[i-1]*ns[j-1];

    return S;
}

Slater_IBS::omls_t Slater_IBS::Grad2() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Grad2(es[i-1],es[j-1],l,l)*ns[i-1]*ns[j-1];

    return S;
}

Slater_IBS::omls_t Slater_IBS::Inv_r1() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Inv_r1(es[i-1],es[j-1],2*l)*ns[i-1]*ns[j-1];

    return S;
}

Slater_IBS::omls_t Slater_IBS::Inv_r2() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Inv_r2(es[i-1],es[j-1],2*l)*ns[i-1]*ns[j-1];

    return S;
}

Slater_IBS::omls_t Slater_IBS::Repulsion() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Repulsion(es[i-1],es[j-1],l,l)*ns[i-1]*ns[j-1];

    return S;
}

Slater_IBS::omlv_t Slater_IBS::Charge() const
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

template <class v> v slater(double r,size_t l,const v& e, const v& n)
{
    return n*uintpow(r,l)*exp(-e*r);
}
template <class v> v grad_slater(double r,size_t l,const v& e, const v& n)
{
    double lr= r==0 ? 0 : l/r;
    return (lr-e)*slater(r,l,e,n);
}

Slater_IBS::Vec    Slater_IBS::operator() (const RVec3& r) const
{
    return convert(slater(norm(r),l,es,ns));
}

Slater_IBS::Vec3Vec Slater_IBS::Gradient(const RVec3& r) const
{
    ds_t grad=grad_slater(norm(r),l,es,ns);
    RVec3 rhat=r/norm(r);
    Vec3Vec ret(size());
    size_t i=0;
    for (auto& g:grad) ret(++i)=g*rhat;
    return ret;
}