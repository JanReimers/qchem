// File: BasisSet1/Imp/Orbital_DFT_IBS.C
module;
#include <cassert>
#include <blaze/Math.h>
module qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.Internal.DB_Cache;

namespace BasisSet
{
template <class T> const ERI3<T>& Orbital_DFT_IBS<T>::Overlap3C  (const Fit_IBS& c) const
{
    auto cache=theGlobalCache;
    assert(cache);
    IntegralsCache_Base::IBS_ID_t abid(RadialID(),AngularID());
    IntegralsCache_Base::IBS_ID_t cid(c.RadialID(),c.AngularID());
    return cache->Has(IntegralsCache_Base::I3C::Overlap,abid,cid)
        ? cache->GetERI3() : cache->Set(MakeOverlap3C(c));
} 
//! 3 centre repulsion used for DFT \f$\left\langle a\left(1\right)b\left(1\right)\left|\frac{1}{r_{12}}\right|c\left(2\right)\right\rangle =\int d^{3}\vec{r}_{1}\:d^{3}\vec{r}_{2}\:g_{a}\left(\vec{r}_{1}\right)g_{b}\left(\vec{r}_{1}\right)\frac{1}{r_{12}}f_{c}\left(\vec{r}_{2}\right) \f$
template <class T> const ERI3<T>& Orbital_DFT_IBS<T>::Repulsion3C(const Fit_IBS& c) const
{
    auto cache=theGlobalCache;
    assert(cache);
    IntegralsCache_Base::IBS_ID_t abid(RadialID(),AngularID());
    IntegralsCache_Base::IBS_ID_t cid(c.RadialID(),c.AngularID());
    return cache->Has(IntegralsCache_Base::I3C::Repulsion,abid,cid)
        ? cache->GetERI3() : cache->Set(MakeRepulsion3C(c));
} 

template <class T> vec_t<T> Orbital_DFT_IBS<T>::Overlap3C(const smat_t<T>& Dcd, const Fit_IBS* c) const
{
    vec_t<T> ret(c->GetNumFunctions());
    auto& S=this->Overlap3C(*c);
    for(auto i:iv_t(0,S.size()))
        ret[i]=sum(Dcd%S[i]);
    return ret;
}

template <class T> vec_t<T> Orbital_DFT_IBS<T>::Repulsion3C(const smat_t<T>& Dcd, const Fit_IBS* c) const
{
    vec_t<T> ret(c->GetNumFunctions());
    auto& R=this->Repulsion3C(*c);
    for(auto i:iv_t(0,R.size()))
        ret[i]=sum(Dcd%R[i]);
    return ret;
}

template class Orbital_DFT_IBS<double>;
}
