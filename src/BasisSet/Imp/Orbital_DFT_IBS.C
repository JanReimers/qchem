// File: BasisSet/Imp/Orbital_DFT_IBS.C
module;
#include <cassert>
module qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.Internal.DB_Cache;
import qchem.Blaze;

namespace BasisSet
{
template <class T> const ERI3<T>& Orbital_DFT_IBS<T>::Overlap3C  (const Fit_IBS& c) const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I3C::Overlap,this,&c,
        [this,&c]{ return MakeOverlap3C(c); });
} 
//! 3 centre repulsion used for DFT \f$\left\langle a\left(1\right)b\left(1\right)\left|\frac{1}{r_{12}}\right|c\left(2\right)\right\rangle =\int d^{3}\vec{r}_{1}\:d^{3}\vec{r}_{2}\:g_{a}\left(\vec{r}_{1}\right)g_{b}\left(\vec{r}_{1}\right)\frac{1}{r_{12}}f_{c}\left(\vec{r}_{2}\right) \f$
template <class T> const ERI3<T>& Orbital_DFT_IBS<T>::Repulsion3C(const Fit_IBS& c) const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I3C::Repulsion,this,&c,
        [this,&c]{ return MakeRepulsion3C(c); });
} 

template <class T> vec_t<T> Orbital_DFT_IBS<T>::Overlap3C(const smat_t<T>& Dcd, const Fit_IBS* c) const
{
    vec_t<T> ret(c->GetNumFunctions());
    auto& S=this->Overlap3C(*c);
    for(auto i:iv_t(0,S.size()))
        ret[i]=blazem::sum(Dcd%S[i]);
    return ret;
}

template <class T> vec_t<T> Orbital_DFT_IBS<T>::Repulsion3C(const smat_t<T>& Dcd, const Fit_IBS* c) const
{
    vec_t<T> ret(c->GetNumFunctions());
    auto& R=this->Repulsion3C(*c);
    for(auto i:iv_t(0,R.size()))
        ret[i]=blazem::sum(Dcd%R[i]);
    return ret;
}

template class Orbital_DFT_IBS<double>;
}
