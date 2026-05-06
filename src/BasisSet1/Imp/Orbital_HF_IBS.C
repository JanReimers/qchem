// File: BasisSet/Imp/Orbital_HF_IBS1.C  Build up 4C contractions for HF calculations.
module;
#include <cassert>
#include <blaze/math/SymmetricMatrix.h>
module qchem.BasisSet1.Orbital_HF_IBS;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet1.DB_Cache;

namespace BasisSet1
{

template <class T> const ERI4& Orbital_HF_IBS<T>::Direct(const Orbital_HF_IBS<T>& c) const
{
    auto cache=theGlobalCache;
    assert(cache);
    IntegralsCache_Base::IBS_ID_t ab(RadialID(),AngularID());
    IntegralsCache_Base::IBS_ID_t cd(c.RadialID(),c.AngularID());
    return cache->Has(IntegralsCache_Base::I4C::Direct,ab,cd)
        ? cache->GetERI4() : cache->SetDirect(MakeDirect(c));
}

template <class T> const   ERI4&  Orbital_HF_IBS<T>::Exchange(const Orbital_HF_IBS<T>& c) const
{
    auto cache=theGlobalCache;
    assert(cache);
    IntegralsCache_Base::IBS_ID_t ab(RadialID(),AngularID());
    IntegralsCache_Base::IBS_ID_t cd(c.RadialID(),c.AngularID());
    return cache->Has(IntegralsCache_Base::I4C::Exchange,ab,cd)
        ? cache->GetERI4() : cache->SetExchange(MakeExchange(c));
}

template <class T> void Orbital_HF_IBS<T>::AccumulateDirect(rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* cd) const
{
    assert(!isnan(Dcd));
    assert(max(abs(Dcd))>0.0);  //Dcd should be pre-screened for zero.
    const Orbital_HF_IBS<T>* ab=this;
    MatMul(Sab,ab->Direct(*cd),Dcd); //ERI4 Jabcd=ab->Direct(*cd);
}

template <class T> void Orbital_HF_IBS<T>::AccumulateExchange(rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* cd) const
{
    assert(!isnan(Dcd));
    assert(max(abs(Dcd))>0.0);  //Dcd should be pre-screened for zero.
    const Orbital_HF_IBS<T>* ab=this;
    MatMul(Sab,ab->Exchange(*cd),Dcd); // ERI4 Kabcd=ab->Exchange(*cd);
}

template class Orbital_HF_IBS<double>;

} //namespace