// File: BasisSet/Imp/Orbital_HF_IBS1.C  Build up 4C contractions for HF calculations.
module;
#include <cassert>
module qchem.BasisSet.Orbital_HF_IBS;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.DB_Cache;
import qchem.Blaze;

namespace qchem::BasisSet
{

template <class T> const ERI4& Orbital_HF_IBS<T>::Direct(const Orbital_HF_IBS<T>& c) const
{
    return theCache<T>().Get(IntegralsCache_Base::I4C::Direct,this,&c,
        [this,&c]{ return MakeDirect(c); });
}

template <class T> const   ERI4&  Orbital_HF_IBS<T>::Exchange(const Orbital_HF_IBS<T>& c) const
{
    return theCache<T>().Get(IntegralsCache_Base::I4C::Exchange,this,&c,
        [this,&c]{ return MakeExchange(c); });
}

template <class T> void Orbital_HF_IBS<T>::AccumulateDirect(rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* cd) const
{
    assert(!blazem::isnan(Dcd));
    assert(blazem::max(blazem::abs(Dcd))>0.0);  //Dcd should be pre-screened for zero.
    const Orbital_HF_IBS<T>* ab=this;
    MatMul(Sab,ab->Direct(*cd),Dcd); //ERI4 Jabcd=ab->Direct(*cd);
}

template <class T> void Orbital_HF_IBS<T>::AccumulateExchange(rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* cd) const
{
    assert(!blazem::isnan(Dcd));
    assert(blazem::max(blazem::abs(Dcd))>0.0);  //Dcd should be pre-screened for zero.
    const Orbital_HF_IBS<T>* ab=this;
    MatMul(Sab,ab->Exchange(*cd),Dcd); // ERI4 Kabcd=ab->Exchange(*cd);
}

template class Orbital_HF_IBS<double>;

} //namespace