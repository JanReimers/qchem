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

// Fetch the ONE canonical Direct block J(this,cd) and scatter it into both irreps' Fock blocks -- the
// bra-ket partner J(cd,this) is never requested, so it is never built or cached (doc/ERI4Rework.md §4).
// Di/Dj may be zero (an empty irrep); ScatterBoth adds an exact zero for that side, so no pre-screen here
// (the caller short-circuits only the both-empty case).
template <class T> void Orbital_HF_IBS<T>::AccumulateDirectBoth(rsmat_t& Ji, rsmat_t& Jj, const smat_t<T>& Di, const smat_t<T>& Dj, const Orbital_HF_IBS<T>* cd) const
{
    assert(!blazem::isnan(Di) && !blazem::isnan(Dj));
    Direct(*cd).ScatterBoth(Ji,Jj,Di,Dj);
}

// Exchange counterpart: fetch the ONE canonical Exchange block K(this,cd) and scatter into both irreps'
// K blocks; K(cd,this) is never requested/built (doc/ERI4Rework.md §4).
template <class T> void Orbital_HF_IBS<T>::AccumulateExchangeBoth(rsmat_t& Ki, rsmat_t& Kj, const smat_t<T>& Di, const smat_t<T>& Dj, const Orbital_HF_IBS<T>* cd) const
{
    assert(!blazem::isnan(Di) && !blazem::isnan(Dj));
    Exchange(*cd).ScatterBoth(Ki,Kj,Di,Dj);
}

template class Orbital_HF_IBS<double>;

} //namespace