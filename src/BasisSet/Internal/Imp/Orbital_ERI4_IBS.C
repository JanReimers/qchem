// File: BasisSet/Internal/Imp/Orbital_ERI4_IBS.C  Build up 4C contractions for HF calculations.
module;
#include <cassert>
#include <stdexcept>
#include <string>   // std::string comparison for the canonical-order (BasisSetID) check
module qchem.BasisSet.Internal.Orbital_ERI4_IBS;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.DB_Cache;
import qchem.Blaze;

namespace qchem::BasisSet
{

// The one place the contraction face is narrowed to the ERI4 substrate.  A null result means a caller
// paired an ERI4 basis with a basis that has no ERI4 blocks -- the SALC decorator is the only such
// basis today, and it overrides the whole contraction face, so the pairing can only arise by mistake.
// Loud, and naming both sides, because the pre-R1.7 behaviour here was an empty ERI4: a zero Fock
// contribution with no diagnostic at all.
template <class T> const Orbital_ERI4_IBS<T>& Orbital_ERI4_IBS<T>::Substrate(const Orbital_HF_IBS<T>* cd) const
{
    assert(cd);
    const Orbital_ERI4_IBS<T>* e=dynamic_cast<const Orbital_ERI4_IBS<T>*>(cd);
    if (!e)
        // Name(), NOT BasisSetID(): the ID is the cache key, and a molecular one is a whole
        // screenful (and currently carries hex addresses -- see R1.9).  Name() is the human face.
        throw std::runtime_error("Orbital_ERI4_IBS: the HF partner basis '"+cd->Name()+
                                 "' has no ERI4 substrate, so no (ab|cd) block spans it and '"+
                                 this->Name()+"'.  An ERI4 basis cannot be contracted against a "
                                 "basis that builds its Fock block some other way (e.g. a SALC decorator, "
                                 "which slices the whole-AO Fock instead).");
    return *e;
}

template <class T> const ERI4& Orbital_ERI4_IBS<T>::Direct(const Orbital_ERI4_IBS<T>& c) const
{
    return theCache<T>().Get(IntegralsCache_Base::I4C::Direct,this,&c,
        [this,&c]{ return MakeDirect(c); });
}

template <class T> const   ERI4&  Orbital_ERI4_IBS<T>::Exchange(const Orbital_ERI4_IBS<T>& c) const
{
    return theCache<T>().Get(IntegralsCache_Base::I4C::Exchange,this,&c,
        [this,&c]{ return MakeExchange(c); });
}

template <class T> void Orbital_ERI4_IBS<T>::AccumulateDirect(rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* cd) const
{
    assert(!blazem::isnan(Dcd));
    assert(blazem::max(blazem::abs(Dcd))>0.0);  //Dcd should be pre-screened for zero.
    const Orbital_ERI4_IBS<T>* ab=this;
    MatMul(Sab,ab->Direct(Substrate(cd)),Dcd); //ERI4 Jabcd=ab->Direct(*cd);
}

template <class T> void Orbital_ERI4_IBS<T>::AccumulateExchange(rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* cd) const
{
    assert(!blazem::isnan(Dcd));
    assert(blazem::max(blazem::abs(Dcd))>0.0);  //Dcd should be pre-screened for zero.
    const Orbital_ERI4_IBS<T>* ab=this;
    MatMul(Sab,ab->Exchange(Substrate(cd)),Dcd); // ERI4 Kabcd=ab->Exchange(*cd);
}

// Fetch the ONE canonical Direct block J(this,cd) and scatter it into both irreps' Fock blocks -- the
// bra-ket partner J(cd,this) is never requested, so it is never built or cached (doc/ERI4Rework.md §4).
// Di/Dj may be zero (an empty irrep); ScatterBoth adds an exact zero for that side, so no pre-screen here
// (the caller short-circuits only the both-empty case).
// Request the CANONICAL (min-BasisSetID) block only -- the ERI4 cache throws on the partner orientation.
// ScatterBoth(Si,Sj,Di,Dj) on a block with Nab=Si-side, Ncd=Sj-side gives Si+=J·Dj and Sj+=Jᵀ·Di; when the
// partner cd sorts first we fetch its block J(cd,this) and swap the two target/density pairs, which yields
// the identical result (J(cd,this)=J(this,cd)^T).  The irrep pair here is always off-diagonal (this!=cd).
template <class T> void Orbital_ERI4_IBS<T>::AccumulateDirectBoth(rsmat_t& Ji, rsmat_t& Jj, const smat_t<T>& Di, const smat_t<T>& Dj, const Orbital_HF_IBS<T>* cd) const
{
    assert(!blazem::isnan(Di) && !blazem::isnan(Dj));
    const Orbital_ERI4_IBS<T>& e=Substrate(cd);
    if (this->BasisSetID() <= e.BasisSetID()) Direct(e).ScatterBoth(Ji,Jj,Di,Dj);        // canonical (this,cd)
    else                                      e.Direct(*this).ScatterBoth(Jj,Ji,Dj,Di);  // canonical (cd,this): swap
}

// Exchange counterpart, same canonical-only fetch + swap.
template <class T> void Orbital_ERI4_IBS<T>::AccumulateExchangeBoth(rsmat_t& Ki, rsmat_t& Kj, const smat_t<T>& Di, const smat_t<T>& Dj, const Orbital_HF_IBS<T>* cd) const
{
    assert(!blazem::isnan(Di) && !blazem::isnan(Dj));
    const Orbital_ERI4_IBS<T>& e=Substrate(cd);
    if (this->BasisSetID() <= e.BasisSetID()) Exchange(e).ScatterBoth(Ki,Kj,Di,Dj);        // canonical (this,cd)
    else                                      e.Exchange(*this).ScatterBoth(Kj,Ki,Dj,Di);  // canonical (cd,this): swap
}

template class Orbital_ERI4_IBS<double>;

} //namespace
