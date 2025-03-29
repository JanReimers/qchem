// File: TBasisSetImplementation.C


#include "Imp/BasisSet/TIBS_Common.H"
#include "Imp/BasisSet/DFT_IBS_Common.H"
#include <Fit_IBS.H>
#include "oml/vector.h"

template <class T> typename Orbital_DFT_IBS_Common<T>::Vec Orbital_DFT_IBS_Common<T>::
Overlap3C(const SMat& Dcd, const fbs_t* ff) const
{
    Vec ret(ff->size());
    const typename Integrals_Base<T>::ERI3& S=this->Overlap3C(*ff);
    for(auto i:ret.indices())
        ret(i)=Dot(Dcd,S[i-1]);
    return ret;
}

template <class T> typename Orbital_DFT_IBS_Common<T>::Vec Orbital_DFT_IBS_Common<T>::
Repulsion3C(const SMat& Dcd, const fbs_t* ff) const
{
    Vec ret(ff->size());
    const typename Integrals_Base<T>::ERI3& repulsion=this->Repulsion3C(*ff);
    for(auto i:ret.indices())
        ret(i)=Dot(Dcd,repulsion[i-1]);
    return ret;
}

template class Orbital_DFT_IBS_Common<double>;

