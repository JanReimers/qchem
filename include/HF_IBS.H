// File: HF_IBS.H  Interface for a Hartree-Fock (HF) Orbital Irrep Basis Set.
#ifndef _HF_IBS_H_
#define _HF_IBS_H_

#include <Irrep_BS.H>
#include <HF_IE.H>

template <class T> class TOrbital_HF_IBS
    : public virtual TOrbital_IBS<T>
    , public virtual Integrals_HF<T> //Two electron integrals used for HF
{
    typedef typename Integrals_Base<T>::SMat SMat;
    typedef typename Integrals_Base<T>::obs_t obs_t;
public:
    using Integrals_HF<T>::Direct;
    using Integrals_HF<T>::Exchange;
    virtual SMat Direct  (const SMat& Dcd, const obs_t* bs_cd) const=0;
    virtual SMat Exchange(const SMat& Dcd, const obs_t* bs_cd) const=0;

};

#endif //_HF_IBS_H_
