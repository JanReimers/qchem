// File: HF_IBS_Common.H  Common implementation for all Hartree-Fock (HF) Irrep Basis Sets.
#ifndef _HF_IBS_Common_H_
#define _HF_IBS_Common_H_

#include <BasisSet/HF_IBS.H>
template <class T> class Orbital_HF_IBS_Common
    : public virtual TOrbital_HF_IBS<T>
{
    typedef typename Integrals_Base<T>::SMat SMat;
    typedef typename TOrbital_HF_IBS<T>::obs_t obs_t;
public:
    virtual SMat Direct  (const SMat& Dcd, const obs_t* bs_cd) const;
    virtual SMat Exchange(const SMat& Dcd, const obs_t* bs_cd) const;
};


#endif //_HF_IBS_Common_H_