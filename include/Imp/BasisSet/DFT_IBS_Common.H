// File: DFT_IBS_Common.H  Common implementation for all DFT Irrep Basis Sets.
#ifndef _DFT_IBS_Common_H_
#define _DFT_IBS_Common_H_

#include <BasisSet/DFT_IBS.H>
template <class T> class Orbital_DFT_IBS_Common
    : public virtual TOrbital_DFT_IBS<T>
{
    typedef typename Integrals_Base<T>::SMat SMat;
    typedef typename Integrals_DFT<T>::Vec  Vec;
public:
    using Integrals_DFT<T>::Overlap3C; //Unhide
    using Integrals_DFT<T>::Repulsion3C; //Unhide
    virtual Vec Overlap3C  (const SMat& Dcd, const Fit_IBS* ff) const;
    virtual Vec Repulsion3C(const SMat& Dcd, const Fit_IBS* ff) const;
};


#endif //_DFT_IBS_Common_H_