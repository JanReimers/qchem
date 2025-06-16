// File: DFT_IBS.H  Interface for a Density Functional Theory (DFT) Orbital Irrep Basis Set.
#ifndef _DFT_IBS_H_
#define _DFT_IBS_H_

#include <BasisSet/Irrep_BS.H>
#include <BasisSet/DFT_IE.H>
template <class T> class TOrbital_DFT_IBS
    : public virtual TOrbital_IBS<T>
    , public virtual Integrals_DFT<T> //DFT integrals
    
{
    typedef typename Integrals_Base<T>::SMat SMat;
    typedef typename Integrals_Base<T>::fbs_t fbs_t;
    typedef Vector<T> Vec;
public:
    virtual fbs_t*    CreateCDFitBasisSet(const BasisSet*,const Cluster*) const=0;
    virtual fbs_t*    CreateVxcFitBasisSet(const BasisSet*,const Cluster*) const=0;
    using Integrals_DFT<T>::Overlap3C; //Unhide
    using Integrals_DFT<T>::Repulsion3C; //Unhide
    virtual Vec Overlap3C  (const SMat& Dcd, const fbs_t* ff) const=0;
    virtual Vec Repulsion3C(const SMat& Dcd, const fbs_t* ff) const=0;
};

#endif //_DFT_IBS_H_
