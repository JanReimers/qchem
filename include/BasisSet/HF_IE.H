// File: HF_IE.H Integrals used for Hartree-Fock calculations.
#ifndef _HF_IE_H_
#define _HF_IE_H_

#include <BasisSet/IntegralEngine.H>

//! \brief Interface for 4-center ERI integrals used in HF calculations.
//! This particular interface if for serving up ERIs between two Irrep Basis Sets (IRBs)
template <class T> class Integrals_HF : public virtual Integrals_Base<T>
{
    typedef Integrals_Base<T> Base;
    typedef typename Base::obs_t obs_t;
public:
    virtual ERI4 Direct  (const obs_t& c) const=0; //! <ab|1/r_12|cd>, this=a.
    virtual ERI4 Exchange(const obs_t& b) const=0; //! <ac|1/r_12|bd>, this=a.

};

//! \brief Interface for 4-center ERI integrals used in HF calculations.
//! This particular interface if for serving up ERIs between two Irrep Basis Sets (IRBs)
//! using IBS IDs.
template <class T> class Integrals_BS_2E
{
public:
    typedef UniqueID::IDtype IDType;   
    virtual ERI4 Direct  (IDType a,IDType c) const=0;
    virtual ERI4 Exchange(IDType a,IDType b) const=0;
};

#endif //_HF_IE_H_
