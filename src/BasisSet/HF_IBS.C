// File: HF_IBS.H  Interface for a Hartree-Fock (HF) Orbital Irrep Basis Set.
module;
#include "BasisSet/fwd.H"

export module qchem.HF_IBS;
import qchem.BasisSet.Integrals;
import qchem.BasisSet.ERI4;

import Common.UniqueID;
export import qchem.Irrep_BS;

//! \brief Interface for 4-center ERI integrals used in HF calculations.
//! This particular interface if for serving up ERIs between two Irrep Basis Sets (IRBs)
export template <class T> class Integrals_HF : public virtual Integrals_Base<T>
{
public:
    typedef TOrbital_IBS<T> obs_t;
    virtual ERI4 Direct  (const obs_t& c) const=0; //! <ab|1/r_12|cd>, this=a.
    virtual ERI4 Exchange(const obs_t& b) const=0; //! <ac|1/r_12|bd>, this=a.

};

//! \brief Interface for 4-center ERI integrals used in HF calculations.
//! This particular interface if for serving up ERIs between two Irrep Basis Sets (IRBs)
//! using IBS IDs.
export template <class T> class Integrals_BS_2E
{
public:
    typedef UniqueID::IDtype IDType;   
    virtual ERI4 Direct  (IDType a,IDType c) const=0;
    virtual ERI4 Exchange(IDType a,IDType b) const=0;
};


export template <class T> class TOrbital_HF_IBS
    : public virtual TOrbital_IBS<T>
    , public virtual Integrals_HF<T> //Two electron integrals used for HF
{
protected:
    typedef typename Integrals_Base<T>::SMat SMat;
    typedef TOrbital_HF_IBS obs_t;
public:
    using Integrals_HF<T>::Direct;
    using Integrals_HF<T>::Exchange;
    virtual SMat Direct  (const SMat& Dcd, const obs_t* bs_cd) const=0;
    virtual SMat Exchange(const SMat& Dcd, const obs_t* bs_cd) const=0;

};

