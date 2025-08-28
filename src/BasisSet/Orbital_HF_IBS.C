// File: BasisSet/Orbital_HF_IBS.C  Interface for a Hartree-Fock (HF) Orbital Irrep Basis Set.
module;

export module qchem.Orbital_HF_IBS;
export import qchem.Orbital_1E_IBS;

import qchem.BasisSet.Internal.ERI4;
import Common.UniqueID;

export template <class T> class Orbital_HF_IBS;
//! \brief Interface for 4-center ERI integrals used in HF calculations.
//! This particular interface if for serving up ERIs between two Irrep Basis Sets (IRBs)
export template <class T> class Integrals_HF 
{
public:
    virtual ERI4 Direct  (const Orbital_HF_IBS<T>& c) const=0; //! <ab|1/r_12|cd>, this=a.
    virtual ERI4 Exchange(const Orbital_HF_IBS<T>& b) const=0; //! <ac|1/r_12|bd>, this=a.

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


export template <class T> class Orbital_HF_IBS
    : public virtual Orbital_IBS<T>
    , public virtual Integrals_HF<T> //Two electron integrals used for HF
{
public:
    using Integrals_HF<T>::Direct;
    using Integrals_HF<T>::Exchange;
    virtual SMatrix<T> Direct  (const SMatrix<T>& Dcd, const Orbital_HF_IBS<T>* bs_cd) const=0;
    virtual SMatrix<T> Exchange(const SMatrix<T>& Dcd, const Orbital_HF_IBS<T>* bs_cd) const=0;

};

