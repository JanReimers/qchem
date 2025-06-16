// File: DHF_IBS.H  Interface for a Dirac-Hartree-Fock (HF) Orbital Irrep Basis Set.
#ifndef _DHF_IBS_H_
#define _DHF_IBS_H_

#include <BasisSet/Irrep_BS.H>
#include <BasisSet/DHF_IE.H>

template <class T> class Orbital_RKB_IBS
    : public virtual TIrrepBasisSet<T>
    , public virtual Integrals_RKB<T> 
{
public:
    
};

template <class T> class Orbital_RKBL_IBS
    : public virtual TIrrepBasisSet<T>
    , public virtual Integrals_RKBL<T> //One electron integrals used for everything
{
public:
    //int GetKappa() const;
};

template <class T> class Orbital_RKBS_IBS
    : public virtual TIrrepBasisSet<T>
    , public virtual Integrals_RKBS<T> //One electron integrals used for everything
{
public:
    // using IrrepBasisSet::size;
    //int GetKappa() const;
    virtual void InsertBasisFunctions(const Orbital_RKBL_IBS<T>* l)=0;
};

#endif //_DHF_IBS_H_
