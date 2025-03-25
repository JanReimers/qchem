// File: DHF_IBS.H  Interface for a Dirac-Hartree-Fock (HF) Orbital Irrep Basis Set.
#ifndef _DHF_IBS_H_
#define _DHF_IBS_H_

#include <Irrep_BS.H>
#include <DHF_IE.H>

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
    //int GetKappa() const;
};

#endif //_DHF_IBS_H_
