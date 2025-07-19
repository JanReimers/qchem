// File: DHF_IBS.H  Interface for a Dirac-Hartree-Fock (HF) Orbital Irrep Basis Set.
#ifndef _DHF_IBS_H_
#define _DHF_IBS_H_

#include <BasisSet/Integrals.H>


//! \brief Interface for one electron integrals used in Dirac-Hartree-Fock (DHF) calculations.
template <class T> class Integrals_RKB
: public virtual Integrals_Overlap<T>
, public virtual Integrals_Kinetic<T>
, public virtual Integrals_Nuclear<T>
, public virtual Integrals_RestMass<T>
{
public:
    
};

//! \brief Interface for Large-Component one electron integrals used in Dirac-Hartree-Fock (DHF) calculations.
template <class T> class Integrals_RKBL
: public virtual Integrals_Overlap<T>
, public virtual Integrals_XKinetic<T>
, public virtual Integrals_Nuclear<T>
{
public:
    
};

//! \brief Interface for Small-Component one electron integrals used in Dirac-Hartree-Fock (DHF) calculations.
template <class T> class Integrals_RKBS 
: public virtual Integrals_Kinetic<T> //Serves as the overlap.
, public virtual Integrals_Nuclear<T>
{
public:
   
 };

 import qchem.Irrep_BS;

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
