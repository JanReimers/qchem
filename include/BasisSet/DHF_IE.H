// File: DHF_IE.H Interface for integrals used in Dirac-Hartree-Fock (DHF) calculations.
#ifndef _DHF_IE_H_
#define _DHF_IE_H_

#include <BasisSet/IntegralEngine.H>


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

 #endif //_DHF_IE_H_
