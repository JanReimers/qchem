// File: Orbital.H  Interface for Orbital functions.
#ifndef _Orbital_H_
#define _Orbital_H_

#include "Orbital/ElectronContainer.H"
#include "Misc/pmstream.h"



//#############################################################
//
//  An orbital is 1) a real space function, 2) it holds some number
//  and type of electrons, 3) it has an eigen energy.  This non
//  templated portion is independant of whether the orbital real
//  or complex valued.
//
class Orbital
    : public virtual PMStreamableObject
    , public virtual ElectronContainer
{
public:
    virtual ~Orbital() {};
    virtual double   GetEigenEnergy(        ) const=0;

    virtual std::ostream& Write  (std::ostream&) const=0;
    virtual std::istream& Read   (std::istream&)      =0;
    static  Orbital* Factory     (std::istream&)        ;

private:
    friend class OrbitalGroup;
    virtual void FixUpPointer(const BasisSet*)=0;
};

#endif //_Orbital_H_
