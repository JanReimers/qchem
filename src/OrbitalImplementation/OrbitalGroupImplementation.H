// File: OrbitalGroupImplementation.H  general orbital group implementation.
#ifndef _OrbitalGroupImplementation_H_
#define _OrbitalGroupImplementation_H_



#include "Orbital.H"
#include "BasisSet.H"
#include "Imp/Containers/ptr_vector.h"

class OrbitalGroupImplementation
    : public virtual OrbitalGroup
{
public:
    OrbitalGroupImplementation();
    OrbitalGroupImplementation(const IrrepBasisSet*);
    virtual ~OrbitalGroupImplementation() {};

    virtual index_t  GetNumOrbitals     (                   ) const;
    virtual double   GetEigenValueChange(const OrbitalGroup&) const;

    virtual const_iterator begin() const {return itsOrbitals.begin();}
    virtual const_iterator end  () const {return itsOrbitals.end  ();} 
    virtual       iterator begin()       {return itsOrbitals.begin();}
    virtual       iterator end  ()       {return itsOrbitals.end  ();} 

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

protected: //Why not private!!!!!??
//    virtual void  FixUpPointer(const rc_ptr<const IrrepBasisSet>&);

    optr_vector1<Orbital*> itsOrbitals;
    const IrrepBasisSet* itsBasisSet;
};

#endif //_OrbitalGroupImplementation_H_
