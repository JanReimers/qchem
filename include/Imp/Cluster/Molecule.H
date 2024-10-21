// File: Molecule.H  Implementation for A cluster of atoms.
#ifndef _Molecule_H_
#define _Molecule_H_

#include "Cluster.H"
#include "Imp/Containers/ptr_vector.h"
#include <iostream>

class Molecule :
    public virtual Cluster
{
public:
    Molecule();

    virtual void   Insert        (Atom*)      ;
    virtual size_t GetNumAtoms        () const;
    virtual int    GetNuclearCharge   () const;
    virtual double GetNetCharge       () const;
    virtual double GetNumElectrons    () const;
    virtual Mesh*  CreateMesh(const MeshParams&) const;
    
    virtual const_iterator begin() const {return itsAtoms.begin();}
    virtual const_iterator end  () const {return itsAtoms.end  ();} 
    virtual       iterator begin()       {return itsAtoms.begin();}
    virtual       iterator end  ()       {return itsAtoms.end  ();} 

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

private:

    double       itsNumElectrons;
    optr_vector1<Atom*> itsAtoms;
};

#endif

