// File: Molecule.C  Implementation for A cluster of atoms.
module;
#include <iostream>

export module qchem.Molecule;

import qchem.Atom;
import qchem.Cluster;
import Common.UniqueIDImp;
import Mesh;

export class Molecule 
    : public virtual Cluster
    , public UniqueIDImp
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

private:
    double       itsNumElectrons;
    av_t         itsAtoms;
};



