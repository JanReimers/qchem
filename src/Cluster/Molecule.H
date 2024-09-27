// File: Molecule.H  Implementation for A cluster of atoms.
#ifndef _Molecule_H_
#define _Molecule_H_

#include "Cluster/Cluster.H"
#include "Misc/ptr_vector.h"

class Molecule :
    public virtual Cluster
{
public:
    Molecule();

    virtual void           Insert        (Atom*)      ;
    virtual int            GetNumAtoms        () const;
    virtual int            GetNuclearCharge   () const;
    virtual double         GetNetCharge       () const;
    virtual double         GetNumElectrons    () const;
    virtual ChargeDensity* GetChargeDensity   () const;

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

private:
    friend class ClusterBrowser;

    double      itsNumElectrons;
    optr_vector<Atom*> itsAtoms;
};

#endif

