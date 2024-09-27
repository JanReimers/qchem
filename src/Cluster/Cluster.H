// File: Cluster.H  A cluster of atoms, Molecule, unit cell, solid.
#ifndef _Cluster_H_
#define _Cluster_H_

#include "Misc/Persistent/Persistent.H"
#include "Functions/ScalarFunction.H"
#include "Cluster/Atom.H"

class Cluster
    : public virtual Persistent
{
public:
    virtual ~Cluster() {};
    virtual bool operator==(const Cluster&) const
    {
        return false;
    }

    virtual void           Insert    (Atom*    )      =0;
    virtual int            GetNumAtoms        () const=0;
    virtual int            GetNuclearCharge   () const=0;
    virtual double         GetNetCharge       () const=0;
    virtual double         GetNumElectrons    () const=0;
    virtual ChargeDensity* GetChargeDensity   () const=0;

    static  Cluster*       Factory    (std::istream&)        ;

private:
};

#endif
