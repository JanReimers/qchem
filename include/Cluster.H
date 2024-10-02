// File: Cluster.H  A cluster of atoms, Molecule, unit cell, solid.
#ifndef _Cluster_H_
#define _Cluster_H_

#include "Misc/Persistent/Persistent.H"
#include "Functions/ScalarFunction.H"
#include "Cluster/Atom.H"
#include <vector>

class Cluster
    : public virtual Persistent
{
public:
    typedef std::vector<Atom*>::      iterator       iterator;
    typedef std::vector<Atom*>::const_iterator const_iterator;

    virtual ~Cluster() {};
    virtual bool operator==(const Cluster&) const
    {
        return false;
    }

    virtual void           Insert    (Atom*    )      =0;
    virtual size_t         GetNumAtoms        () const=0;
    virtual int            GetNuclearCharge   () const=0;
    virtual double         GetNetCharge       () const=0;
    virtual double         GetNumElectrons    () const=0;
    virtual ChargeDensity* GetChargeDensity   () const=0;
    
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
    virtual       iterator begin()      =0;
    virtual       iterator end  ()      =0;

    static  Cluster*       Factory    (std::istream&)        ;

private:
};

#endif
