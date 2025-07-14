// File: Atom.H  A Atom like data type.
#ifndef _Atom_H_
#define _Atom_H_

#include <omlfwd.H>
#include <Mesh/fwd.H>
#include "Common/pmstream.h"
#include "oml/vector3d.h"
#include <string>

import Common.UniqueIDImp;

class Atom
    : public virtual PMStreamableObject
    , public UniqueIDImp
{
public:
    Atom();
    Atom(int Z, double charge);
    Atom(int Z, double charge, const RVec3& R);

    virtual double GetNumElectrons      () const;
    virtual Mesh*  CreateMesh(const MeshParams&) const;

    virtual std::ostream& Write  (std::ostream&) const;
    static  Atom*    Factory(std::istream&);

    int     itsZ;      //Atomic number.
    double  itsCharge; //Net charge. Z-numElectrons.
    RVec3   itsR;      //Spatial position.

private:
    Atom& operator=(const Atom&);
};

#endif
