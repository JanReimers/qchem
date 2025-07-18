// File: Atom.C  A Atom like data type.
module;
#include <iostream>
#include <Common/pmstream.h>

export module qchem.Atom;
import Common.UniqueIDImp;
import Mesh;
import oml;

export class Atom
    : public virtual PMStreamableObject
    , public UniqueIDImp
{
public:
    using RVec3=Vector3D<double>;
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




