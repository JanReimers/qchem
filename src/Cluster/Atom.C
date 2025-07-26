// File: Atom.C  A Atom like data type.
module;
#include <iosfwd>

export module qchem.Atom;
export import qchem.Types;
import Common.UniqueIDImp;
import qchem.Mesh;
import qchem.Streamable;

export class Atom
    : public virtual Streamable
    , public UniqueIDImp
{
public:
    Atom();
    Atom(int Z, double charge);
    Atom(int Z, double charge, const RVec3& R);

    virtual double GetNumElectrons      () const;
    virtual Mesh*  CreateMesh(const MeshParams&) const;

    virtual std::ostream& Write  (std::ostream&) const;

    int     itsZ;      //Atomic number.
    double  itsCharge; //Net charge. Z-numElectrons.
    RVec3   itsR;      //Spatial position.

private:
    Atom& operator=(const Atom&); //Why?  Rule of 5/6
};




