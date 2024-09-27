// File: Atom.H  A Atom like data type.
#ifndef _Atom_H_
#define _Atom_H_

#include "types.H"
#include "Misc/pmstream.h"
#include "oml/vector3d.h"
#include <string>

class Atom
    : public virtual PMStreamableObject
{
public:
    Atom();
    Atom(int Z, double charge, const RVec3& R);

    void  SetMesh             (Mesh*);
    void  SetMeshFile         (const char*);
    void  SetChargeDensityFile(const char*);

    Mesh*          GetIntegrationMesh() const;
    ChargeDensity* GetChargeDensity  () const;
    double         GetNumElectrons   () const;

    virtual std::ostream& Write  (std::ostream&) const;
    virtual std::istream& Read   (std::istream&)      ;
    static  Atom*    Factory(std::istream&);

    int     itsZ;      //Atomic number.
    double  itsCharge; //Net charge. Z-numElectrons.
    RVec3   itsR;      //Spatial position.

private:
    Atom& operator=(const Atom&);

    std::string  itsMeshFileName; //Integration Mesh.
    std::string  itsChargeDensityFileName;
    mutable Mesh* itsMesh;
};

#endif
