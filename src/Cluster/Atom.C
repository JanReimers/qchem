// File: Atom.C  A Atom like data type.



#include "Cluster/Atom.H"
#include "AtomMesh.H"
#include <Mesh/MeshParams.H>
#include "oml/io3d.h"
#include <iostream>
#include <cassert>
#include "Mesh/MHLRadialMesh.H"
#include "Mesh/GaussAngularMesh.H"


Atom::Atom()
    : itsZ(0)
    , itsCharge(0)
    , itsR( )
{};

Atom::Atom(int Z, double charge)
    : itsZ(Z)
    , itsCharge(charge)
    , itsR(0,0,0)
{
    assert(itsZ>0);
    assert(itsZ<150); //Maybe there is an island of stability at Z=140!!!!
};

Atom::Atom(int Z, double charge, const RVec3& R)
    : itsZ(Z)
    , itsCharge(charge)
    , itsR(R)
{
    assert(itsZ>0);
    assert(itsZ<150); //Maybe there is an island of stability at Z=140!!!!
};

Mesh* Atom::CreateMesh(const MeshParams& mp) const
{
    MHLRadialMesh    rm(mp.Nradial,mp.MHL_m,mp.MHL_alpha);
    GaussAngularMesh am(mp.Nangle);  
    return new AtomMesh(rm,am,itsR); 
}


double Atom::GetNumElectrons() const
{
    return itsZ-itsCharge;
}

std::ostream& Atom::Write  (std::ostream& os) const
{
    if (Binary())
    {
        BinaryWrite(itsZ,os);
        os << itsR << " ";
    }
    if (Ascii())
    {
        os << itsZ << " " << itsR  << " ";
    }
    if (Pretty())
    {
        os.setf(std::ios::fixed,std::ios::floatfield);
        os << std::setw(4) << itsZ << "    "
        << std::setw(5) << std::setprecision(2) << itsR << "     ";
    }
    if (!Binary()) os << std::endl;
    return os;
}

std::istream& Atom::Read   (std::istream& is)
{
    if (StreamableObject::Binary())
    {
        BinaryRead(itsZ,is);
        is >> itsR >> std::ws;
    }
    else
    {
        is >> itsZ >> itsR >>  std::ws;
    }

    return is;
}

