// File: Atom.C  A Atom like data type.



#include "Imp/Cluster/Atom.H"
#include "Imp/Cluster/AtomMesh.H"
#include <Mesh.H>
#include <ChargeDensity.H>
#include "oml/imp/binio.h"
#include "oml/io3d.h"
//#include "Misc/Unpickle.H"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <stdlib.h>
#include "Imp/Mesh/MHLRadialMesh.H"
#include "Imp/Mesh/GaussAngularMesh.H"


Atom::Atom()
    : itsZ(0)
    , itsCharge(0)
    , itsR( )
{};

Atom::Atom(int Z, double charge, const RVec3& R)
    : itsZ(Z)
    , itsCharge(charge)
    , itsR(R)
{
    assert(itsZ>0);
    assert(itsZ<150); //Maybe there is an island of stability at Z=140!!!!
};

Mesh* Atom::Create_MHL_G_Mesh(size_t Nradial, size_t Nangle) const
{
    MHLRadialMesh    rm(Nradial,3U,2.0);
    GaussAngularMesh am(Nangle);  
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

