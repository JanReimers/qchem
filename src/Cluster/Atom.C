// File: Atom.C  A Atom like data type.



#include "Cluster/Atom.H"
#include "Mesh/Mesh.H"
#include "ChargeDensity.H"
#include "oml/imp/binio.h"
#include "oml/io3d.h"
#include "Misc/Unpickle.H"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <stdlib.h>


Atom::Atom()
    : itsZ(0)
    , itsCharge(0)
    , itsR( )
    , itsMeshFileName         ("NULL") //This helps with IO.
    , itsChargeDensityFileName("NULL") //This helps with IO.
    , itsMesh(0)
{};

Atom::Atom(int Z, double charge, const RVec3& R)
    : itsZ(Z)
    , itsCharge(charge)
    , itsR(R)
    , itsMeshFileName         ("NULL") //This helps with IO.
    , itsChargeDensityFileName("NULL") //This helps with IO.
    , itsMesh(0)
{
    assert(itsZ>0);
    assert(itsZ<150); //Maybe there is an island of stability at Z=140!!!!
};

void  Atom::SetMesh(Mesh* m)
{
    itsMesh=m->Clone();
    itsMesh->ShiftOrigin(itsR);
}
void Atom::SetMeshFile(const char* filename)
{
    if(!std::ifstream(filename))
    {
        std::cerr << "Atom::SetMeshFile Can't open mesh file '" << filename << "'" << std::endl;
        exit(-1);
    }
    itsMeshFileName=filename;
}

void Atom::SetChargeDensityFile(const char*  filename)
{
    if(!std::ifstream(filename))
    {
        std::cerr << "Atom::SetChargeDensityFile Can't open charge density file '" << filename << "'" << std::endl;
        exit(-1);
    }
    itsChargeDensityFileName=filename;
}

Mesh* Atom::GetIntegrationMesh() const
{
    if (itsMesh) return itsMesh;
    if(!UnPickle(itsMesh,itsMeshFileName.c_str(),"integration mesh"))
    {
        std::cerr << "Atom::GetIntegrationMesh mesh filename not initialized" << std::endl;
        exit(-1);
    }
    itsMesh->ShiftOrigin(itsR);
    return itsMesh;
}

ChargeDensity* Atom::GetChargeDensity() const
{
    ChargeDensity* ret=0;
    if(!UnPickle(ret,itsChargeDensityFileName.c_str(),"charge density"))
    {
        std::cerr << "Atom::GetChargeDensity charge density filename not initialized" << std::endl;
        exit(-1);
    }
    ret->ShiftOrigin(itsR);
    return ret;
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
        os << itsR << " " << itsMeshFileName << " " << itsChargeDensityFileName << " ";
    }
    if (Ascii())
    {
        os << itsZ << " " << itsR << " " << itsMeshFileName << " " << itsChargeDensityFileName << " ";
    }
    if (Pretty())
    {
        os.setf(std::ios::fixed,std::ios::floatfield);
        os << std::setw(4) << itsZ << "    "
        << std::setw(5) << std::setprecision(2) << itsR << "     "
        << itsMeshFileName << "     "
        << itsChargeDensityFileName << "  ";
    }
    if (!Binary()) os << std::endl;
    return os;
}

std::istream& Atom::Read   (std::istream& is)
{
    if (StreamableObject::Binary())
    {
        BinaryRead(itsZ,is);
        is >> itsR >> itsMeshFileName >> std::ws >> itsChargeDensityFileName >> std::ws;
    }
    else
    {
        is >> itsZ >> itsR >> itsMeshFileName >> std::ws >> itsChargeDensityFileName >>  std::ws;
    }

    return is;
}

