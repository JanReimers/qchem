// File: Molecule.C  Implementation for A cluster of atoms.



#include "Imp/Cluster/Molecule.H"
#include "Imp/Cluster/Atom.H"
//#include "Imp/ChargeDensity/CompositeCD.H"
#include "Imp/Cluster/MoleculeMesh.H"
#include "Imp/Containers/ptr_vector_io.h"
#include "oml/imp/binio.h"
#include <iostream>
#include <iomanip>
#include <cassert>

Molecule::Molecule()
    : itsNumElectrons(0)
    , itsAtoms    ( )
{};

void Molecule::Insert(Atom* a)
{
    itsAtoms.push_back(a);
    itsNumElectrons+=a->GetNumElectrons();
}

size_t Molecule::GetNumAtoms() const
{
    return itsAtoms.size();
}

int Molecule::GetNuclearCharge() const
{
    int chg=0;
    for(auto b:*this) chg+=b->itsZ;
    return chg;
}

double Molecule::GetNetCharge() const
{
    return GetNuclearCharge()-itsNumElectrons;
}

double Molecule::GetNumElectrons() const
{
    return itsNumElectrons;
}

Mesh*  Molecule::CreateMesh(const MeshParams& mp) const
{
    return new MoleculeMesh(*this,mp);
}

//ChargeDensity* Molecule::GetChargeDensity() const
//{
//    CompositeCD* ret=new CompositeCD;
//    for(auto b:*this) ret->Insert(b->GetChargeDensity());
//    return ret;
//}

std::ostream& Molecule::Write(std::ostream& os) const
{
    if (!StreamableObject::Pretty())
    {
        UniqueID::Write(os);
        if (StreamableObject::Binary())
        {
            BinaryWrite(GetNumElectrons(),os);
        }
        else
        {
            os << GetNumElectrons() << " ";
        }
        os << itsAtoms;
    }
    else
    {
        os << "Molecule with " << GetNumAtoms() << " atoms"
        << ", nuclear charge " << GetNuclearCharge() << "(e)"
        << ", net charge "<< GetNetCharge() << "(e)" << std::endl;
        os << "Atom #  Element  Position vector     Mesh file    Charge density file" << std::endl;
        int i=1;
        for (auto b:*this) os << std::setw(5) << i++ << "   " << *b;
        os << std::endl;
    }

    return os;
}

std::istream& Molecule::Read(std::istream& is)
{
    UniqueID::Read(is);
    if (StreamableObject::Binary())
    {
        BinaryRead(itsNumElectrons,is);
    }
    else
    {
        is >> itsNumElectrons >> std::ws;
    }
    is >> itsAtoms;
    return is;
}


