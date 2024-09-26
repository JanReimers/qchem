// File: Molecule.C  Implementation for A cluster of atoms.



#include "Cluster/Molecule.H"
#include "Cluster/ClusterBrowser.H"
#include "ChargeDensityImplementation/CompositeCD/CompositeCD.H"
#include "Misc/ptrvector_io.h"
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

int Molecule::GetNumAtoms() const
{
    return itsAtoms.size();
}

int Molecule::GetNuclearCharge() const
{
    int chg=0;
    for(ClusterBrowser b(*this); b; b++) chg+=b->itsZ;
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

ChargeDensity* Molecule::GetChargeDensity() const
{
    CompositeCD* ret=new CompositeCD;
    for(ClusterBrowser b(*this); b; b++) ret->Insert(b->GetChargeDensity());
    return ret;
}

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
        ClusterBrowser b(*this);
        for(int i=1; b; b++,i++) os << std::setw(5) << i << "   " << *b;
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


