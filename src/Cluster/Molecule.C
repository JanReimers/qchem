// File: Molecule.C  Implementation for A cluster of atoms.



#include <Common/pmstream.h>
#include "Common/stl_io.h"
#include <iostream>
#include <iomanip>
#include <cassert>

#include <Cluster/Molecule.H>
#include "MoleculeMesh.H"
import qchem.Atom;

Molecule::Molecule()
    : itsNumElectrons(0)
    , itsAtoms    ( )
{};

void Molecule::Insert(Atom* a)
{
    itsAtoms.push_back(std::unique_ptr<Atom>(a));
    itsNumElectrons+=a->GetNumElectrons();
}

size_t Molecule::GetNumAtoms() const
{
    return itsAtoms.size();
}

int Molecule::GetNuclearCharge() const
{
    int chg=0;
    for(auto& b:*this) chg+=b->itsZ;
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

std::ostream& Molecule::Write(std::ostream& os) const
{
    os << "Molecule with " << GetNumAtoms() << " atoms"
    << ", nuclear charge " << GetNuclearCharge() << "(e)"
    << ", net charge "<< GetNetCharge() << "(e)" << std::endl;
    os << "Atom #  Element  Position vector     Mesh file    Charge density file" << std::endl;
    int i=1;
    for (auto& b:*this) os << std::setw(5) << i++ << "   " << *b;
    os << std::endl;

    return os;
}



