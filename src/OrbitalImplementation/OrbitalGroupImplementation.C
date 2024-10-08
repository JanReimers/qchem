// File: OrbitalGroupImplementation.C  general orbital group implementation.


#include "OrbitalImplementation/OrbitalGroupImplementation.H"
#include "QuantumNumber.H"
#include "Misc/DFTDefines.H"
#include "Misc/ptr_vector1_io.h"
#include <cmath>
#include <iostream>

//-----------------------------------------------------------------
//
//  Construction zone
//
OrbitalGroupImplementation::OrbitalGroupImplementation()
    : itsBasisSet( )
    , itsRCBasisSet(0)
{};

OrbitalGroupImplementation::OrbitalGroupImplementation(const rc_ptr<const IrrepBasisSet>& bs)
    : itsBasisSet(&*bs)
    , itsRCBasisSet(bs)
{};

//-----------------------------------------------------------------
//
//  BasisSet reference injection for unpickling.
//
void OrbitalGroupImplementation::FixUpPointer(const rc_ptr<const IrrepBasisSet>& bs)
{
    itsBasisSet.FixUpPointer(&*bs);
    itsRCBasisSet=bs;
    OrbitalGroup::FixUpPointer(itsRCBasisSet);
}

//-----------------------------------------------------------------
//
//  Orbital group stuff.
//
index_t  OrbitalGroupImplementation::GetNumOrbitals() const
{
    return itsOrbitals.size();
}

double OrbitalGroupImplementation::GetEigenValueChange(const OrbitalGroup& og) const
{
    // No UT coverage
    // TODO: OrbitalGroup should return a vector of energies.
    double del=0;
    auto b2=og.begin();
    for (auto b1:*this)
    {
        del+=Square(b1->GetEigenEnergy()-(*b2)->GetEigenEnergy());
        b2++;
    }
    return sqrt(del);
}

//-----------------------------------------------------------------
//
//  Streamable stuff.
//
std::ostream& OrbitalGroupImplementation::Write(std::ostream& os) const
{
    if (!StreamableObject::Pretty())
    {
        os  << itsOrbitals;
        if (!StreamableObject::Binary()) os << std::endl;
    }
    else
    {
        os << "        Orbital group with " << GetNumOrbitals() << " " << itsBasisSet->GetQuantumNumber() << "orbitals:" << std::endl;
        os << "            Occupation      Energy      Eigenvector" << std::endl;
        os << itsOrbitals;
    }
    if (!StreamableObject::Pretty())
    {
        os  << itsBasisSet;
        if (StreamableObject::Ascii()) os << std::endl;
    }
    return os;
}

std::istream& OrbitalGroupImplementation::Read(std::istream& is)
{
    is >> itsOrbitals;
    if (!StreamableObject::Binary()) is >> std::ws;
    is  >> itsBasisSet;
    if (!StreamableObject::Binary()) is >> std::ws;
    return is;
}

