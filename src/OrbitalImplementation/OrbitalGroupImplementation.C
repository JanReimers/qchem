// File: OrbitalGroupImplementation.C  general orbital group implementation.


#include "OrbitalImplementation/OrbitalGroupImplementation.H"
#include "Orbital/OrbitalGroupBrowser.H"
#include "Misc/DFTDefines.H"
#include "Misc/ptrvector_io.h"
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

OrbitalGroupImplementation::OrbitalGroupImplementation(const rc_ptr<const BasisSet>& bs)
    : itsBasisSet(&*bs)
    , itsRCBasisSet(bs)
{};

//-----------------------------------------------------------------
//
//  BasisSet reference injection for unpickling.
//
void OrbitalGroupImplementation::FixUpPointer(const rc_ptr<const BasisSet>& bs)
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
    double del=0;
    OrbitalGroupBrowser b1(*this);
    OrbitalGroupBrowser b2(og);
    for(; b1&&b2; b1++,b2++) del+=Square((*b1).GetEigenEnergy()-(*b2).GetEigenEnergy());
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
        os << "Orbital group with " << GetNumOrbitals() << " orbitals" << std::endl;
        os << "Occupation      Energy      Eigenvector" << std::endl;
        os << itsOrbitals << std::endl;
    }
    if (!StreamableObject::Pretty())
    {
        os  << itsBasisSet;
        if (!StreamableObject::Binary()) os << std::endl;
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

