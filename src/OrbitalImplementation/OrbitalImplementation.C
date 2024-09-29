// File: OrbitalImplementation.C  General implementation of an orbital, non-template part.



#include "OrbitalImplementation/OrbitalImplementation.H"
#include "oml/imp/binio.h"
#include <iostream>
#include <iomanip>

#define TYPE_STRING "BasisSet"
#define TYPE BasisSet
#include "Misc/Persistent/IDRef.Ci"

template class IDRef<const BasisSet>;

OrbitalImplementation::OrbitalImplementation(const IDRef<const BasisSet>& bs,double e, const Spin& S)
    : ElectronContainerImplementation(S,bs->GetQuantumNumber())
    , itsBasisSet   (bs)
    , itsEigenEnergy(e)
{};

OrbitalImplementation::OrbitalImplementation()
    : itsEigenEnergy(0)
{};

double OrbitalImplementation::GetEigenEnergy() const
{
    return itsEigenEnergy;
}

std::ostream& OrbitalImplementation::Write(std::ostream& os) const
{
    ElectronContainerImplementation::Write(os);
    if (Pretty())
        os << "              " << GetOccupation() << "/" << GetDegeneracy() << "       " << std::setw(12) << itsEigenEnergy << "      ";
    else
        os << itsBasisSet;
    
    if (Binary())
        BinaryWrite(itsEigenEnergy,os);
    if (Ascii ())
        os << itsEigenEnergy << " ";

    return os;
}

std::istream& OrbitalImplementation::Read (std::istream& is)
{
    ElectronContainerImplementation::Read(is);
    is >> itsBasisSet;
    if (Binary())
        BinaryRead(itsEigenEnergy,is);
    else
        is >> itsEigenEnergy >> std::ws;

    return is;
}

