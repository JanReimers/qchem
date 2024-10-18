// File: OrbitalImplementation.C  General implementation of an orbital, non-template part.



#include "Imp/Orbitals//Orbital.H"
#include "oml/imp/binio.h"
#include <iostream>
#include <iomanip>

#define TYPE_STRING "IrrepBasisSet"
#define TYPE IrrepBasisSet
#include "Misc/Persistent/IDRef.Ci"

template class IDRef<const IrrepBasisSet>;

OrbitalImp::OrbitalImp(const IDRef<const IrrepBasisSet>& bs,double e, const Spin& S)
    : ElectronContainerImp(S,bs->GetQuantumNumber())
    , itsBasisSet   (bs)
    , itsEigenEnergy(e)
{};

OrbitalImp::OrbitalImp()
    : itsEigenEnergy(0)
{};

double OrbitalImp::GetEigenEnergy() const
{
    return itsEigenEnergy;
}

std::ostream& OrbitalImp::Write(std::ostream& os) const
{
    ElectronContainerImp::Write(os);
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

std::istream& OrbitalImp::Read (std::istream& is)
{
    ElectronContainerImp::Read(is);
    is >> itsBasisSet;
    if (Binary())
        BinaryRead(itsEigenEnergy,is);
    else
        is >> itsEigenEnergy >> std::ws;

    return is;
}

