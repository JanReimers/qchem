// File: ElectronContainerImplementation.C  General implementation Electron container.



#include "Imp/Orbitals/ElectronContainer.H"
#include <Symmetry.H>
#include "oml/imp/binio.h"
#include <iostream>
#include <cassert>
#include <stdlib.h>

ElectronContainerImp::ElectronContainerImp(const Spin& S,const Symmetry& qn)
    : itsOccupation       (0)
    , itsSpin             (S)
    , itsOrbitalDegeneracy(qn.GetDegeneracy())
{
//    std::cout << "itsOrbitalDegeneracy=" << itsOrbitalDegeneracy << std::endl;
};

ElectronContainerImp::ElectronContainerImp()
    : itsOccupation       (0)
    , itsSpin             ( )
    , itsOrbitalDegeneracy(1)
{};

Spin ElectronContainerImp::GetSpin() const
{
    return itsSpin;
}

int ElectronContainerImp::GetDegeneracy() const
{
    return itsSpin.GetDegeneracy() * itsOrbitalDegeneracy;
}

double ElectronContainerImp::GetOccupation() const
{
    return itsOccupation;
}

bool ElectronContainerImp::IsOccupied() const
{
    return itsOccupation>0;
}

void ElectronContainerImp::Empty()
{
    itsOccupation=0;
}

double ElectronContainerImp::TakeElectrons(double n)
{
    assert(n>=0);
    double g=GetDegeneracy();
    itsOccupation=std::min(g,n);
    return n-itsOccupation;
}


std::ostream& ElectronContainerImp::Write(std::ostream& os) const
{
    if (StreamableObject::Binary())
    {
        BinaryWrite(itsOccupation   ,os);
        BinaryWrite(itsOrbitalDegeneracy            ,os);
    }
    if (StreamableObject::Ascii())
        os << itsOccupation << " "  << itsOrbitalDegeneracy << " ";
    if (!StreamableObject::Pretty()) os << itsSpin;
    return os;
}

std::istream& ElectronContainerImp::Read (std::istream& is)
{
    if (StreamableObject::Binary())
    {
        BinaryRead(itsOccupation,is);
        BinaryRead(itsOrbitalDegeneracy         ,is);
    }
    else
    {
        is >> itsOccupation >> itsOrbitalDegeneracy;
    }
    return is >> itsSpin;
}
