// File: ElectronContainerImplementation.C  General implementation Electron container.



#include "Imp/Orbitals/ElectronContainer.H"
#include <QuantumNumber.H>
#include "oml/imp/binio.h"
#include <iostream>
#include <cassert>
#include <stdlib.h>

ElectronContainerImp::ElectronContainerImp(const Spin& S,const QuantumNumber& qn)
    : itsOccupation       (0)
    , itsSpin             (S)
    , itsOrbitalDegeneracy(qn.GetDegeneracy())
{};

ElectronContainerImp::ElectronContainerImp()
    : itsOccupation       (0)
    , itsSpin             ( )
    , itsOrbitalDegeneracy(1)
{};

double ElectronContainerImp::GetSpin() const
{
    double ret=0;
    switch (itsSpin.itsState)
    {
    case Spin::Down :
        ret=-0.5;
        break;
    case Spin::None :
        ret= 0.0;
        break;
    case Spin::Up   :
        ret= 0.5;
        break;
    };
    return ret;
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

void  ElectronContainerImp::SetOccupation(double n)
{
    assert(n>=0);
    if (n>GetDegeneracy())
    {
        std::cerr << "Trying to load " << n << " electrons into an orbital with degeneracy "
                  << GetDegeneracy() << std::endl;
//        exit(-1);
    }
    itsOccupation=n;
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
