// File: OrbitalIMp.C  Non template part of hte Orbital implemtaton.



#include <WaveFunction/EnergyLevel.H>
#include <Symmetry/Symmetry.H>
#include "TOrbital.H"
// #include <iostream>
// #include <cassert>
// #include <stdlib.h>

OrbitalImp::OrbitalImp(double e, const Orbital_QNs& qns)
    : itsEigenEnergy(e)
    , itsOccupation       (0)
    , itsQNs(qns)
{};

OrbitalImp::OrbitalImp()
    : itsOccupation       (0)
{};

double OrbitalImp::GetEigenEnergy() const
{
    return itsEigenEnergy;
}

Orbital_QNs OrbitalImp::GetQNs() const
{
    return itsQNs;
}

std::string OrbitalImp::GetLabel() const
{
    std::ostringstream oss;
    oss << GetQNs();
    std::string s(oss.str());
    return s;
}

int OrbitalImp::GetDegeneracy() const
{
    return itsQNs.GetDegeneracy();
}

double OrbitalImp::GetOccupation() const
{
    return itsOccupation;
}

bool OrbitalImp::IsOccupied() const
{
    return itsOccupation>0;
}

void OrbitalImp::Empty()
{
    itsOccupation=0;
}

double OrbitalImp::TakeElectrons(double n)
{
    assert(n>=0);
    double g=GetDegeneracy();
    itsOccupation=std::min(g,n);
    return n-itsOccupation;
}


std::ostream& OrbitalImp::Write(std::ostream& os) const
{
    return os << itsQNs;
}

