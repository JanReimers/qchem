// File: EnergyLevel.C  Energy level with degeneracy and orbital list.



#include <EnergyLevel.H>
#include <Orbital.H>
#include <cmath>
#include <iostream>
#include <cassert>

EnergyLevel::EnergyLevel()
    : itsEnergy(0)
    , itsTolerance(0.000001)
    , itsOrbitals()
{};

EnergyLevel::EnergyLevel(double tol)
    : itsEnergy(0)
    , itsTolerance(tol)
    , itsOrbitals()
{};

bool EnergyLevel::IsOccupied() const
{
    bool occ=false;
    for (auto o:itsOrbitals) occ = occ || o->IsOccupied();
    return occ;
}

double EnergyLevel::GetOccupation() const
{
    double occ=0;
    for (auto o:itsOrbitals)  occ += o->GetOccupation();
    return occ;
}

double EnergyLevel::GetSpin() const
{
    double spin=0;
    for (auto o:itsOrbitals)  spin += o->GetSpin();
    return spin;
}

void EnergyLevel::Empty()
{
    for (auto o:itsOrbitals) o->Empty();
}

void EnergyLevel::SetOccupation(double OccupationFactor)
{
    assert(OccupationFactor>=0);
    assert(OccupationFactor<=1);
    double NumElectronsPerOrbital=GetDegeneracy()*OccupationFactor/GetNumOrbitals();
    for (auto o:itsOrbitals) o->SetOccupation(NumElectronsPerOrbital);
}

int  EnergyLevel::GetDegeneracy() const
{
    int  deg=0;
    for (auto o:itsOrbitals) deg+=o->GetDegeneracy();
    return deg;
};

index_t EnergyLevel::GetNumOrbitals () const
{
    return itsOrbitals.size();
}

double EnergyLevel::GetEnergy() const
{
    return itsEnergy;
}

bool EnergyLevel::Add(Orbital* orbital)
{
    bool empty=(GetNumOrbitals()==0);
    if (empty) itsEnergy=orbital->GetEigenEnergy();

    bool matchEnergy = (fabs(itsEnergy - orbital->GetEigenEnergy()) <= itsTolerance);
    if (matchEnergy) itsOrbitals.push_back(orbital);

    return matchEnergy;
}


std::ostream& EnergyLevel::Write(std::ostream& os) const
{
    return os << itsEnergy << " " << itsTolerance << " "
           << GetDegeneracy() << " " << GetNumOrbitals() << std::endl;
}

std::istream& EnergyLevel::Read(std::istream& is)
{
    return is >> itsEnergy >> itsTolerance;
}









