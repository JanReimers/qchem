// File: EnergyLevel.C  Energy level with degeneracy and orbital list.



#include "Orbital/EnergyLevel.H"
#include "Orbital/Orbital.H"
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
    ptr_vector<Orbital*>::const_iterator b(itsOrbitals.begin());
    for (; b!=itsOrbitals.end(); b++) occ = occ || b->IsOccupied();
    return occ;
}

double EnergyLevel::GetOccupation() const
{
    double occ=0;
    ptr_vector<Orbital*>::const_iterator b(itsOrbitals.begin());
    for (; b!=itsOrbitals.end(); b++)  occ += b->GetOccupation();
    return occ;
}

double EnergyLevel::GetSpin() const
{
    double spin=0;
    ptr_vector<Orbital*>::const_iterator b(itsOrbitals.begin());
    for (; b!=itsOrbitals.end(); b++)  spin += b->GetSpin();
    return spin;
}

void EnergyLevel::Empty()
{
    for (ptr_vector<Orbital*>::iterator i(itsOrbitals.begin()); i!=itsOrbitals.end(); i++) i->Empty();
}

void EnergyLevel::SetOccupation(double OccupationFactor)
{
    assert(OccupationFactor>=0);
    assert(OccupationFactor<=1);
    double NumElectronsPerOrbital=GetDegeneracy()*OccupationFactor/GetNumOrbitals();
    for (ptr_vector<Orbital*>::iterator i(itsOrbitals.begin()); i!=itsOrbitals.end(); i++) i->SetOccupation(NumElectronsPerOrbital);
}

int  EnergyLevel::GetDegeneracy() const
{
    int  deg=0;
    ptr_vector<Orbital*>::const_iterator b(itsOrbitals.begin());
    for (; b!=itsOrbitals.end(); b++) deg+=b->GetDegeneracy();
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









