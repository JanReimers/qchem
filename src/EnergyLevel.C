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

void EnergyLevel::SetOccupation(double ne)
{
    for (auto o:itsOrbitals) 
    {
        int degen=o->GetDegeneracy();
        if (degen>ne)
        {
            o->SetOccupation(ne);
            ne=0.0;
        }
        else
        {
            o->SetOccupation(degen);
            ne-=degen;
        }
        if (ne<=0.0) break;
    }
    //for (auto o:itsOrbitals) std::cout << *o << std::endl;
    
}

int  EnergyLevel::GetDegeneracy() const
{
    int  deg=0;
    int i=1;
    //std::cout << "  Energy level E=" << GetEnergy() << " with " << GetNumOrbitals() << " orbitals." << std::endl;
    for (auto o:itsOrbitals) 
    {
        //std::cout << "    Orbital " << i++ << " degen=" << o->GetDegeneracy() << " e=" << o->GetEigenEnergy() << std::endl;
        deg+=o->GetDegeneracy();
        //if (i>5) break;
    }
    //std::cout << "  ---------------------------------------" << std::endl;
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









