//  file: ElectronDumper.cpp  Takes a list of orbitals, sorts them in terms of energy,
//                      makes a list of energy levels, and finally fills the
//                      levels with electrons.

#include "Imp/WaveFunction/ElectronDumper.H"
#include <Orbital.H>
#include <EnergyLevel.H>
#include "Imp/WaveFunction/FermiThermalizer.H"
#include <algorithm> //sort
#include <numeric> //iota
#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>

ElectronDumper::ElectronDumper(double tol, double kT)
    : IsDirty        (true)
    , itsOrbitals    (    )
    , itsEnergyLevels(    )
    , itsTolerance   (tol )
    , itskT          (kT  )
    , itsFermiEnergy (0   )
{};

void ElectronDumper::Add(Orbital* o)
{
    itsOrbitals.push_back(o);
    IsDirty=true;
}

void ElectronDumper::Add(Orbitals* og)
{
    for (auto o:*og) itsOrbitals.push_back(o);
    IsDirty=true;
}

template <class T> void ReIndex(std::vector<T*>& a, std::vector<int> index)
{
    assert(a.size()==index.size());
    std::vector<T*> dest;
    for (auto b:index)  dest.push_back(a[b]);
    a=dest;
}

void ElectronDumper::MakeEnergyLevels()
{
    //
    //  Sort orbitals by energy.
    //
    std::vector<double> Energies;
    for(auto b:itsOrbitals) Energies.push_back(b->GetEigenEnergy());
    
    std::vector<int> indexs(itsOrbitals.size());
    std::iota(indexs.begin(),indexs.end(),0); //Initialize indexs=[0,1,2,3...]
    std::sort( indexs.begin(),indexs.end(), [&](int i,int j){return Energies[i]<Energies[j];} );
    ReIndex(itsOrbitals,indexs);

    //
    //  Update all the EnergyLevel objects.
    //
    itsEnergyLevels.clear();
    itsEnergyLevels.push_back(new EnergyLevel(itsTolerance));
    {
        for(auto i:itsOrbitals)
        {
            bool OrbitalWasAdded = itsEnergyLevels.back()->Add(i); //Try to add orbital to energy level.
            if(!OrbitalWasAdded)
            {
                itsEnergyLevels.push_back(new EnergyLevel(itsTolerance));
                itsEnergyLevels.back()->Add(i);
            }
        }
    }
//    std::cout << "Sorted energy levels:" << std::endl;
//    for (auto e : itsEnergyLevels)
//        std::cout << "  Level i E=" << e->GetEnergy() << " with " << e->GetNumOrbitals() << " orbitals." << std:: endl;
//    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    IsDirty=false;
}

void ElectronDumper::DumpInElectrons(double NumElectrons)
{
    if (IsDirty) MakeEnergyLevels();
    FermiThermalizer ft(itsEnergyLevels,itskT,NumElectrons);
    itsFermiEnergy=ft.GetFermiEnergy();

    double Ntotal=NumElectrons;
    for (auto i:itsEnergyLevels)
    {
        double n=ft.GetOccupation(i->GetEnergy());
        double degen=i->GetDegeneracy();
//        std:: cout << "1 Ntotal,n,degen:  " << Ntotal << " " << n << " " << degen << std::endl; 
        if (degen>Ntotal)
            n=Ntotal;  //finish off the remaining electrons.
        else
            n=degen;  //Fill the orbital(s)
        i->SetOccupation(n); //Dump electron into orbitals at this E level.
        Ntotal-=n;
//        std::cout << "2 Ntotal,n,degen:  " << Ntotal << " " << n << " " << degen << std::endl; 
//        std::cout << "===================================================" << std::endl; 
        if (Ntotal<=0.0) break;
    }

}

double ElectronDumper::GetFermiEnergy() const
{
    return itsFermiEnergy;
}

std::ostream& operator<<(std::ostream& os,const ElectronDumper& ed)
{
    int nLumo=2; //Show two LUMO ;eve;
    for (auto b:ed.itsEnergyLevels)
    {
        os.setf(std::ios::fixed,std::ios::floatfield);
        os << std::setw(8) << std::setprecision(4) << b->GetEnergy();
        os << "(";
//        os.setf(std::ios::scientific,std::ios::floatfield);
        os << std::setw(4) << std::setprecision(1);
        if (b->GetOccupation() != b->GetDegeneracy()) os << b->GetOccupation() << "/";
        os << b->GetDegeneracy() << ") ";
        os << b->GetQuantumNumbers() << " ";
        os << std::endl;
        if (b->GetOccupation() ==0 && --nLumo==0) break;
    }
    return os;
}


