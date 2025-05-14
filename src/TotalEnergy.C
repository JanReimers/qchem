// File: TotalEnergy.C  Store and display a breakdown of the total energy.



#include "TotalEnergy.H"
#include <iostream>

EnergyBreakdown::EnergyBreakdown()
    : Kinetic  (0)
    , Enn      (0)
    , Een      (0)
    , Eee      (0)
    , EeeFit   (0)
    , EeeFitFit(0)
    , Exc      (0)
    , ExcFit   (0)
    , ExcFitFit(0)
{};

double EnergyBreakdown::GetVirial() const
{
    double V=GetPotentialEnergy();
    return V/Kinetic;
}

EnergyBreakdown& EnergyBreakdown::operator+=(const EnergyBreakdown& e1)
{
    Kinetic    += e1.Kinetic;
    Enn        += e1.Enn;
    Een        += e1.Een;
    Eee        += e1.Eee;
    EeeFit     += e1.EeeFit;
    EeeFitFit  += e1.EeeFitFit;
    Exc       += e1.Exc;
    ExcFit    += e1.ExcFit;
    ExcFitFit += e1.ExcFitFit;
    return *this;
}

using std::cout;
using std::endl;

void EnergyBreakdown::Display() const
{
    cout << endl;
    cout << "Total energy breakdown :" << endl;
    cout << "------------------------" << endl;
    cout << "Kinetic   :" << Kinetic << endl;
    cout << "Enn       :" << Enn << endl;
    cout << "Een       :" << Een << endl;
    cout << "Eee       :" << Eee << endl;
    cout << "EeeFit    :" << EeeFit << endl;
    cout << "EeeFitFit :" << EeeFitFit << endl;
    cout << "Exc       :" << Exc << endl;
    cout << "ExcFit    :" << ExcFit << endl;
    cout << "ExcFitFit :" << ExcFitFit << endl;
    cout << "------------------------" << endl << endl;
}
