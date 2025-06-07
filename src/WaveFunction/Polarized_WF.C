// File: Polarized_WF.C  Wave function for an unpolarized atom.

#include "Imp/WaveFunction/Polarized_WF.H"
#include "Imp/ChargeDensity/PolarizedCD.H"
#include <cassert>
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;


Polarized_WF::Polarized_WF(const BasisSet* bs,const ElectronConfiguration* ec,SCFAccelerator& acc)
    : Composite_WF(bs,ec) 
{
    MakeIrrep_WFs(acc,Spin::Up);
    MakeIrrep_WFs(acc,Spin::Down);
};

DM_CD* Polarized_WF::GetChargeDensity() const
{
    return new Polarized_CDImp(GetChargeDensity(Spin::Up),GetChargeDensity(Spin::Down));
}

WaveFunction::sf_t* Polarized_WF::GetSpinDensity() const
{
    return new SpinDensity(GetChargeDensity(Spin::Up),GetChargeDensity(Spin::Down));
}


void Polarized_WF::DisplayEigen() const
{
    StreamableObject::SetToPretty();
    cout << "                         |       Spin up           |          Spin Down      |    Spin " << endl;
    cout << " Orbital                 |  Occ/g  |     E         |  Occ/g  |       E       |  Splitting" << endl;
    EnergyLevels els_up=GetEnergyLevels(Spin::Up), els_dn=GetEnergyLevels(Spin::Down);
    for (auto iup:els_up)
    {
        auto up=iup.second;
        assert(iup.first == up.e);
        Orbital_QNs dnqns(up.qns.n,Spin::Down,up.qns.sym);
        auto dn=els_dn.find(dnqns); 

        cout  << up.qns << " (" << std::fixed << std::setw(2) << std::setprecision(0) << up.occ  << "/"  << std::setw(2) << up.degen << ") |";
        cout  << std::fixed << std::setw(14) << std::setprecision(8) << up.e << " |";
        if ((dn.e<=0.0 && dn.occ>0) || up.e>0.0) 
        {
            cout  << " (" << std::fixed << std::setw(2) << std::setprecision(0) << dn.occ  << "/"  << std::setw(2) << dn.degen << ") |";
            cout << std::setw(14)  << std::setprecision(8) << dn.e << " |";
            if (fabs(up.e-dn.e)>1e-4) 
            {
                cout.precision(4);
                cout << std::fixed << std::setw(8) << up.e-dn.e;
            }
        }
        else
        {
            cout << "         |               |";
        }
        
        
        cout << endl;
        if (up.e>0.0) break;
    }
   
}

