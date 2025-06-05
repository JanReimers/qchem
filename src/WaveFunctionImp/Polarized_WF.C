// File: Polarized_WF.C  Wave function for an unpolarized atom.



#include "Imp/WaveFunction/Polarized_WF.H"
#include "Imp/ChargeDensity/PolarizedCD.H"
#include "Imp/ChargeDensity/CompositeCD.H"
#include "Imp/WaveFunction/Irrep_WF.H"
#include "Imp/Orbitals/TOrbitals.H"
#include "Imp/SCFAccelerator.H"
#include <BasisSet.H>
#include <Symmetry.H>
#include <Irrep_BS.H>
#include <cassert>
#include <Spin.H>
#include <Orbital_QNs.H>
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;


Polarized_WF::Polarized_WF(const BasisSet* bs,const ElectronConfiguration* ec,const SCFAccelerator& acc)
    : itsBS(bs) //Basis set
    , itsEC(ec) //Electron cofiguration
{
    assert(itsEC);
    assert(itsBS->GetNumFunctions()>0);
    for (auto b:itsBS->Iterate<TOrbital_IBS<double> >())
    {
        uiwf_t wfup(new Irrep_WF(b,Spin(Spin::Up  ),acc.Create(b)));
        uiwf_t wfdn(new Irrep_WF(b,Spin(Spin::Down),acc.Create(b)));
        itsQN_WFs[wfup->GetQNs()]=wfup.get();
        itsQN_WFs[wfdn->GetQNs()]=wfdn.get();
        // Do tte move last.
        itsSpinUpIWFs.push_back(std::move(wfup));
        itsSpinDnIWFs.push_back(std::move(wfdn));
    }
};

//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the FillOrbitals member function
//  to fill up the orbitals with electrons.
//
void Polarized_WF::DoSCFIteration(Hamiltonian& ham,const DM_CD* cd)
{
    for (auto& w:itsSpinUpIWFs) w->DoSCFIteration(ham,cd);
    for (auto& w:itsSpinDnIWFs) w->DoSCFIteration(ham,cd);
}

DM_CD* Polarized_WF::GetChargeDensity() const
{
    Composite_CD* up = new Composite_CD();
    Composite_CD* dn = new Composite_CD();
    for (auto& w:itsSpinUpIWFs) up->Insert(w->GetChargeDensity());
    for (auto& w:itsSpinDnIWFs) dn->Insert(w->GetChargeDensity());
    return new Polarized_CDImp(up,dn);
}

WaveFunction::sf_t* Polarized_WF::GetSpinDensity() const
{
    Composite_CD* up = new Composite_CD();
    Composite_CD* dn = new Composite_CD();
    for (auto& w:itsSpinUpIWFs) up->Insert(w->GetChargeDensity());
    for (auto& w:itsSpinDnIWFs) dn->Insert(w->GetChargeDensity());
    return new SpinDensity(up,dn);
}


const Orbitals* Polarized_WF::GetOrbitals(const Irrep_QNs& qns) const
{
    return const_cast<Polarized_WF*>(this)->GetOrbitals(qns);
}

Orbitals* Polarized_WF::GetOrbitals(const Irrep_QNs& qns) 
{
    auto i=itsQN_WFs.find(qns);
    assert(i!=itsQN_WFs.end());
    return i->second->GetOrbitals();
}

EnergyLevels  Polarized_WF::GetEnergyLevels () const
{
    EnergyLevels ret(itsUpELevels);
    ret.merge(itsDnELevels);
    return ret;
} 

Polarized_WF::iqns_t Polarized_WF::GetQNs() const
{
    iqns_t iqns;
    for (auto q:itsQN_WFs) iqns.push_back(q.first);
    return iqns;
}


void Polarized_WF::FillOrbitals()
{
    itsUpELevels.clear();
    itsDnELevels.clear();

    for (auto& w:itsSpinUpIWFs) 
         itsUpELevels.merge(w->FillOrbitals(itsEC));//,0.000001);
    for (auto& w:itsSpinDnIWFs) 
         itsDnELevels.merge(w->FillOrbitals(itsEC));//,0.000001);

    // cout << "FillOrbitals Up:" << endl;
    // itsUpELevels.Report(cout);
    // cout << "FillOrbitals Down:" << endl;
    // itsDnELevels.Report(cout);

}


void Polarized_WF::DisplayEigen() const
{
    StreamableObject::SetToPretty();
    cout << "                         |       Spin up           |          Spin Down      |    Spin " << endl;
    cout << " Orbital                 |  Occ/g  |     E         |  Occ/g  |       E       |  Splitting" << endl;
    for (auto iup:itsUpELevels)
    {
        auto up=iup.second;
        assert(iup.first == up.e);
        Orbital_QNs dnqns(up.qns.n,Spin::Down,up.qns.sym->Clone());
        auto dn=itsDnELevels.find(dnqns); 

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

