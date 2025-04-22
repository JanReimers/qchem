// File: Polarized_WF.C  Wave function for an unpolarized atom.



#include "Imp/WaveFunction/Polarized_WF.H"
#include "Imp/ChargeDensity/PolarizedCD.H"
#include "Imp/ChargeDensity/CompositeCD.H"
#include "Imp/WaveFunction/Irrep_WF.H"
#include "Imp/Orbitals/TOrbitals.H"
#include <BasisSet.H>
#include <Symmetry.H>
#include <Irrep_BS.H>
#include <cassert>
#include <Spin.H>
#include <Orbital_QNs.H>
#include <iostream>



Polarized_WF::Polarized_WF(const BasisSet* bs,const ElectronConfiguration* ec)
    : itsBS(bs) //Basis set
    , itsEC(ec) //Electron cofiguration
{
    assert(itsEC);
    assert(itsBS->GetNumFunctions()>0);
    for (auto b:itsBS->Iterate<TOrbital_IBS<double> >())
    {
        uiwf_t wfup(new Irrep_WF(b,Spin(Spin::Up)));
        uiwf_t wfdn(new Irrep_WF(b,Spin(Spin::Down)));
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
    auto i=itsQN_WFs.find(qns);
    assert(i!=itsQN_WFs.end());
    return i->second->GetOrbitals();
}


void Polarized_WF::FillOrbitals()
{
    itsUpELevels.clear();
    itsDnELevels.clear();

    for (auto& w:itsSpinUpIWFs) 
         itsUpELevels.merge(w->FillOrbitals(itsEC),0.0001);
    for (auto& w:itsSpinDnIWFs) 
         itsDnELevels.merge(w->FillOrbitals(itsEC),0.0001);
}


void Polarized_WF::DisplayEigen() const
{
    std::cout << "Spin:         up                 down            avg" << std::endl;
    auto iup=itsUpELevels.begin();
    auto idn=itsDnELevels.begin();
    while (iup!=itsUpELevels.end() && idn!=itsDnELevels.end())
    {
        bool valid_up = iup!=itsUpELevels.end();
        bool valid_dn = idn!=itsDnELevels.end();
        bool print_up = valid_up && iup->first<=0.0;
        bool print_dn = valid_dn && idn->first<=0.0;
        bool qn_match = valid_up && valid_dn && iup->second.qns.MatchNoSpin(idn->second.qns);
        if (print_up)
            iup->second.Report(std::cout);
        else
            std::cout << "                                     ";
        
        if (print_dn && qn_match)
            idn->second.Report(std::cout);
        
        
            
        std::cout << std::endl;
        if (valid_up) iup++;
        if (valid_dn && qn_match) idn++; //If no match, let iup catch up.
        if (iup->first>0.0 && idn->first>0.0) break;
    }
   
}

