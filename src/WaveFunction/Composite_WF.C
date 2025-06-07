// File: Composite_WF.H  Wave function as a list of Irrep wave functions.

#include "Imp/WaveFunction/Composite_WF.H"
#include "Imp/WaveFunction/Irrep_WF.H"
#include "Imp/ChargeDensity/CompositeCD.H"
#include "Imp/Orbitals/TOrbitals.H"
#include "Imp/SCFAccelerator.H"
#include <BasisSet.H>
#include <Irrep_BS.H>
#include <cassert>
#include <Spin.H>
#include <Orbital_QNs.H>
#include <iostream>

Composite_WF::Composite_WF()
    : itsEC(0)
{};

Composite_WF::Composite_WF(const BasisSet* bs,const ElectronConfiguration* ec)
    : itsBS(bs)
    , itsEC(ec)
{
    assert(itsBS);
    assert(itsEC);
    assert(itsBS->GetNumFunctions()>0);
    
};

void Composite_WF::MakeIrrep_WFs(SCFAccelerator& acc, Spin s)
{

    for (auto b:itsBS->Iterate<TOrbital_IBS<double> >())
    {
        uiwf_t wf(new Irrep_WF(b,s,acc.Create(b)));
        itsQN_WFs[wf->GetQNs()]=wf.get();
        itsSpin_WFs[s].push_back(wf.get());
        itsIWFs.push_back(std::move(wf)); //Do the move last.
    }
}

//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the FillOrbitals member function
//  to fill up the orbitals with electrons.
//
void Composite_WF::DoSCFIteration(Hamiltonian& ham,const DM_CD* cd)
{
    for (auto& w:itsIWFs) w->DoSCFIteration(ham,cd);
}

DM_CD* Composite_WF::GetChargeDensity(Spin s) const
{
    auto i = itsSpin_WFs.find(s);
    assert(i!=itsSpin_WFs.end());
    Composite_CD* cd = new Composite_CD();
    for (auto& w:i->second) cd->Insert(w->GetChargeDensity());
    return cd;
}

const Orbitals* Composite_WF::GetOrbitals(const Irrep_QNs& qns) const
{
    return const_cast<Composite_WF*>(this)->GetOrbitals(qns);
}
Orbitals* Composite_WF::GetOrbitals(const Irrep_QNs& qns) 
{
    assert(qns.ms==Spin::None);
    auto i=itsQN_WFs.find(qns);
    assert(i!=itsQN_WFs.end());
    return i->second->GetOrbitals();

}

Composite_WF::iqns_t Composite_WF::GetQNs() const
{
    iqns_t iqns;
    for (auto q:itsQN_WFs) iqns.push_back(q.first);
    return iqns;
}



void Composite_WF::FillOrbitals()
{
    itsELevels.clear();
    for (auto& w:itsIWFs) 
        itsELevels.merge(w->FillOrbitals(itsEC),0.0001);
}




