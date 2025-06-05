// File: UnPolarized_WF.C  Wave function for an unpolarized atom.

#include "Imp/WaveFunction/UnPolarized_WF.H"
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

UnPolarized_WF::UnPolarized_WF()
    : itsEC(0)
{};

UnPolarized_WF::UnPolarized_WF(const BasisSet* bs,const ElectronConfiguration* ec,SCFAccelerator& acc)
    : itsBS(bs)
    , itsEC(ec)
{
    assert(itsBS);
    assert(itsEC);
    assert(itsBS->GetNumFunctions()>0);
    for (auto b:itsBS->Iterate<TOrbital_IBS<double> >())
    {
        uiwf_t wf(new Irrep_WF(b,Spin(Spin::None),acc.Create(b)));
        itsQN_WFs[wf->GetQNs()]=wf.get();
        itsIWFs.push_back(std::move(wf)); //Do the move last.
    }
};

UnPolarized_WF::~UnPolarized_WF()
{
    
}
//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the FillOrbitals member function
//  to fill up the orbitals with electrons.
//
void UnPolarized_WF::DoSCFIteration(Hamiltonian& ham,const DM_CD* cd)
{
    for (auto& w:itsIWFs) w->DoSCFIteration(ham,cd);
}

DM_CD* UnPolarized_WF::GetChargeDensity() const
{
    Composite_CD* cd = new Composite_CD();
    for (auto& w:itsIWFs) cd->Insert(w->GetChargeDensity());
    return cd;
}

const Orbitals* UnPolarized_WF::GetOrbitals(const Irrep_QNs& qns) const
{
    return const_cast<UnPolarized_WF*>(this)->GetOrbitals(qns);
}
Orbitals* UnPolarized_WF::GetOrbitals(const Irrep_QNs& qns) 
{
    assert(qns.ms==Spin::None);
    auto i=itsQN_WFs.find(qns);
    assert(i!=itsQN_WFs.end());
    return i->second->GetOrbitals();

}

UnPolarized_WF::iqns_t UnPolarized_WF::GetQNs() const
{
    iqns_t iqns;
    for (auto q:itsQN_WFs) iqns.push_back(q.first);
    return iqns;
}



void UnPolarized_WF::FillOrbitals()
{
    itsELevels.clear();
    for (auto& w:itsIWFs) 
        itsELevels.merge(w->FillOrbitals(itsEC),0.0001);
}

void UnPolarized_WF::DisplayEigen() const
{
    StreamableObject::SetToPretty();

    std::cout << "Alpha+Beta spin :" << std::endl;
    itsELevels.Report(std::cout);
}



