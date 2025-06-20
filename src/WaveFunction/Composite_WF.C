// File: Composite_WF.H  Wave function as a list of Irrep wave functions.

#include "Composite_WF.H"
#include "Irrep_WF.H"
#include <ChargeDensity/CompositeCD.H>
#include <BasisSet/BasisSet.H>
#include <BasisSet/Irrep_BS.H>
#include <Symmetry/ElectronConfiguration.H>
#include <SCFAccelerator/SCFAccelerator.H>
#include <cassert>

Composite_WF::Composite_WF(const BasisSet* bs,const ElectronConfiguration* ec,SCFAccelerator* acc )
    : itsBS(bs)
    , itsEC(ec)
    , itsAccelerator(acc)
{
    assert(itsBS);
    assert(itsEC);
    assert(itsAccelerator);
    assert(itsBS->GetNumFunctions()>0);
    
};

void Composite_WF::MakeIrrep_WFs(Spin s)
{

    for (auto b:itsBS->Iterate<TOrbital_IBS<double> >())
    {
        auto las=b->CreateSolver(); //Irrep_WF will own delete this thing.
        Irrep_QNs qns(s,b->GetSymmetry());
        SCFIrrepAccelerator* acc=itsAccelerator->Create(las,qns,itsEC->GetN(qns));
        
        uiwf_t wf(new Irrep_WF(b,las,qns,acc));
        itsQN_WFs[qns]=wf.get();
        itsSpin_WFs[s].push_back(wf.get());
        itsIWFs.push_back(std::move(wf)); //Do the move last. wf is invalid after the move.
    }
}

Composite_WF::~Composite_WF() 
{
    // delete itsAccelerator; NO!!!! SCFiterator deletes the accelerator.
};

//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the FillOrbitals member function
//  to fill up the orbitals with electrons.
//
void Composite_WF::DoSCFIteration(Hamiltonian& ham,const DM_CD* cd)
{
    for (auto& w:itsIWFs) w->CalculateH(ham,cd); //Feed F,D' into all the irre eccelerators.
    itsAccelerator->CalculateProjections();
    for (auto& w:itsIWFs) w->DoSCFIteration();
}

DM_CD* Composite_WF::GetChargeDensity(Spin s) const
{
    auto i = itsSpin_WFs.find(s);
    assert(i!=itsSpin_WFs.end());
    Composite_CD* cd = new Composite_CD();
    for (auto& w:i->second) cd->Insert(w->GetChargeDensity());
    return cd;
}

EnergyLevels Composite_WF::GetEnergyLevels (Spin s) const 
{
    auto i = itsSpin_ELevels.find(s);
    assert(i!=itsSpin_ELevels.end());
    return i->second;
} 

const Orbitals* Composite_WF::GetOrbitals(const Irrep_QNs& qns) const
{
    return const_cast<Composite_WF*>(this)->GetOrbitals(qns);
}
Orbitals* Composite_WF::GetOrbitals(const Irrep_QNs& qns) 
{
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



void Composite_WF::FillOrbitals(double mergeTol)
{
    itsELevels.clear();
    itsSpin_ELevels.clear();
    for (auto& w:itsIWFs) 
    {
        EnergyLevels els=w->FillOrbitals(itsEC);
        itsELevels.merge(els,mergeTol);
        Spin s=w->GetQNs().ms;
        itsSpin_ELevels[s].merge(els,mergeTol);
    }
        
}





