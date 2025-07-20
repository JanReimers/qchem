// File: CompositeWF.H  Wave function as a list of Irrep wave functions.
module;
#include <cassert>
#include <memory>
#include <SCFAccelerator/SCFAccelerator.H>
module qchem.WaveFunction.Internal.CompositeWF;
import qchem.BasisSet;
import qchem.Irrep_BS;
import qchem.CompositeCD;
import qchem.Symmetry.ElectronConfiguration;

CompositeWF::CompositeWF(const BasisSet* bs,const ElectronConfiguration* ec,SCFAccelerator* acc )
    : itsBS(bs)
    , itsEC(ec)
    , itsAccelerator(acc)
{
    assert(itsBS);
    assert(itsEC);
    assert(itsAccelerator);
    assert(itsBS->GetNumFunctions()>0);
    
};

void CompositeWF::MakeIrrepWFs(Spin s)
{

    for (auto b:itsBS->Iterate<TOrbital_IBS<double> >())
    {
        auto las=b->CreateSolver(); //IrrepWF will own delete this thing.
        Irrep_QNs qns(s,b->GetSymmetry());
        SCFIrrepAccelerator* acc=itsAccelerator->Create(las,qns,itsEC->GetN(qns));
        
        uiwf_t wf(new IrrepWF(b,las,qns,acc));
        itsQNWFs[qns]=wf.get();
        itsSpinWFs[s].push_back(wf.get());
        itsIWFs.push_back(std::move(wf)); //Do the move last. wf is invalid after the move.
    }
}

CompositeWF::~CompositeWF() 
{
    // delete itsAccelerator; NO!!!! SCFiterator deletes the accelerator.
};

//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the FillOrbitals member function
//  to fill up the orbitals with electrons.
//
void CompositeWF::DoSCFIteration(Hamiltonian& ham,const DM_CD* cd)
{
    for (auto& w:itsIWFs) w->CalculateH(ham,cd); //Feed F,D' into all the irre eccelerators.
    itsAccelerator->CalculateProjections();
    for (auto& w:itsIWFs) w->DoSCFIteration();
}

DM_CD* CompositeWF::GetChargeDensity(Spin s) const
{
    auto i = itsSpinWFs.find(s);
    assert(i!=itsSpinWFs.end());
    Composite_CD* cd = new Composite_CD();
    for (auto& w:i->second) cd->Insert(w->GetChargeDensity());
    return cd;
}

EnergyLevels CompositeWF::GetEnergyLevels (Spin s) const 
{
    auto i = itsSpin_ELevels.find(s);
    assert(i!=itsSpin_ELevels.end());
    return i->second;
} 

const Orbitals* CompositeWF::GetOrbitals(const Irrep_QNs& qns) const
{
    return const_cast<CompositeWF*>(this)->GetOrbitals(qns);
}
Orbitals* CompositeWF::GetOrbitals(const Irrep_QNs& qns) 
{
    auto i=itsQNWFs.find(qns);
    assert(i!=itsQNWFs.end());
    return i->second->GetOrbitals();

}

CompositeWF::iqns_t CompositeWF::GetQNs() const
{
    iqns_t iqns;
    for (auto q:itsQNWFs) iqns.push_back(q.first);
    return iqns;
}



void CompositeWF::FillOrbitals(double mergeTol)
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





