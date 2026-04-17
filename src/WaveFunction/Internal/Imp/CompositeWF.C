// File: CompositeWF.H  Wave function as a list of Irrep wave functions.
module;
#include <cassert>
#include <memory>
#include "tabulate/table.hpp"

module qchem.WaveFunction.Internal.CompositeWF;
import qchem.SCFAccelerator;
import qchem.BasisSet;
import qchem.IrrepBasisSet;
import qchem.CompositeCD;
import qchem.Symmetry.ElectronConfiguration;
import qchem.LASolver;


LAParams DefaultLAP({qchem::Cholsky,1e-12});


CompositeWF::CompositeWF(const BasisSet* bs,const ElectronConfiguration* ec,SCFAccelerator* acc )
    : itsBS(bs)
    , itsEC(ec)
    , itsAccelerator(acc)
    , itsLAParams(DefaultLAP) //gcc-15.0.1 segfault here

{
    assert(itsBS);
    assert(itsEC);
    assert(itsAccelerator);
    assert(itsBS->GetNumFunctions()>0);
    
};

void CompositeWF::MakeIrrepWFs(Spin s)
{

    for (auto b:itsBS->Iterate<Orbital_IBS<double> >())
    {
        LASolver<double>* lasb=LASolver<double>::Factory(itsLAParams.BasisOrthoAlgorithm,itsLAParams.TruncationTolerance);
        lasb->SetBasisOverlap(b->Overlap());
    // std::cout << "Minimum singular value for basis set overlap= " << Min(las->Get_BS_Diagonal()) << std::endl;
        Irrep_QNs qns(b->GetIrrep(s));
        SCFIrrepAccelerator* acc=itsAccelerator->Create(lasb,qns,itsEC->GetN(qns));
        
        uiwf_t wf(new IrrepWF(b,lasb,qns,acc));
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





