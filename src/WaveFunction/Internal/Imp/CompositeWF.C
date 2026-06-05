// File: CompositeWF.H  Wave function as a list of Irrep wave functions.
module;
#include <cassert>
#include <memory>
#include "tabulate/table.hpp"
// #include <blaze/Math.h>
module qchem.WaveFunction.Internal.CompositeWF;
import qchem.WaveFunction.Types;
import qchem.SCFAccelerator;
import qchem.CompositeCD;
import qchem.Symmetry.ElectronConfiguration;
import qchem.LASolver;

namespace qchem::WaveFunction
{

using std::cerr;
using std::endl;

LAParams DefaultLAP({qchem::Cholsky,1e-12});


CompositeWF::CompositeWF(const bs_t* bs,const ElectronConfiguration* ec,SCFAccelerator* acc )
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

    for (auto b:itsBS->Iterate<obs_t>())
    {
        LASolver<double>* lasb=LASolver<double>::Factory(itsLAParams.BasisOrthoAlgorithm,itsLAParams.TruncationTolerance);
        lasb->SetBasisOverlap(b->Overlap());
        // std::cout << "Minimum singular value for basis set overlap= " << blaze::min(lasb->Get_BS_Diagonal()) << std::endl;
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
    using qchem::ChargeDensity::Composite_CD;
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
    if (i==itsQNWFs.end())
    {
        cerr << "CompositeWF::GetOrbitals cannot find orbital: " << qns << endl;
        cerr << "  Known orbitals are:" << endl;
        for (auto i:itsQNWFs ) cerr << "    " << i.first << endl;
        assert(false);
    }
    // assert(i!=itsQNWFs.end());
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

} //namespace
