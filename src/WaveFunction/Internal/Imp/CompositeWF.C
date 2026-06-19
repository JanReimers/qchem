// File: CompositeWF.H  Wave function as a list of Irrep wave functions.
module;
#include <cassert>
#include <memory>
#include <vector>
#include <map>
#include <algorithm>
#include "tabulate/table.hpp"
module qchem.WaveFunction.Internal.CompositeWF;
import qchem.WaveFunction.Types;
import qchem.SCFAccelerator;
import qchem.CompositeCD;
import qchem.ElectronConfiguration;
import qchem.LASolver;
import qchem.Orbitals;          // Orbital (eigen-energy, degeneracy) for the aufbau

namespace qchem::WaveFunction
{

using std::cerr;
using std::endl;

LAParams DefaultLAP({qchem::Cholsky,1e-12});


CompositeWF::CompositeWF(const bs_t* bs,const ElectronConfiguration* ec,SCFAccelerator* acc )
    : itsBS(bs)
    , itsEC(ec)
    , itsAufbau(ec->UsesAufbau())
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
        Irrep qns(b->GetIrrep(s));
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

// Build the Fock and have each irrep accelerator compute its (un-taken) step.  Returns true
// only if every irrep produced a geodesic step; false means at least one wants to diagonalize
// (the seed step) -- the caller should fall back to DoSCFIteration().
bool CompositeWF::BuildFockAndComputeSteps(Hamiltonian& ham,const DM_CD* cd)
{
    for (auto& w:itsIWFs) w->CalculateH(ham,cd);
    bool allStepped=true;
    for (auto& w:itsIWFs) allStepped &= w->ComputeStep();
    return allStepped;
}

// Move every irrep's orbitals to geodesic fraction t (commit=false for a line-search trial)
// and refill, so GetChargeDensity() reflects the trial/updated orbitals.
void CompositeWF::MoveOrbitals(double t, bool commit, double mergeTol)
{
    for (auto& w:itsIWFs) w->MoveOrbitals(t,commit);
    FillOrbitals(mergeTol);
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

const Orbitals* CompositeWF::GetOrbitals(const Irrep& qns) const
{
    return const_cast<CompositeWF*>(this)->GetOrbitals(qns);
}
Orbitals* CompositeWF::GetOrbitals(const Irrep& qns) 
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
    if (itsAufbau) { FillOrbitalsAufbau(mergeTol); return; }
    for (auto& w:itsIWFs)                              // fixed per-irrep occupation (atoms etc.)
    {
        EnergyLevels els=w->FillOrbitals(itsEC);
        itsELevels.merge(els,mergeTol);
        Spin s=w->GetQNs().ms;
        itsSpin_ELevels[s].merge(els,mergeTol);
    }
}

// Molecular aufbau: per spin channel, occupy the globally-lowest orbitals across all irreps,
// then fill each irrep with its resulting electron count.  The point-group irrep an occupied
// MO lands in is an OUTPUT of the SCF (unlike an atom's fixed l-occupation), so the per-irrep
// counts are recomputed every iteration from the current eigenvalues.
void CompositeWF::FillOrbitalsAufbau(double mergeTol)
{
    struct Slot { double e; IrrepWF* w; double cap; };
    for (auto& [s, wfs] : itsSpinWFs)
    {
        if (wfs.empty()) continue;
        double Nc = (double)itsEC->GetN(wfs.front()->GetQNs());   // total electrons in this spin channel

        std::vector<Slot> slots;                                  // every orbital across the channel
        for (auto w : wfs)
            for (auto o : w->GetOrbitals()->Iterate<qchem::Orbitals::Orbital>())
                slots.push_back({o->GetEigenEnergy(), w, (double)o->GetDegeneracy()});
        std::sort(slots.begin(), slots.end(), [](const Slot& a, const Slot& b){return a.e<b.e;});

        std::map<IrrepWF*,double> ne;                             // per-irrep electron count
        for (auto w : wfs) ne[w]=0.0;
        double rem = Nc;
        for (const auto& sl : slots)                              // fill lowest-first
        {
            if (rem<=0.0) break;
            double take = std::min(sl.cap, rem);
            ne[sl.w] += take;
            rem -= take;
        }

        for (auto w : wfs)                                        // occupy + collect energy levels
        {
            EnergyLevels els = w->FillOrbitals(ne[w]);
            itsELevels.merge(els, mergeTol);
            itsSpin_ELevels[s].merge(els, mergeTol);
        }
    }
}

} //namespace
