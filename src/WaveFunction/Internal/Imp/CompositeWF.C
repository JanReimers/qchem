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

// MOM (Maximum Overlap Method) occupation tracking is implemented but parked: for the current
// closed-shell cases the empty-irrep DIIS discriminator already gives clean convergence with no
// occupation flip, so MOM is not load-bearing (see doc/SCF_DIIS_SALC_notes.md).  The machinery is
// kept compiled; this flag gates only its call sites.  Flip to true for excited-state / hard cases.
constexpr bool EnableMOM=false;


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
    // Once the accelerator extrapolates, switch the molecular aufbau from eigenvalue order to MOM
    // (overlap): re-running a plain aufbau on the (non-physical) extrapolated Fock can flip the
    // occupation (e.g. a near-degenerate B2<->A1 swap in H2O), wrecking convergence.  MOM keeps the
    // occupied orbitals that match the previous (settled) ones by shape.  Sticky for the rest of run.
    [[maybe_unused]] const bool engaged = itsAccelerator->CalculateProjections();
    if constexpr (EnableMOM) if (engaged) itsMOMActive=true;
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

// Molecular aufbau: per spin channel, pick which orbitals across all irreps are occupied, then
// fill each irrep with its resulting electron count.  The point-group irrep an occupied MO lands
// in is an OUTPUT of the SCF (unlike an atom's fixed l-occupation).  Two selection rules:
//   * plain aufbau (default): occupy the globally-lowest eigenvalues;
//   * MOM (once the accelerator engages): occupy the orbitals with the largest overlap onto the
//     previous iteration's occupied subspace, so a near-degenerate cross-irrep pair on the
//     (non-physical) extrapolated Fock cannot flip the configuration.
// Either way the per-irrep references are re-captured at the end for the next iteration's MOM.
void CompositeWF::FillOrbitalsAufbau(double mergeTol)
{
    struct Slot { double key; IrrepWF* w; double cap; };
    for (auto& [s, wfs] : itsSpinWFs)
    {
        if (wfs.empty()) continue;
        double Nc = (double)itsEC->GetN(wfs.front()->GetQNs());   // total electrons in this spin channel

        std::map<IrrepWF*,rvec_t> mom;                            // per-irrep MOM scores (empty if no ref)
        if constexpr (EnableMOM) if (itsMOMActive) for (auto w : wfs) mom[w]=w->MOMScores();

        std::vector<Slot> slots;                                  // every orbital across the channel
        for (auto w : wfs)
        {
            size_t idx=0;
            const rvec_t& sc = mom[w];                            // empty unless MOM active & referenced
            for (auto o : w->GetOrbitals()->Iterate<qchem::Orbitals::Orbital>())
            {
                // MOM: higher overlap = occupy first (unreferenced/empty irrep scores 0).
                // Aufbau: lower eigenvalue = occupy first, so key on -energy for a common "bigger wins".
                double key = itsMOMActive ? (idx<sc.size() ? sc[idx] : 0.0) : -o->GetEigenEnergy();
                slots.push_back({key, w, (double)o->GetDegeneracy()});
                ++idx;
            }
        }
        std::sort(slots.begin(), slots.end(), [](const Slot& a, const Slot& b){return a.key>b.key;});

        std::map<IrrepWF*,double>& ne = itsAufbauNe[s];
        for (auto w : wfs) ne[w]=0.0;
        double rem = Nc;
        for (const auto& sl : slots)                              // fill highest-priority first
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
        if constexpr (EnableMOM) for (auto w : wfs) w->CaptureMOMReference(); // reference for next iteration's MOM
    }
}

} //namespace
