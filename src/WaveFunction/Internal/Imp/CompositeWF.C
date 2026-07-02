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
using SCFAccelerators::tSCFIrrepAccelerator;

// MOM (Maximum Overlap Method) occupation tracking is implemented but parked: for the current
// closed-shell cases the empty-irrep DIIS discriminator already gives clean convergence with no
// occupation flip, so MOM is not load-bearing (see doc/SCF_DIIS_SALC_notes.md).  The machinery is
// kept compiled; this flag gates only its call sites.  Flip to true for excited-state / hard cases.
constexpr bool EnableMOM=false;


template <class T> tCompositeWF<T>::tCompositeWF(const tbs_t<T>* bs,const ElectronConfiguration* ec,tSCFAccelerator<T>* acc )
    : itsBS(bs)
    , itsEC(ec)
    , itsAufbau(ec->UsesAufbau())
    , itsAccelerator(acc)
{
    assert(itsBS);
    assert(itsEC);
    assert(itsAccelerator);
    assert(itsBS->GetNumFunctions()>0);

};

template <class T> void tCompositeWF<T>::MakeIrrepWFs(Spin s)
{

    for (auto b:itsBS->template Iterate<tobs_t<T>>())
    {
        LASolver<T>* lasb=LASolver<T>::Factory(qchem::Cholesky);
        lasb->SetBasisOverlap(b->Overlap());
        // std::cout << "Minimum singular value for basis set overlap= " << blaze::min(lasb->Get_BS_Diagonal()) << std::endl;
        Irrep qns(b->GetIrrep(s));
        tSCFIrrepAccelerator<T>* acc=itsAccelerator->Create(lasb,qns,itsEC->GetN(qns));

        uiwf_t wf(new iwf_t(b,lasb,qns,acc));
        itsQNWFs[qns]=wf.get();
        itsSpinWFs[s].push_back(wf.get());
        itsIWFs.push_back(std::move(wf)); //Do the move last. wf is invalid after the move.
    }
}

template <class T> tCompositeWF<T>::~tCompositeWF()
{
    // delete itsAccelerator; NO!!!! SCFiterator deletes the accelerator.
};

//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the FillOrbitals member function
//  to fill up the orbitals with electrons.
//
// The cross-irrep view threaded to the dynamic Hamiltonian terms (doc/ERI4Rework.md §5.4): every
// participating irrep's orbital basis, in one place, so a term CAN exploit beyond-one-irrep structure.
// Currently inert -- no term reads it yet (stage 3a is the plumbing; the Coulomb/exchange consumer lands
// in 3b).  The basis list is spatial-irrep only (spin is applied per IrrepWF), matching MakeIrrepWFs.
template <class T> tHamiltonianContext<T> tCompositeWF<T>::MakeContext() const
{
    tHamiltonianContext<T> ctx;
    for (auto b:itsBS->template Iterate<tobs_t<T>>()) ctx.irrepBases.push_back(b);
    return ctx;
}

template <class T> void tCompositeWF<T>::DoSCFIteration(tHamiltonian<T>& ham,const tChargeDensity<T>* cd)
{
    tHamiltonianContext<T> ctx=MakeContext();
    for (auto& w:itsIWFs) w->CalculateH(ham,cd,ctx); //Feed F,D' into all the irre eccelerators.
    // Once the accelerator extrapolates, switch the molecular aufbau from eigenvalue order to MOM
    // (overlap): re-running a plain aufbau on the (non-physical) extrapolated Fock can flip the
    // occupation (e.g. a near-degenerate B2<->A1 swap in H2O), wrecking convergence.  MOM keeps the
    // occupied orbitals that match the previous (settled) ones by shape.  Sticky for the rest of run.
    [[maybe_unused]] const bool engaged = itsAccelerator->CalculateProjections();
    if constexpr (EnableMOM) if (engaged) itsMOMActive=true;
    for (auto& w:itsIWFs) w->DoSCFIteration();
}

// Iteration-0 seed: build the Fock from \a seed, diagonalize, fill, and hand back the first real density.
template <class T> tDM_CD<T>* tCompositeWF<T>::Init(tHamiltonian<T>& ham,const tChargeDensity<T>* seed, double mergeTol)
{
    DoSCFIteration(ham, seed);
    FillOrbitals(mergeTol);
    return GetChargeDensity();
}

// Build the Fock and have each irrep accelerator compute its (un-taken) step.  Returns true
// only if every irrep produced a geodesic step; false means at least one wants to diagonalize
// (the seed step) -- the caller should fall back to DoSCFIteration().
template <class T> bool tCompositeWF<T>::BuildFockAndComputeSteps(tHamiltonian<T>& ham,const tChargeDensity<T>* cd)
{
    tHamiltonianContext<T> ctx=MakeContext();
    for (auto& w:itsIWFs) w->CalculateH(ham,cd,ctx);
    bool allStepped=true;
    for (auto& w:itsIWFs) allStepped &= w->ComputeStep();
    return allStepped;
}

// Move every irrep's orbitals to geodesic fraction t (commit=false for a line-search trial)
// and refill, so GetChargeDensity() reflects the trial/updated orbitals.
template <class T> void tCompositeWF<T>::MoveOrbitals(double t, bool commit, double mergeTol)
{
    for (auto& w:itsIWFs) w->MoveOrbitals(t,commit);
    FillOrbitals(mergeTol);
}

template <class T> tDM_CD<T>* tCompositeWF<T>::GetChargeDensity(Spin s) const
{
    using qchem::ChargeDensity::tComposite_CD;
    auto i = itsSpinWFs.find(s);
    assert(i!=itsSpinWFs.end());
    tComposite_CD<T>* cd = new tComposite_CD<T>();
    for (auto& w:i->second) cd->Insert(w->GetChargeDensity());
    return cd;
}

template <class T> EnergyLevels tCompositeWF<T>::GetEnergyLevels (Spin s) const
{
    auto i = itsSpin_ELevels.find(s);
    assert(i!=itsSpin_ELevels.end());
    return i->second;
}

template <class T> const Orbitals* tCompositeWF<T>::GetOrbitals(const Irrep& qns) const
{
    return const_cast<tCompositeWF*>(this)->GetOrbitals(qns);
}
template <class T> Orbitals* tCompositeWF<T>::GetOrbitals(const Irrep& qns)
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

template <class T> typename tCompositeWF<T>::iqns_t tCompositeWF<T>::GetQNs() const
{
    iqns_t iqns;
    for (auto q:itsQNWFs) iqns.push_back(q.first);
    return iqns;
}



template <class T> void tCompositeWF<T>::FillOrbitals(double mergeTol)
{
    itsELevels.clear();
    itsSpin_ELevels.clear();
    if (itsAufbau) { FillOrbitalsAufbau(mergeTol); return; }
    for (auto& w:itsIWFs)                              // fixed per-irrep occupation (atoms etc.)
    {
        EnergyLevels els=w->FillOrbitals(itsEC);
        itsELevels.merge(els,mergeTol);
        Spin s=w->GetIrrep().ms;
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
template <class T> void tCompositeWF<T>::FillOrbitalsAufbau(double mergeTol)
{
    struct Slot { double key; iwf_t* w; double cap; };
    for (auto& [s, wfs] : itsSpinWFs)
    {
        if (wfs.empty()) continue;
        double Nc = (double)itsEC->GetN(wfs.front()->GetIrrep());   // total electrons in this spin channel

        std::map<Irrep,rvec_t> mom;                              // per-irrep MOM scores (empty if no ref), keyed by irrep
        if constexpr (EnableMOM) if (itsMOMActive) for (auto w : wfs) mom[w->GetIrrep()]=w->MOMScores();

        std::vector<Slot> slots;                                  // every orbital across the channel
        for (auto w : wfs)
        {
            size_t idx=0;
            const rvec_t& sc = mom[w->GetIrrep()];               // empty unless MOM active & referenced
            for (auto o : w->GetOrbitals()->Iterate())
            {
                // MOM: higher overlap = occupy first (unreferenced/empty irrep scores 0).
                // Aufbau: lower eigenvalue = occupy first, so key on -energy for a common "bigger wins".
                double key = itsMOMActive ? (idx<sc.size() ? sc[idx] : 0.0) : -o->GetEigenEnergy();
                slots.push_back({key, w, (double)o->GetDegeneracy()});
                ++idx;
            }
        }
        std::sort(slots.begin(), slots.end(), [](const Slot& a, const Slot& b){return a.key>b.key;});

        std::map<Irrep,double>& ne = itsAufbauNe[s];
        for (auto w : wfs) ne[w->GetIrrep()]=0.0;
        assert(ne.size()==wfs.size() && "aufbau key collision: irreps must be unique within a spin channel");
        double rem = Nc;
        for (const auto& sl : slots)                              // fill highest-priority first
        {
            if (rem<=0.0) break;
            double take = std::min(sl.cap, rem);
            ne[sl.w->GetIrrep()] += take;
            rem -= take;
        }

        for (auto w : wfs)                                        // occupy + collect energy levels
        {
            EnergyLevels els = w->FillOrbitals(ne[w->GetIrrep()]);
            itsELevels.merge(els, mergeTol);
            itsSpin_ELevels[s].merge(els, mergeTol);
        }
        if constexpr (EnableMOM) for (auto w : wfs) w->CaptureMOMReference(); // reference for next iteration's MOM
    }
}

template class tCompositeWF<double>;
template class tCompositeWF<dcmplx>;

} //namespace
