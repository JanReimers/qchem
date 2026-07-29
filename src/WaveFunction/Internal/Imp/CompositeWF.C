// File: CompositeWF.H  Wave function as a list of Irrep wave functions.
module;
#include <cassert>
#include <cmath>
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include "tabulate/table.hpp"
module qchem.WaveFunction.Internal.CompositeWF;
import qchem.WaveFunction.Types;
import qchem.SCFAccelerator;
import qchem.CompositeCD;
import qchem.ElectronConfiguration;
import qchem.LASolver;
import qchem.Orbitals;          // Orbital (eigen-energy, degeneracy) for the aufbau
import qchem.Reporting;         // the run report -- this loop owns the basis.perIrrep / basis.removed rows

namespace qchem::WaveFunction
{

using std::cerr;
using std::endl;
using SCFAccelerators::tSCFIrrepAccelerator;

// MOM (Maximum Overlap Method) occupation tracking is configured per-run from SCFParams (UseMOM/
// MOMStartIter), pushed in by the SCFIterator via SetMOM (which forwards to every irrep WF).  For the
// closed-shell molecular cases the empty-irrep DIIS discriminator already gives clean convergence with no
// occupation flip, so MOM is OFF by default; it is turned on for the periodic NaF Γ map, where a
// giant-response diffuse virtual periodically dives below the Fermi edge and plain aufbau captures it
// (an occupation swap → energy spike, doc/GPWPlan §0b″).  The ACTIVE, tested path is the crystal's
// WITHIN-irrep MOM in tIrrepWF::FillOrbitals; this file's cross-irrep aufbau MOM is the parked molecular one.


template <class T> tCompositeWF<T>::tCompositeWF(const tbs_t<T>* bs,const ElectronConfiguration* ec,tSCFAccelerator<T>* acc,
                                                 qchem::Ortho basisOrtho, double basisOrthoTol )
    : itsBS(bs)
    , itsEC(ec)
    , itsBasisOrtho(basisOrtho)
    , itsBasisOrthoTol(basisOrthoTol)
    , itsAufbau(ec->UsesAufbau())
    , itsGlobalFermi(ec->UsesGlobalFermi())
    , itsAccelerator(acc)
{
    assert(itsBS);
    assert(itsEC);
    assert(itsAccelerator);
    assert(itsBS->GetNumFunctions()>0);
    assert(!(itsAufbau && itsGlobalFermi) && "aufbau and global-Fermi fills are mutually exclusive");
};

// Compact irrep label for the report (symmetry symbol + spin arrow), via the Streamable Write().
static std::string IrrepLabel(const Irrep& q) { std::ostringstream os; q.Write(os); return os.str(); }

template <class T> void tCompositeWF<T>::MakeIrrepWFs(Spin s)
{
    namespace rpt = qchem::report;
    for (auto b:itsBS->template Iterate<tobs_t<T>>())
    {
        const hmat_t<T>& S = b->Overlap();
        LASolver<T>* lasb=LASolver<T>::Factory(itsBasisOrtho, itsBasisOrthoTol);
        Irrep qns(b->GetIrrep(s));

        // Emit the report's basis.perIrrep row ONLY when an ancestor opened a "basis" section (the
        // orchestrator's choice).  The molecular Calculation opens it around this construction; a GPW
        // orchestrator that vets conditioning in a PRE-FLIGHT (before this WF build) leaves it closed, so
        // this stays silent and does not duplicate.  When open, the LASolver's conditioning write (five
        // layers down in SetBasisOverlap) lands in THIS row with no irrep identity threaded to it.
        const bool reporting = rpt::InSection("basis");
        if (reporting)
        {
            rpt::Row row("perIrrep");
            rpt::Set("irrep",      IrrepLabel(qns));
            rpt::Set("nFunctions", (long)b->GetNumFunctions());
            lasb->SetBasisOverlap(S);
        }
        else
            lasb->SetBasisOverlap(S);   // ortho only -- no basis context, no report
        // basis.exponents (report-only, Verbose console): the basis SERIALIZES its own radial parameters into
        // this row -- {irrep, values:[...]} -- so a consumer (CLIapps/valgen) can label the basis.usage rows by
        // exponent instead of index.  The evaluator writes the values (a no-op for a non-exponent basis); the
        // exponents stay encapsulated (json sink, never a getter).
        if (reporting)
        {
            rpt::Row r("exponents");
            rpt::Set("irrep", IrrepLabel(qns));
            b->EmitRadialReport();
        }
        // basis.removed (report-only detector, doc/GPWPlan1.md §4a): the redundant AO functions this irrep
        // carries.  {irrep, index} for now (exponent/atom naming awaits a per-function metadata accessor).
        if (reporting)
            for (size_t idx : qchem::PivotedCholeskyDrops<T>(S))
            {
                rpt::Row r("removed");
                rpt::Set("irrep", IrrepLabel(qns));
                rpt::Set("index", (long)idx);
            }

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
template <class T> void tCompositeWF<T>::DoSCFIteration(tHamiltonian<T>& ham,const tChargeDensity<T>* cd)
{
    // itsBS (the whole/composite basis) IS the cross-irrep view a dynamic term may exploit: Iterate<tobs_t>()
    // over it yields every irrep block (doc/ERI4Rework.md §5.4).  Static terms and most dynamic terms ignore it.
    for (auto& w:itsIWFs) w->CalculateH(ham,cd,itsBS); //Feed F,D' into all the irre eccelerators.
    // CalculateProjections() has the DIIS side effect (it accumulates the extrapolation history) so it
    // must run every iteration regardless of MOM -- keep the call.  MOM activation is NO LONGER gated on
    // the accelerator engaging (that was the parked molecular heuristic, and NaF's Null accelerator never
    // engages); it is switched on in FillOrbitalsAufbau as soon as a reference occupied subspace exists.
    itsAccelerator->CalculateProjections();
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
    for (auto& w:itsIWFs) w->CalculateH(ham,cd,itsBS);
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

// Emit the run report's `basis.usage` section: per-function occupation-weighted Mulliken populations across
// every irrep WF (both spins for a polarized run).  Rows are {irrep, index, pop} -- the WF reports only what
// it knows (populations); a consumer (CLIapps/valgen) joins the exponent from basis.exponents and draws the
// heat bar.  A LATE addendum to the already-rendered `basis` section (EmitAt, absolute path), gated to
// Detail::Verbose on the console but ALWAYS recorded in the json.  Called by the SCFIterator post-convergence.
template <class T> void tCompositeWF<T>::EmitBasisUsage() const
{
    namespace rpt = qchem::report;
    if (rpt::Depth()==0) return;                           // no run open -> nothing to record or render

    rpt::json rows=rpt::json::array();
    for (const auto& w : itsIWFs)
    {
        const rvec_t P = w->GetBasisPopulations();
        std::ostringstream os; w->GetIrrep().Write(os);   // "s"/"p"/... (+ spin arrow when polarized)
        const std::string lbl=os.str();
        for (std::size_t i=0;i<P.size();++i)
            rows.push_back(rpt::json{{"irrep",lbl},{"index",(long)i},{"pop",P[i]}});
    }
    rpt::EmitAt("basis","usage",rows,rpt::Detail::Verbose);
}



// Configure MOM for this run (from SCFParams, via the SCFIterator).  Store our own copy (used by the
// molecular cross-irrep aufbau) and forward to every irrep WF (the crystal's within-irrep MOM lives there).
template <class T> void tCompositeWF<T>::SetMOM(bool useMOM, int startIter)
{
    itsUseMOM = useMOM;
    for (auto& w : itsIWFs) w->SetMOM(useMOM, startIter);
}

// Fermi smearing (doc/GPWPlan1.md 4b): push the same kT (and MOM-mask penalty) into every irrep block, and keep
// our own copy of kT (the global-μ metal fill solves ONE μ at the composite level, so it needs kT here too).
// Per-block μ (insulator/Γ) is solved inside each block's FillOrbitals; global μ across the mesh in
// FillOrbitalsGlobalFermi.  No-op at kT=0.
template <class T> void tCompositeWF<T>::SetSmearing(double kT, double momPenalty)
{
    itsSmearingkT=kT;
    for (auto& w : itsIWFs) w->SetSmearing(kT, momPenalty);
}

// The run's Mermin free-energy term −TS: the SUM over blocks, each already BZ-weighted (tIrrepWF::GetEntropyTerm
// applies w_k), so this is Σ_k w_k(−T S_k) -- consistent with the BZ-weighted E.  0 with no smearing.  Stamped
// into EnergyBreakdown by the SCFIterator.
template <class T> double tCompositeWF<T>::GetEntropyTerm() const
{
    double minusTS=0.0;
    for (auto& w : itsIWFs) minusTS += w->GetEntropyTerm();
    return minusTS;
}

// Grid-continuation (doc/GPWPlan §0e): per irrep, hand this WF's irrep block the CONVERGED source WF's
// occupied orbitals for the SAME irrep as its fixed MOM reference.  from.GetOrbitals(irr) is the source's
// physical occupied subspace; the orthonormal metric matches (grid-independent Bloch overlap), so the C'
// transfer verbatim.  Call after construction (Init) and before Iterate; effective with SCFParams::UseMOM.
template <class T> void tCompositeWF<T>::AdoptMOMReference(const tWaveFunction<T>& from)
{
    for (auto& w : itsIWFs) w->AdoptMOMReference(*from.GetOrbitals(w->GetIrrep()));
}

// 0h MOM guard actuator: forward the release to every irrep WF (each drops its reference + re-arms capture).
template <class T> void tCompositeWF<T>::ReleaseMOMReference()
{
    for (auto& w : itsIWFs) w->ReleaseMOMReference();
}

template <class T> void tCompositeWF<T>::FillOrbitals(double mergeTol)
{
    itsELevels.clear();
    itsSpin_ELevels.clear();
    if (itsGlobalFermi) { FillOrbitalsGlobalFermi(mergeTol); return; }  // metal: one μ across the k-mesh
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
    // Snapshot MOM state at ENTRY: activating it below (after the reference capture) must not leak into
    // this same call's LATER spin channels (each channel's reference is only captured at its own end).
    const bool useMOM = itsMOMActive;
    for (auto& [s, wfs] : itsSpinWFs)
    {
        if (wfs.empty()) continue;
        double Nc = (double)itsEC->GetN(wfs.front()->GetIrrep());   // total electrons in this spin channel

        std::map<Irrep,rvec_t> mom;                              // per-irrep MOM scores (empty if no ref), keyed by irrep
        if (useMOM) for (auto w : wfs) mom[w->GetIrrep()]=w->MOMScores();

        std::vector<Slot> slots;                                  // every orbital across the channel
        for (auto w : wfs)
        {
            size_t idx=0;
            const rvec_t& sc = mom[w->GetIrrep()];               // empty unless MOM active & referenced
            for (auto o : w->GetOrbitals()->Iterate())
            {
                // MOM: higher overlap = occupy first (unreferenced/empty irrep scores 0).
                // Aufbau: lower eigenvalue = occupy first, so key on -energy for a common "bigger wins".
                double key = useMOM ? (idx<sc.size() ? sc[idx] : 0.0) : -o->GetEigenEnergy();
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
        // (The per-irrep MOM reference is captured inside tIrrepWF::FillOrbitals above, so no capture here.)
    }
    // Molecular CROSS-irrep MOM (this function) is the PARKED path: no molecular test enables MOM, and it
    // has NOT been re-validated against the new delayed reference capture (tIrrepWF captures at MOMStartIter,
    // so useMOM would score against empty references for the first few fills).  The ACTIVE, tested path is
    // the crystal's WITHIN-irrep MOM in tIrrepWF::FillOrbitals.  Kept here for the molecular hard cases.
    if (itsUseMOM) itsMOMActive=true;
}

// GLOBAL-μ METAL FILL (doc/GPWPlan1.md item 3): one chemical potential across the Bloch k-mesh, so charge
// sloshes between k-points (a partially-filled band -- the metal case per-block fixed occupation cannot
// represent).  Gather every orbital across the channel's k-blocks tagged with its BZ weight w_k, solve ONE μ
// on Σ_k w_k Σ_i g_i f((ε_i,k−μ)/kT)=N_total (the whole-mesh total), then set every block at that μ.  The
// weight rides in the level capacity g so the SAME unweighted qchem::Orbitals::FermiLevel bisector serves the
// per-block (w=1) and cross-k cases.  Reduces EXACTLY to the per-block Fermi at a single k (one block, w=1).
template <class T> void tCompositeWF<T>::FillOrbitalsGlobalFermi(double mergeTol)
{
    assert(itsSmearingkT>0.0 && "global-μ metal fill needs SmearingkT>0 (a partially-filled band has no integer fill)");
    for (auto& [s, wfs] : itsSpinWFs)
    {
        if (wfs.empty()) continue;
        const double Ntot=(double)itsEC->GetN(wfs.front()->GetIrrep());   // whole-mesh total for this spin channel

        // Aggregate the channel's levels: effective (bare) energy ε and BZ-weighted capacity w_k·g_i.
        size_t ntot=0;
        for (auto w : wfs) ntot += w->GetOrbitals()->GetNumOrbitals();
        rvec_t e(ntot), g(ntot);
        size_t idx=0;
        for (auto w : wfs)
        {
            const double wk=w->GetIrrep().sym->GetWeight();              // SAME weight GetChargeDensity applies to D
            for (auto o : w->GetOrbitals()->Iterate())
            {
                e[idx]=o->GetEigenEnergy();
                g[idx]=wk*(double)o->GetDegeneracy();
                ++idx;
            }
        }
        const double mu=qchem::Orbitals::FermiLevel(e,g,Ntot,itsSmearingkT);   // one μ across the whole mesh

        // Instrument GPW_METALTRACE (cf. GPW_GDMTRACE): dump μ and the per-k electron count n_k so the charge
        // sloshing is visible (n_k varies; Σ_k w_k n_k = N_total).  Off by default; per-iteration when set.
        const bool trace=(bool)std::getenv("GPW_METALTRACE");
        if (trace) std::cout<<"[metal] μ="<<mu<<" kT="<<itsSmearingkT<<" Ntot="<<Ntot<<"  per-k n_k:\n";
        double wsum=0.0;
        for (auto w : wfs)                                               // set every block at the shared μ
        {
            EnergyLevels els=w->FillOrbitalsAtMu(mu);
            if (trace)
            {
                double nk=0.0; for (auto o:w->GetOrbitals()->Iterate()) nk+=o->GetOccupation();
                const double wk=w->GetIrrep().sym->GetWeight();
                wsum+=wk*nk;
                std::cout<<"[metal]   k="<<w->GetIrrep()<<" w="<<wk<<" n_k="<<nk<<" (w·n="<<wk*nk<<")\n";
            }
            itsELevels.merge(els,mergeTol);
            itsSpin_ELevels[s].merge(els,mergeTol);
        }
        if (trace) std::cout<<"[metal] Σ w_k n_k = "<<wsum<<"\n";
    }
}

template class tCompositeWF<double>;
template class tCompositeWF<dcmplx>;

} //namespace
