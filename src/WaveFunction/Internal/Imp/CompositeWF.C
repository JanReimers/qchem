// File: CompositeWF.H  Wave function as a list of Irrep wave functions.
module;
#include <cassert>
#include <cmath>
#include <memory>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <variant>
#include <stdexcept>
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
// MOMStartIter) on the SCFIterator's OccupationPolicy slot, which every fill consults (V1.11 inc 3).  For the
// closed-shell molecular cases the empty-irrep DIIS discriminator already gives clean convergence with no
// occupation flip, so MOM is OFF by default; it is turned on for the periodic NaF Γ map, where a
// giant-response diffuse virtual periodically dives below the Fermi edge and plain aufbau captures it
// (an occupation swap → energy spike, doc/GPWPlan §0b″).  The ACTIVE, tested path is the crystal's
// WITHIN-irrep MOM in tIrrepWF::FillOrbitals; this file's cross-irrep aufbau MOM is the parked molecular one.


template <class T> tCompositeWF<T>::tCompositeWF(const tbs_t<T>* bs,const ElectronConfiguration* ec,SCFAccelerator* acc,
                                                 qchem::Ortho basisOrtho, double basisOrthoTol )
    : itsBS(bs)
    , itsEC(ec)
    , itsBasisOrtho(basisOrtho)
    , itsBasisOrthoTol(basisOrthoTol)
    , itsPartition(ec->GetPartition())
    , itsAccelerator(acc)
{
    assert(itsBS);
    assert(itsEC);
    assert(itsAccelerator);
    assert(itsBS->GetNumFunctions()>0);
    // (The old aufbau-vs-globalFermi mutual exclusion is now UNREPRESENTABLE: ranked-vs-prescribed is one
    //  flag on one partition, not two competing modes.)
    assert(!(itsPartition.ranksIntegerFill && !itsPartition.spansSpatial)
           && "a ranked integer fill IS a cross-block aufbau -- it needs a spatial-spanning reservoir");
    assert(!(itsPartition.ranksIntegerFill && itsPartition.spansSpin)
           && "the ranked (aufbau) fill fills each spin channel to its own count -- it cannot pool the channels");
};

// Compact irrep label for the report (symmetry symbol + spin arrow), via the Streamable Write().
static std::string IrrepLabel(const Irrep& q) { std::ostringstream os; q.Write(os); return os.str(); }

// ---- Child-slot access (Step 2b; same pattern as tComposite_CD's 2a) --------------------------------
// Scalar-independent views, single-source generic lambdas over both alternatives:
static iwf_ref_t     RefOf     (const iwf_child_t& c) {return std::visit([](const auto& p){return iwf_ref_t(p.get());}, c);}
static const Irrep&  IrrepOf   (const iwf_ref_t&   w) {return std::visit([](auto p)->const Irrep&{return p->GetIrrep();}, w);}
static Orbitals*     OrbitalsOf(const iwf_ref_t&   w) {return std::visit([](auto p)->Orbitals*   {return p->GetOrbitals();}, w);}
// Same-scalar view for the T-TYPED calls (the Hamiltonian and the OccupationPolicy carry the composite
// FACE's T).  The cross arm is unreachable today (children are built same-T); it becomes reachable -- and
// gets the genuine narrowing -- when Step 2c/3 land.  A throw, not an assert: loud in Release too.
template <class T> static tIrrepWF<T>& SameT(const iwf_child_t& c)
{
    if (auto* p=std::get_if<std::unique_ptr<tIrrepWF<T>>>(&c)) return **p;
    throw std::logic_error("tCompositeWF: a child block's scalar differs from the composite face -- "
                           "T-typed operations cannot cross scalars until RealComplexPlan Step 2c/3 land");
}
template <class T> static tIrrepWF<T>& SameT(const iwf_ref_t& c)
{
    if (auto* p=std::get_if<tIrrepWF<T>*>(&c)) return **p;
    throw std::logic_error("tCompositeWF: a child block's scalar differs from the composite face -- "
                           "T-typed operations cannot cross scalars until RealComplexPlan Step 2c/3 land");
}
// The Fock-build dispatch over the child slot (Step 3c-2): a same-scalar child takes the native path; a
// REAL child inside a complex run drives the Hamiltonian's Ham_RealBlock assembly face -- 2b's reserved
// cross arm, now live.  (A complex child inside a real-faced run stays impossible.)
template <class T> static void CalcH(const iwf_child_t& w, tHamiltonian<T>& ham,
                                     const ChargeDensity::tChargeDensity<T>* cd, const tbs_t<T>* bs)
{
    std::visit([&](const auto& p)
    {
        using WT=std::decay_t<decltype(*p)>;
        if constexpr (std::is_same_v<WT,tIrrepWF<T>>)
            p->CalculateH(ham,cd,bs);
        else if constexpr (std::is_same_v<T,dcmplx>)
        {
            auto* rb=dynamic_cast<qchem::Hamiltonian::Ham_RealBlock*>(&ham);
            if (!rb) throw std::logic_error("tCompositeWF: a real child needs the Hamiltonian's "
                                            "real-block assembly face (RealComplexPlan Step 3c-2)");
            p->CalculateH(*rb,cd,bs);
        }
        else
            throw std::logic_error("tCompositeWF: a complex child inside a real-faced run is impossible");
    }, w);
}

// The REAL child's FILL POLICY (Step 3c-3; simplified by R2.21).  A real TRIM block inside a
// complex-faced run fills under a genuine OccupationPolicy<double> built over the RUN'S OWN
// OccupationState -- so the run keeps ONE entropy ledger, ONE fill clock, ONE arming flag, and the real
// block's MOM reference lands in that same ledger under its own scalar.  Before the R2.21 state split
// this had to be a forwarding VIEW whose reference capture THREW ("run with forceComplex"), because the
// references were typed by the RUN; a mixed run with MOM on therefore failed mid-SCF.  Now there is
// nothing to special-case: the Factory rebuilds the run's OWN policy at the other scalar from the same
// OccupationConfig value, so this is literally the policy its complex siblings use, one scalar over.
// (A held direct-min leg carries a default config -- Integer x Bare -- which is exactly what it decides,
// so its real block holds the stored block too.)
static std::unique_ptr<OccupationPolicy<double>> RealBlockPolicy(OccupationPolicy<dcmplx>& run)
{
    return MakeOccupationPolicy<double>(run.Config(), run.State());
}

// The fill dispatch over the child slot (Step 3c-3, the CalcH pattern): a same-scalar child consults
// the run policy natively; a REAL child inside a complex run fills under RealBlockPolicy (R2.21: a real policy, shared state).  Returns
// BY VALUE (every caller merges into its own ledger anyway).
template <class T> static EnergyLevels FillChild(const iwf_ref_t& w, OccupationPolicy<T>& pol, double ne)
{
    return std::visit([&](auto p) -> EnergyLevels
    {
        using WT=std::decay_t<decltype(*p)>;
        if constexpr (std::is_same_v<WT,tIrrepWF<T>>) return p->FillOrbitals(pol,ne);
        else if constexpr (std::is_same_v<T,dcmplx>)
        {
            auto v=RealBlockPolicy(pol);
            return p->FillOrbitals(*v,ne);
        }
        else
            throw std::logic_error("tCompositeWF: a complex child inside a real-faced run is impossible");
    }, w);
}
template <class T> static EnergyLevels FillChildAtMu(const iwf_ref_t& w, OccupationPolicy<T>& pol, double mu)
{
    return std::visit([&](auto p) -> EnergyLevels
    {
        using WT=std::decay_t<decltype(*p)>;
        if constexpr (std::is_same_v<WT,tIrrepWF<T>>) return p->FillOrbitalsAtMu(pol,mu);
        else if constexpr (std::is_same_v<T,dcmplx>)
        {
            auto v=RealBlockPolicy(pol);
            return p->FillOrbitalsAtMu(*v,mu);
        }
        else
            throw std::logic_error("tCompositeWF: a complex child inside a real-faced run is impossible");
    }, w);
}

// The mixed walk (Step 3c-2): the cross-scalar view first (a real TRIM block inside a complex-faced
// set -- non-null only after the 3c-3 factory decision), the same-scalar face otherwise.  A homogeneous
// set takes the second branch for every block, bit-identical to the pre-3c behaviour.
template <class T> void tCompositeWF<T>::MakeIrrepWFs(Spin s)
{
    for (size_t i=0;i<itsBS->GetNumIBS();++i)
        if (const auto* rb=itsBS->GetRealIBS(i)) MakeOneIrrepWF<double>(rb,s);
        else                                     MakeOneIrrepWF<T>((*itsBS)[i],s);
}

template <class T> template <class U> void tCompositeWF<T>::MakeOneIrrepWF(const tobs_t<U>* b, Spin s)
{
    namespace rpt = qchem::report;
    {
        const hmat_t<U>& S = b->Overlap();
        LASolver<U>* lasb=LASolver<U>::Factory(itsBasisOrtho, itsBasisOrthoTol);
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
            rpt::Set("real",       b->IsReal());   // basis-type fact per block (doc/RealComplexPlan.md Step 1)
            rpt::Set("runsReal",   std::is_same_v<U,double>);   // the block's CONSTRUCTED scalar (3c-3 evidence)
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
            for (size_t idx : qchem::PivotedCholeskyDrops<U>(S))
            {
                rpt::Row r("removed");
                rpt::Set("irrep", IrrepLabel(qns));
                rpt::Set("index", (long)idx);
            }

        tSCFIrrepAccelerator<U>* acc=itsAccelerator->Create(lasb,qns,itsEC->GetN(qns));   // §6 typed Create

        std::unique_ptr<tIrrepWF<U>> wf(new tIrrepWF<U>(b,lasb,qns,acc));
        itsQNWFs[qns]=iwf_ref_t(wf.get());
        itsSpinWFs[s].push_back(iwf_ref_t(wf.get()));
        itsIWFs.push_back(iwf_child_t(std::move(wf))); //Do the move last. wf is invalid after the move.
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
    for (auto& w:itsIWFs) CalcH<T>(w,ham,cd,itsBS); //Feed F,D' into all the irrep accelerators.
    // CalculateProjections() has the DIIS side effect (it accumulates the extrapolation history) so it
    // must run every iteration regardless of MOM -- keep the call.  MOM activation is NO LONGER gated on
    // the accelerator engaging (that was the parked molecular heuristic, and NaF's Null accelerator never
    // engages); it is armed by the ranked reservoir fill as soon as a reference occupied subspace exists.
    itsAccelerator->CalculateProjections();
    for (auto& w:itsIWFs) std::visit([](const auto& p){p->DoSCFIteration();}, w);
}

// Iteration-0 seed: build the Fock from \a seed, diagonalize, fill, and hand back the first real density.
template <class T> std::unique_ptr<tDM_CD<T>> tCompositeWF<T>::Init(tHamiltonian<T>& ham,const tChargeDensity<T>* seed,
                                                                    OccupationPolicy<T>& pol, double mergeTol)
{
    DoSCFIteration(ham, seed);
    FillOrbitals(pol, mergeTol);
    return GetChargeDensity();
}

// Build the Fock and have each irrep accelerator compute its (un-taken) step.  Returns true
// only if every irrep produced a geodesic step; false means at least one wants to diagonalize
// (the seed step) -- the caller should fall back to DoSCFIteration().
template <class T> bool tCompositeWF<T>::BuildFockAndComputeSteps(tHamiltonian<T>& ham,const tChargeDensity<T>* cd)
{
    for (auto& w:itsIWFs) CalcH<T>(w,ham,cd,itsBS);
    bool allStepped=true;
    for (auto& w:itsIWFs) allStepped &= std::visit([](const auto& p){return p->ComputeStep();}, w);
    return allStepped;
}

// Move every irrep's orbitals to geodesic fraction t (commit=false for a line-search trial)
// and refill, so GetChargeDensity() reflects the trial/updated orbitals.
template <class T> void tCompositeWF<T>::MoveOrbitals(OccupationPolicy<T>& pol, double t, bool commit,
                                                      double mergeTol)
{
    for (auto& w:itsIWFs) std::visit([&](const auto& p){p->MoveOrbitals(t,commit);}, w);
    FillOrbitals(pol,mergeTol);
}

template <class T> std::unique_ptr<tDM_CD<T>> tCompositeWF<T>::GetChargeDensity(Spin s) const
{
    using qchem::ChargeDensity::tComposite_CD;
    auto i = itsSpinWFs.find(s);
    assert(i!=itsSpinWFs.end());
    // IBZ point group ctor-injected straight from the basis (empty {} for molecules / unfolded crystals = a
    // trivial no-op).  The basis owns the crystal symmetry, so no setter is threaded through the SCF.
    auto cd = std::make_unique<tComposite_CD<T>>(itsBS->GetReciprocalPointOps());
    // Each block BUILDS its own density, typed by ITS scalar; the composite CD's typed Insert
    // overloads (Step 2a) take either alternative -- this seam is already mixed-ready.
    for (auto& w:i->second) std::visit([&](auto p){cd->Insert(p->GetChargeDensity());}, w);
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
    return OrbitalsOf(i->second);

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
        const rvec_t P = std::visit([](const auto& p){return p->GetBasisPopulations();}, w);
        std::ostringstream os; IrrepOf(RefOf(w)).Write(os);   // "s"/"p"/... (+ spin arrow when polarized)
        const std::string lbl=os.str();
        for (std::size_t i=0;i<P.size();++i)
            rows.push_back(rpt::json{{"irrep",lbl},{"index",(long)i},{"pop",P[i]}});
    }
    rpt::EmitAt("basis","usage",rows,rpt::Detail::Verbose);
}



// (SetMOM/SetSmearing/GetEntropyTerm/AdoptMOMReference/ReleaseMOMReference bodies are GONE -- the
// SCFIterator's OccupationPolicy slot owns that configuration and state, V1.11 inc 3.  The iterator's
// AdoptMOMReference loops this WF's GetQNs() against the source's GetOrbitals(irrep), feeding the policy.)

// ONE loop over the RESERVOIR PARTITION (V1.11 increment 2, SCFStrategyPlan §5b): group the blocks by the
// EC's partition, then fill each reservoir by the rule it supports.  This replaced a three-way decision
// tree over three mode bools + two hand-rolled strategy methods.  The partition DEGRADES to per-block
// singletons in two cases:
//  * a policy that HOLDS the stored blocks (the direct minimiser's HeldOccupationPolicy): any cross-block
//    re-decision -- ranked aufbau or a shared μ -- is exactly what a held fill exists to prevent;
//  * a μ-pooling (non-ranked) reservoir at kT=0.  At kT=0 the Fermi count is a STAIRCASE in μ: the
//    requested N generally falls between two steps, and the bisector then returns a μ whose count is
//    simply wrong -- silently, because the guarding assert inside is compiled out in Release.  Measured
//    2026-08-10: the SEED fill of every global-μ metal ran this way and lost charge (Al: Σw·n = 2.25
//    against Ntot=3; a shared-spin Mn sextet: 6 electrons in ↑, 0 in ↓, against Ntot=7 -- and the mixers
//    were then seeded with those).  The seed fill hits it because the occupation RULE is configured in
//    Iterate (SetSmearing) while the seed fill runs in the CONSTRUCTOR -- the two-phase hazard V1.11
//    increment 3 closes (the policy will exist before the WF does).  Filling the prescribed per-block
//    counts is also what a seed SHOULD do: there is no self-consistent spectrum yet to redistribute over.
template <class T> void tCompositeWF<T>::FillOrbitals(OccupationPolicy<T>& pol, double mergeTol)
{
    itsELevels.clear();
    itsSpin_ELevels.clear();
    pol.BeginFill();                      // reset the run's per-fill aggregates (the -TS accumulator)
    const ReservoirPartition p=itsPartition;
    const bool degrade = pol.HoldsStoredBlocks() || (!p.ranksIntegerFill && pol.SmearingkT()<=0.0);
    // Reservoir key: the spatial symmetry's SequenceIndex (0 when the spatial blocks pool) + the spin
    // (0 when the channels pool).  Not an Irrep, and not the sym POINTER: Irrep's ordering dereferences
    // sym, so the "no spatial identity" half of this key has no Irrep spelling.
    std::map<std::pair<size_t,int>,std::vector<iwf_ref_t>> reservoirs;
    for (auto& w : itsIWFs)
    {
        const iwf_ref_t r=RefOf(w);
        const Irrep& irr=IrrepOf(r);
        reservoirs[{ (p.spansSpatial && !degrade) ? 0 : irr.sym->SequenceIndex(),
                     (p.spansSpin    && !degrade) ? 0 : (int)irr.ms+1 }].push_back(r);
    }
    // Snapshot MOM state ONCE at entry: a ranked reservoir activating it below (after its reference
    // capture) must not leak into this same call's LATER reservoirs (each channel's reference is only
    // captured at its own end).
    const bool useMOM = pol.CrossIrrepMOMArmed();
    for (auto& [key,wfs] : reservoirs)
    {
        assert(!wfs.empty());
        if      (!degrade && p.ranksIntegerFill)                    FillReservoirRanked(pol, wfs, mergeTol, useMOM);
        else if (!degrade && (p.spansSpatial || p.spansSpin))       FillReservoirAtSharedMu(pol, wfs, mergeTol);
        else
            for (auto w : wfs)   // per-block prescribed count: atoms; held fills; μ-partitions at kT=0
            {
                EnergyLevels els=FillChild<T>(w, pol, (double)itsEC->GetN(IrrepOf(w)));
                itsELevels.merge(els,mergeTol);
                itsSpin_ELevels[IrrepOf(w).ms].merge(els,mergeTol);
            }
    }
    // Molecular CROSS-irrep MOM (the ranked fill's) is the PARKED path: no molecular test enables MOM, and
    // it has NOT been re-validated against the new delayed reference capture (tIrrepWF captures at
    // MOMStartIter, so useMOM would score against empty references for the first few fills).  The ACTIVE,
    // tested path is the crystal's WITHIN-irrep MOM in tIrrepWF::FillOrbitals.  Kept for the hard cases.
    if (!degrade && p.ranksIntegerFill) pol.ArmCrossIrrepMOM();   // no-op unless the run enables MOM
}

// RANKED integer fill of one reservoir -- the molecular aufbau, one spin channel (a ranked reservoir never
// pools spins; asserted at construction).  Pick which orbitals across the reservoir's blocks are occupied,
// then fill each block with its resulting electron count.  The point-group irrep an occupied MO lands in
// is an OUTPUT of the SCF (unlike an atom's fixed l-occupation).  Two selection rules:
//   * plain aufbau (default): occupy the globally-lowest eigenvalues;
//   * MOM (once armed): occupy the orbitals with the largest overlap onto the previous iteration's
//     occupied subspace, so a near-degenerate cross-irrep pair on the (non-physical) extrapolated Fock
//     cannot flip the configuration.
// Either way the per-irrep references are re-captured at the end for the next iteration's MOM.
template <class T> void tCompositeWF<T>::FillReservoirRanked(OccupationPolicy<T>& pol,
                                                             const std::vector<iwf_ref_t>& wfs, double mergeTol,
                                                             bool useMOM)
{
    struct Slot { double key; iwf_ref_t w; double cap; };
    assert(!wfs.empty());
    const Spin s = IrrepOf(wfs.front()).ms;
    double Nc = (double)itsEC->GetN(IrrepOf(wfs.front()));   // total electrons in this spin channel

    std::map<Irrep,rvec_t> mom;                              // per-irrep MOM scores (empty if no ref), keyed by irrep
    if (useMOM) for (auto w : wfs)
    {
        // Cross-cast the base Orbitals to the policy's OrbitalView (abstract->abstract, the sanctioned
        // direction; TOrbitals implements the view).  T-typed (the policy carries the face's T) -> SameT.
        auto* v=dynamic_cast<const qchem::OrbitalView<T>*>(SameT<T>(w).GetOrbitals());
        assert(v && "ranked reservoir fill: the orbitals must implement OrbitalView");
        mom[IrrepOf(w)]=pol.Scores(IrrepOf(w),*v);
    }

    std::vector<Slot> slots;                                  // every orbital across the channel
    for (auto w : wfs)
    {
        size_t idx=0;
        const rvec_t& sc = mom[IrrepOf(w)];               // empty unless MOM active & referenced
        for (auto o : OrbitalsOf(w)->Iterate())
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
    for (auto w : wfs) ne[IrrepOf(w)]=0.0;
    assert(ne.size()==wfs.size() && "aufbau key collision: irreps must be unique within a spin channel");
    double rem = Nc;
    for (const auto& sl : slots)                              // fill highest-priority first
    {
        if (rem<=0.0) break;
        double take = std::min(sl.cap, rem);
        ne[IrrepOf(sl.w)] += take;
        rem -= take;
    }

    for (auto w : wfs)                                        // occupy + collect energy levels
    {
        EnergyLevels els = SameT<T>(w).FillOrbitals(pol, ne[IrrepOf(w)]);
        itsELevels.merge(els, mergeTol);
        itsSpin_ELevels[s].merge(els, mergeTol);
    }
    // (The per-irrep MOM reference is captured inside tIrrepWF::FillOrbitals above, so no capture here.)
}

// SMEARED fill of one reservoir: ONE chemical potential over the blocks that share the electron reservoir.
//
// A μ is a property of a reservoir, not of a spin channel, and two things can widen the reservoir:
//   * spansSpatial (doc/GPWPlan1.md item 3) -- the Bloch k-blocks share it, so charge sloshes between
//     k-points (a partially-filled band, which per-block fixed occupation cannot represent);
//   * spansSpin -- the two SPIN channels share it, so the MOMENT is an output rather than a constraint.
//     Without it the fill conserves nUp and nDown separately, which means μ↑ ≠ μ↓ and an occupation is
//     monotone in ε only WITHIN a channel: MnO run 29 ended with a ↓ level 27 mHa BELOW an ↑ level and
//     LESS occupied (f↓=0.51 at ε=0.2644 against f↑=0.52 at ε=0.2914).  Nothing is wrong with the
//     arithmetic -- Δμ=27 mHa is simply an unrelieved driving force to move charge ↓→↑ that the constraint
//     holds open, so the converged state is not the free minimum (user, 2026-08-10; §7 step 7).
//
// The BZ weight rides in the level capacity g ONLY when the reservoir actually spans the mesh, because
// that is exactly when GetN is a whole-mesh total rather than a per-block count.  The SAME unweighted
// qchem::Orbitals::FermiLevel bisector serves every case, and one block with w=1 reduces to the per-block
// Fermi exactly.
template <class T> void tCompositeWF<T>::FillReservoirAtSharedMu(OccupationPolicy<T>& pol,
                                                                 const std::vector<iwf_ref_t>& wfs, double mergeTol)
{
    const double kT=pol.SmearingkT();
    assert(kT>0.0 && "a shared-μ fill needs SmearingkT>0: integer fills have nothing to relax with");
    assert(!wfs.empty());
    const bool trace=(bool)std::getenv("GPW_METALTRACE");
    // The reservoir's electron count: each spin CHANNEL present contributes its count once.  (With a
    // shared k-mesh GetN is already the whole-mesh total, so summing per block would over-count by n_k.)
    double Ntot=0.0;
    std::set<Spin> counted;
    for (auto w : wfs) if (counted.insert(IrrepOf(w).ms).second) Ntot+=(double)itsEC->GetN(IrrepOf(w));

    // Aggregate the reservoir's levels: bare energy ε and capacity g (BZ-weighted only across the mesh).
    size_t ntot=0;
    for (auto w : wfs) ntot += OrbitalsOf(w)->GetNumOrbitals();
    rvec_t e(ntot), g(ntot);
    size_t idx=0;
    for (auto w : wfs)
    {
        const double wk=itsPartition.spansSpatial ? IrrepOf(w).sym->GetWeight() : 1.0;  // the weight D carries too
        for (auto o : OrbitalsOf(w)->Iterate())
        {
            e[idx]=o->GetEigenEnergy();
            g[idx]=wk*(double)o->GetDegeneracy();
            ++idx;
        }
    }
    const double mu=qchem::Orbitals::FermiLevel(e,g,Ntot,kT);   // ONE μ for this reservoir

    // Instrument GPW_METALTRACE (cf. GPW_GDMTRACE): μ and the per-block electron count, so the charge
    // sloshing is visible -- between k-points, and between SPIN CHANNELS.  Off by default.
    if (trace) std::cout<<"[metal] μ="<<mu<<" kT="<<kT<<" Ntot="<<Ntot<<"  per-block n:\n";
    double wsum=0.0;
    for (auto w : wfs)                                               // set every block at the shared μ
    {
        EnergyLevels els=FillChildAtMu<T>(w, pol, mu);
        if (trace)
        {
            double nk=0.0; for (auto o:OrbitalsOf(w)->Iterate()) nk+=o->GetOccupation();
            const double wk=itsPartition.spansSpatial ? IrrepOf(w).sym->GetWeight() : 1.0;
            wsum+=wk*nk;
            std::cout<<"[metal]   block="<<IrrepOf(w)<<" w="<<wk<<" n="<<nk<<" (w·n="<<wk*nk<<")\n";
        }
        itsELevels.merge(els,mergeTol);
        itsSpin_ELevels[IrrepOf(w).ms].merge(els,mergeTol);
    }
    if (trace) std::cout<<"[metal] Σ w·n = "<<wsum<<"\n";
}

template class tCompositeWF<double>;
template class tCompositeWF<dcmplx>;

} //namespace
