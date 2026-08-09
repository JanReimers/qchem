// File: SCFIterator/Imp/SCFIterator.C  Partial common implementation for an object that manages SCF convergence.
module;
#include <iostream>
#include <iomanip>
#include <complex>
#include <cassert>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <utility>
#include <cctype>
#include <memory>
#include <type_traits>

module qchem.SCFIterator;
import qchem.SCFParams;
import qchem.SCFAccelerator;

import qchem.WaveFunction;
import qchem.WaveFunction.Factory;

import qchem.Hamiltonian;
import qchem.Hamiltonian.Factory;  // build a DFT sibling for the HF/DHF SAD bootstrap (molecular, double-only)
import qchem.Mesh;                 // qcMesh::MeshParams (defaulted -- a seed-quality mesh for the sibling)
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.ChargeDensity.Seed;   // SeedStrategy / MakeSeedDensity
// The Kerker/Pulay G-space machinery (FourierDensity, FourierMixCD, Band_FT_IBS, ReciprocalLattice)
// is entirely the periodic mixer factory's business -- the iterator only holds a tDensityMixer, so
// those imports are gone.  What survives is one periodicity QUESTION for the cell snapshot handed to
// that factory (CleanupCandidates V1.10b retires even that, by building the mixer above the iterator).
import qchem.ReciprocalLattice;              // isPeriodicCell (the structure-neutral periodicity probe)

import qchem.ElectronConfiguration;
import qchem.Math;

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::ios;


namespace qchem::SCFIterator
{

// The band-gap instrument (doc/GPWPlan §0b″).  Default OFF; the NaF Γ test flips it around Iterate.
bool& ReportBandGap() { static bool on = false; return on; }

// Frontier spectrum of the current orbitals.  The energy levels are energy-ordered (multimap key), so
// ε_HOMO is the highest-energy level still carrying occupation and ε_LUMO the LOWEST empty level -- over
// ALL unoccupied levels, NOT just those above ε_HOMO (doc/GPWPlan 0h: the old above-the-HOMO-index scan
// printed gap=0.67 while a −0.36 Ha virtual sat BELOW the occupied set -- it MASKED the hole that is the
// whole diagnostic).  A NON-AUFBAU configuration (ε_LUMO < ε_HOMO) is flagged as \c hole and the gap goes
// honestly negative.  \c metallic flags a partially-occupied frontier level (occ not near a full shell) --
// a genuine (near-)degenerate crossing where "the gap" is ill-defined; itself the finding for NaF at Γ.
struct GapInfo { double eHomo=0, eLumo=0, gap=0; bool haveHomo=false, haveLumo=false, metallic=false, hole=false; };

template <class T> static GapInfo HomoLumo(const qchem::WaveFunction::tWaveFunction<T>* wf)
{
    const double occTol=1e-6, holeTol=1e-6;
    GapInfo g;
    auto els=wf->GetEnergyLevels();
    // First pass: ε_HOMO = last (highest) level with occ>occTol; note any partial occupation on it.
    for (auto it=els.begin(); it!=els.end(); ++it)
    {
        const auto& lvl=it->second;
        if (lvl.occ>occTol) { g.eHomo=lvl.e; g.haveHomo=true;
                              g.metallic = lvl.occ < (double)lvl.degen - occTol; }   // shell not full => crossing
    }
    // Second pass: ε_LUMO = the FIRST empty level in energy order, wherever it sits.
    for (auto it=els.begin(); it!=els.end(); ++it)
    {
        const auto& lvl=it->second;
        if (lvl.occ<=occTol) { g.eLumo=lvl.e; g.haveLumo=true; break; }
    }
    if (g.haveHomo && g.haveLumo)
    {
        g.gap =g.eLumo-g.eHomo;
        g.hole=g.eLumo < g.eHomo-holeTol;   // an empty level below an occupied one: non-aufbau
    }
    return g;
}

// Frontier-window spectrum: the \a nOcc highest occupied and \a nVirt lowest virtual levels around the
// Fermi edge, "ε(occ)" each, with a "|" marking the gap.  Discriminates the NaF Γ mechanism (doc/GPWPlan
// §0b″): a CLUSTER of near-degenerate low virtuals is the wide-diffuse-band signature (a diffuse Gaussian
// has large inter-cell overlap → a wide band whose Γ minimum sits low and responds giantly to the slosh),
// whereas an ISOLATED LUMO well above the pack is a clean conduction state.  Watching the window per
// iteration shows the diving virtual pull DOWN out of (or through) the pack one step before each spike.
template <class T> static std::string FrontierWindow(const qchem::WaveFunction::tWaveFunction<T>* wf,
                                                      int nOcc, int nVirt)
{
    const double occTol=1e-6;
    std::vector<std::pair<double,double>> occ, virt;   // {energy, occupation}, energy-ordered
    auto els=wf->GetEnergyLevels();
    for (auto it=els.begin(); it!=els.end(); ++it)
        (it->second.occ>occTol ? occ : virt).push_back({it->second.e, it->second.occ});
    std::ostringstream os;
    os.setf(std::ios::fixed,std::ios::floatfield);
    int oFrom=std::max(0,(int)occ.size()-nOcc);
    for (int i=oFrom; i<(int)occ.size(); ++i)
        os << setw(9) << setprecision(4) << occ[i].first << "(" << setprecision(1) << occ[i].second << ") ";
    os << "| ";
    for (int i=0; i<std::min(nVirt,(int)virt.size()); ++i)
        os << setw(9) << setprecision(4) << virt[i].first << "(" << setprecision(1) << virt[i].second << ") ";
    return os.str();
}

// Compact textbook electron configuration of the current orbitals, e.g. (1a₁)²(2a₁)²(1b₂)²...:
// every occupied level (energy-ordered) as (n + Mulliken-label)^occupation, the label lowercased
// with its digits subscripted and the occupation superscripted.  Lets the per-iteration trace
// show which irreps the electrons sit in (and catch any occupation flips under acceleration).
template <class T> static std::string ConfigString(const qchem::WaveFunction::tWaveFunction<T>* wf)
{
    static const char* SUB[]={"₀","₁","₂","₃","₄","₅","₆","₇","₈","₉"};
    static const char* SUP[]={"⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹"};
    auto script=[](int v,const char** t){ std::string d; if(v<=0) return std::string(t[0]);
        while(v>0){ d=std::string(t[v%10])+d; v/=10; } return d; };
    std::string cfg;
    auto els = wf->GetEnergyLevels();
    for (auto it=els.begin(); it!=els.end(); ++it)
    {
        const auto& lvl = it->second;
        int occ=(int)(lvl.occ+0.5);
        if (occ<=0) continue;
        std::string lab;
        for (char c : lvl.qns.sym->GetLabel())
            if (std::isdigit((unsigned char)c)) lab += SUB[c-'0'];
            else                                lab += (char)std::tolower((unsigned char)c);
        cfg += "(" + std::to_string(lvl.qns.n) + lab + ")" + script(occ,SUP);
    }
    return cfg;
}


// The SeedStrategy ctor resolves the strategy into a concrete (heap, owned) density -- nullptr for
// CoreGuess -- then DELEGATES to the explicit-seed ctor below.  \a st (the structure) is consumed only by
// the SAD seeds; bs/st are also forwarded (by the target ctor) for the HF/DHF bootstrap.
template <class T> tSCFIterator<T>::tSCFIterator(const tbs_t<T>* bs, const ElectronConfiguration* ec,ham_t* H,acc_t* acc,ChargeDensity::SeedStrategy seed,const Structure* st,qchem::Ortho basisOrtho,double basisOrthoTol)
    : tSCFIterator(bs,ec,H,acc, ChargeDensity::MakeSeedDensity<T>(seed,bs,st,ec, H&&H->IsPolarized()), st, basisOrtho, basisOrthoTol)
{}

// The explicit-seed ctor (grid-continuation): \a seedDensity is a pre-built density (owned; consumed in
// Init) instead of a strategy enum.  Shared construction body -- the SeedStrategy ctor delegates here.
template <class T> tSCFIterator<T>::tSCFIterator(const tbs_t<T>* bs, const ElectronConfiguration* ec,ham_t* H,acc_t* acc,tChargeDensity<T>* seedDensity,const Structure* st,qchem::Ortho basisOrtho,double basisOrthoTol)
    : itsHamiltonian (H )
    , itsAccelerator (acc)
    , itsWaveFunction(qchem::WaveFunction::Factory(itsHamiltonian,bs,ec,itsAccelerator,basisOrtho,basisOrthoTol) )
    , itsCD          (nullptr)
    , itsOldCD       (nullptr)
    , itsIterationCount(0)
    , itsConverged(false)
{
    assert(itsHamiltonian);
    assert(itsWaveFunction);
    itsBS=bs;   // kept for the optional Kerker G-space fit basis (bs outlives the iterator; safe to hold raw)
    // Deep-copy the periodic cell NOW (st is valid here but comes from a temporary Lattice_3D::GetStructure,
    // so it dangles by Iterate time) -- only for the dcmplx periodic path, so molecular runs pay nothing.
    if constexpr (std::is_same_v<T,dcmplx>)
    {
        // The dcmplx SCF lineage IS the periodic one (there is no complex molecular path), so periodicity
        // is an INVARIANT to assert here -- not a condition to branch on.  Branching on it was the
        // periodic-vs-molecular smell: the `if constexpr` above already answers that question at compile
        // time.  Clone() is polymorphic, so a derived cell (e.g. FCC) is not sliced.
        assert((!st || isPeriodicCell(st)) && "the dcmplx SCF lineage is periodic: Structure must be a UnitCell");
        if (st) itsKerkerCell = st->Clone();
    }
    Initialize(seedDensity, bs, st);   // Init owns the seed transiently (see Initialize)
}


template <class T> void tSCFIterator<T>::Initialize(tChargeDensity<T>* seed, const tbs_t<T>* bs, const Structure* st)
{
    // The seed is only needed to build the iteration-0 Fock; it is NOT a working density (it may be a fit
    // with no density matrix), so own it transiently here rather than storing it as itsOldCD.
    std::unique_ptr<tChargeDensity<T>> seedOwner(seed);
    itsOldCD=nullptr;    //set in the SCF loop; the seed is not a working density
    itsIterationCount=0;
    itsConverged=false;
    itsLineage=std::make_shared<qchem::ChargeDensity::Lineage>();   // one lineage per SCF run (see SetWorkingCD)

    // HF/DHF can't build a Fock from a matrix-FREE seed: their exact-exchange K needs the density MATRIX
    // (Vxc::CalcMatrix asserts the density is a rDM_CD).  So if the seed has no matrix (a SAD fit) and this
    // Hamiltonian RequiresDensityMatrix(), bootstrap: run the seed through a default Dirac-Xalpha (LDA) DFT
    // SIBLING -- built right here from bs + st (H is basis-agnostic and never held them) + a defaulted seed
    // mesh -- for one iteration-0 step to manufacture a real D0, then let the SCF loop continue with the real
    // HF/DHF Hamiltonian seeded by D0.  HF/DHF thus get a SAD-quality start (better than the core guess) for
    // free.  Molecular (double) only: no dcmplx Hamiltonian requires a density matrix today.  Mirrors the
    // core guess (whose seed step also builds the first orbitals from a non-Fock operator).
    if constexpr (std::is_same_v<T,double>)
    {
        namespace H = qchem::Hamiltonian;
        if (seed && itsHamiltonian->RequiresDensityMatrix() && !dynamic_cast<const tDM_CD<T>*>(seed))
        {
            assert(bs && st && "HF/DHF SAD bootstrap needs the orbital basis + structure to build a DFT sibling");
            // The sibling is a NON-relativistic LDA Hamiltonian, so D0 is meaningful only for non-rel HF.  A
            // relativistic DHF basis (large+small components) has no such sibling -> a matrix-free seed is a
            // user error there (use CoreGuess).  This is exactly why Ham_DHF_* report RequiresDensityMatrix().
            assert(!itsHamiltonian->IsRelativistic() &&
                   "DHF cannot seed from a matrix-free (SAD) density: the LDA sibling is non-relativistic -- use CoreGuess");
            const H::Pol pol = itsHamiltonian->IsPolarized() ? H::Pol::Polarized : H::Pol::UnPolarized;
            // Non-owning: the sibling lives only for this Init call, and st outlives it (it is the ctor arg).
            H::st_t stView(st, [](const Structure*){});
            std::unique_ptr<H::rHamiltonian> dftSibling(
                H::Factory(pol, stView, 2.0/3.0, qcMesh::MeshParams{}, bs));   // Dirac exchange (alpha=2/3)
            SetWorkingCD(cd_t(itsWaveFunction->Init(*dftSibling, seed, 0.0001)));  // D0: a real (matrix-backed) density
            assert(itsCD);
            return;
        }
    }
    // Iteration-0: the WaveFunction builds the Fock from the seed (or core operator if null), diagonalizes,
    // fills, and returns the first real (matrix-backed) density.
    SetWorkingCD(cd_t(itsWaveFunction->Init(*itsHamiltonian, seed, 0.0001))); //first real (matrix-backed) density
    assert(itsCD);
}
//
//  Recall that the wavefunction is not owned buy this.
//

template <class T> tSCFIterator<T>::~tSCFIterator()
{
    delete itsHamiltonian;
    delete itsAccelerator;
    delete itsWaveFunction;
    // itsCD / itsOldCD are shared_ptr -- freed automatically.
}

template <class T> bool tSCFIterator<T>::Iterate(const SCFParams& ipar)
{
    assert(itsWaveFunction);
    assert(itsHamiltonian);
    assert(itsCD);
    itsWaveFunction->SetMOM(ipar.UseMOM, ipar.MOMStartIter);   // occupation strategy for this run (SCFParams)
    itsWaveFunction->SetSmearing(ipar.SmearingkT, ipar.MOMSmearPenalty);  // Fermi smearing (0=off) + MOM mask; doc/GPWPlan1.md 4b
    size_t idealVirial=itsHamiltonian->IsRelativistic() ? 1 : 2;
    if (ipar.Verbose) DisplayColumnHeaders(cout, ipar, idealVirial);   // per-system (item 2): base=molecular, Solid overrides

    double ChargeDensityChange=1;
    double  E=0, Eold=0, dE=1e10;
    double FD=0,FDold=1e10,dFD=1e10;
    // double Eoldold=0;
    EnergyBreakdown eb;
    itsConverged=false;
    // Density-face mixer for this run (the density-mixing policy + state; doc/SCFStrategyPlan.md):
    // Kerker ρ̃-mixing when KerkerG0>0 AND the basis/cell/seed are periodic, else linear D-mixing.
    // α=StartingRelaxRo (1.0 default = passthrough -- there is no NullMixer).
    itsMixer = CreateMixer(ipar, itsBS, itsKerkerCell.get(), itsCD.get());

    // The order parameter at the STARTING point (the seed's diagonalized density, "iteration 0") -- the
    // reference every later value is judged against.  A probe that already reads ~0 here means the SEED never
    // carried the order, which is a DIFFERENT bug from the loop losing it.
    if (itsOrderProbe && ipar.Verbose)
        cout << "[order] iteration 0 (seed): " << itsOrderName << " = " << itsOrderProbe(*itsCD) << endl;

    int holeRun=0, momReleases=0;   // 0h MOM-guard state (per run): consecutive hole iterations + releases
    std::string prevConfig;         // previous occupied configuration (the cfg '*' change flag; item 2)
    for (itsIterationCount=1;
        itsIterationCount   <= ipar.NMaxIter && !itsConverged;
        itsIterationCount++)
    {
        // LOOP-FACE (doc/SCFStrategyPlan.md): the accelerator's mode selects the driver -- direct-min
        // (GDM/OT: geodesic line search, no mixing) or fixed-point (diagonalize + density-mix).  Queried
        // every iteration so a ladder tail hand-off flips the loop the moment it switches rungs.  The step
        // BODY is virtual (was a mode `if`); the density LIFECYCLE stays here behind the context callbacks.
        LoopContext<T> lc{ itsHamiltonian, itsWaveFunction, itsMixer.get(), &itsCD, &itsOldCD, ipar.MergeTol, Eold,
                           [this](cd_t x){ itsOldCD=itsCD; SetWorkingCD(std::move(x)); },
                           [this](double e,double tol){ return DirectMinStep(e,tol); } };
        // Direct-min (GDM/OT) OWNS the density update via its geodesic line search, so it must DISABLE the
        // density mixer entirely -- both the per-step Mix() (already: DirectMinDriver never calls it) AND the
        // post-step adaptive re-damp/UpdateRelax below (a LinearMixer's WantsReDamp would otherwise re-mix
        // AFTER a geodesic step, corrupting it).  One query drives the driver choice, the mixer bypass, and
        // the honest ρ_mix="----" display.  CanLineSearch() gates the ACTUALITY: a minimizer that WantsLineSearch
        // but is not yet READY (still seeding, or [F,D] above its FDMax engage threshold) takes a MIXED
        // fixed-point step this iteration instead of the unmixed diagonalize that runs away when ill-conditioned.
        const bool lineSearch = itsAccelerator->WantsLineSearch() && itsAccelerator->CanLineSearch();
        const tLoopDriver<T>& driver = lineSearch
                                     ? static_cast<const tLoopDriver<T>&>(itsDirectDriver)
                                     : static_cast<const tLoopDriver<T>&>(itsFixedDriver);
        ChargeDensityChange = driver.Step(lc);
        // cout << "Total charge=" << itsCD->GetTotalCharge() << endl;

        eb=TotalEnergy(itsCD.get());        // includes the Mermin −TS => E is the free energy A when smearing on
        E=eb.GetTotalEnergy();
        dE=(E-Eold)/fabs(E);
        itsAccelerator->SetEnergy(E); //the ladder gates its hand-off on the energy change
        FD=itsAccelerator->GetError(); //i.e. [F,D]
        dFD=(FD-FDold);
        // The caller's ORDER PARAMETER on THIS iteration's working density (§9 diagnostic metric).  Measured
        // before the display so the trace and the observer see the same number, and unconditionally (not just
        // under Verbose): a headless client watching the observer needs it too.  No probe => no cost.
        const double order = itsOrderProbe ? itsOrderProbe(*itsCD) : 0.0;
        if (ipar.Verbose)
        {
            // Build the full per-iteration trace (all fields, all cheap) and let the per-system display
            // virtual render its honest columns (item 2).  The cfg '*' flags an occupation change vs the
            // previous iteration (blank on iteration 1 -- there is no prior config to differ from).
            std::string config = ConfigString(itsWaveFunction);
            const GapInfo g=HomoLumo(itsWaveFunction);   // frontier spectrum for the gap column
            const double N=itsCD->GetTotalCharge();      // Tr(DS); normalise the grid-charge leak per electron
            IterationTrace tr{ itsIterationCount, eb, itsMixer->GetRelax(),
                               itsMixer->Tag(), itsAccelerator->Tag(), itsAccelerator->Count(),
                               FD, dFD, ChargeDensityChange, dE, idealVirial,
                               itsIterationCount>1 && config!=prevConfig, lineSearch,
                               (N!=0.0 ? eb.GridChargeLost/N : 0.0),
                               g.eHomo, g.eLumo, g.gap, g.haveHomo, g.haveLumo, g.metallic, g.hole,
                               itsOrderProbe ? itsOrderName.c_str() : nullptr, order };
            DisplayColumns(cout, tr);
            prevConfig=std::move(config);
        }
        if (itsObserver) itsObserver({itsIterationCount, E, fabs(E-Eold), FD, ChargeDensityChange, order});
        // Adaptive [F,D]-keyed density-mixing policy (LinearMixer only; Kerker takes the no-op defaults).  The
        // re-damp re-fetches the fresh density + recomputes the energy -- the density LIFECYCLE stays here.
        // SKIPPED under direct-min: GDM/OT own the density update, so no post-step re-mix (see lineSearch above).
        if (!lineSearch && itsMixer->WantsReDamp({E,FD,FDold}))
        {
            SetWorkingCD(cd_t(itsWaveFunction->GetChargeDensity())); //Get new charge density.
            ChargeDensityChange = itsMixer->ReDampMix(*itsCD, *itsOldCD);
            eb=TotalEnergy(itsCD.get());
        }
        if (!lineSearch) itsMixer->UpdateRelax({E,FD,FDold});

        // Eoldold=Eold;
        Eold=E;
        FDold=FD;
        // cout << "ChargeDensityChange    < ipar.MinΔρ " << (ChargeDensityChange < ipar.MinΔρ) << endl;
        // cout << "fabs(dFD)              < ipar.MinΔFD    " << (fabs(dFD)           < ipar.MinΔFD) << endl;
        // cout << "FD                     < ipar.MinFD   " << (FD                 < ipar.MinFD) << endl;
        // cout << "fabs(eb.GetVirial()+2) < ipar.MinVirial  " << (fabs(eb.GetVirial()+2) < ipar.MinVirial) << endl;
        itsConverged=  ChargeDensityChange < ipar.MinΔρ
                 && fabs(dFD)              < ipar.MinΔFD
                 && fabs(dE)               < ipar.MinΔE      // relative total-energy change (default off)
                 && FD                     < ipar.MinFD
                 && fabs(eb.GetVirial()+idealVirial) < ipar.MinVirial
                  ;
        // 0h MOM GUARD (doc/GPWPlan 0h): a MOM reference is TRUSTED, never verified -- a stale/wrong one
        // (the grid-continuation transfer; a reference captured mid-transient) pins an EXCITED state whose
        // signature is a PERSISTENT HOLE: an unoccupied ε sitting below an occupied ε, iteration after
        // iteration (measured: +0.754 Ha on NaF, the hole 0.36 Ha deep).  On 3 consecutive hole iterations,
        // RELEASE the reference (drop + re-arm the delayed-IMOM capture: aufbau fills for the settling
        // window, then a fresh reference from the now-physical occupied set) and VETO this iteration's
        // convergence so the run gets to relax into the recovered state.  Capped at 2 releases per run --
        // a hole that survives them is reported loudly below, never silently.
        if (ipar.UseMOM)
        {
            GapInfo g=HomoLumo(itsWaveFunction);
            holeRun = g.hole ? holeRun+1 : 0;
            if (holeRun>=ipar.Guard.HolePersistence && momReleases<ipar.Guard.MaxReleases)
            {
                std::cerr << "[MOM guard] PERSISTENT HOLE: unoccupied ε=" << g.eLumo << " sits "
                          << (g.eHomo-g.eLumo) << " Ha below occupied ε=" << g.eHomo
                          << " for " << ipar.Guard.HolePersistence
                          << " iterations -- the MOM reference pins a non-aufbau state. "
                          << "Releasing the reference (aufbau + delayed re-capture)." << std::endl;
                itsWaveFunction->ReleaseMOMReference();
                ++momReleases; holeRun=0;
                itsConverged=false;                          // the recovery needs further iterations
            }
        }
        // DisplayEigen();
    }
    // NEVER SILENT: whatever the recipe, a run that ENDS non-aufbau is reported (the honest instrument --
    // the old εH/εL line masked exactly this; doc/GPWPlan 0h).
    if (GapInfo g=HomoLumo(itsWaveFunction); g.hole)
        std::cerr << "[MOM guard] WARNING: run ended NON-AUFBAU (unoccupied ε=" << g.eLumo << " below occupied ε="
                  << g.eHomo << ") -- excited-state energy; check the occupation recipe (MOM reference/smearing)."
                  << std::endl;
//             Etotal       2+V/K    Del(E)  Del(Ro) [F,D]   Nproj    SVMin   Bail      relax
// │ │                      (1e+05)  (1e-02) (2e-05) (2e-06) 
    itsIterationCount--;
    size_t nprec=12,ndigits=log10(-eb.Een)+1,w=1+ndigits+1+nprec;
    nprec-=ndigits;
    if (ipar.Verbose)
    {
        cout << "----------------------------------------------------------------------------------------------------------" << endl;
        // cout << "ndigits=" << ndigits << " nprec=" << nprec << endl;
        cout << "Energy    Breakdown  "
        << "Total: " << std::fixed << setw(w) << setprecision(nprec) << eb.GetTotalEnergy() << "  "
        << "Kinetic: " << std::fixed << setw(w) << setprecision(nprec) << eb.Kinetic << "  "
        << "Potential: " << std::fixed << setw(w) << setprecision(nprec) << eb.GetPotentialEnergy() << endl;
        cout << "Potential Breakdown  "
        << "Een  : " << std::fixed << setw(w) << setprecision(nprec) << eb.Een << "  "
        << "Eee    : " << std::fixed << setw(w) << setprecision(nprec) << eb.Eee << "  "
        << "Eex      : " << std::fixed << setw(w) << setprecision(nprec) << eb.Exc << endl;
        cout << "Virial               V/K  : " << std::fixed << setw(w) << setprecision(11) << eb.GetVirial() << "  ";
        if (eb.Exc!=0.0)
            cout << "Eee/Exc: " << std::fixed << setw(w) << setprecision(nprec) << eb.Eee/eb.Exc << "  " ;
        if (eb.RestMass!=0.0)
            cout << "RestMass : " << std::fixed << setw(w) << setprecision(nprec) << eb.RestMass ;
        cout << endl;
        DisplayEigen();
    }

    // Basis-usage heat map (doc/GPWPlan1 §1): after convergence, record per-function occupation-weighted
    // populations into the run report's basis.usage (Verbose-only console, always in the json).  Self-guards
    // when no run is open, so non-reporting callers pay nothing.
    itsWaveFunction->EmitBasisUsage();

    return ChargeDensityChange <= ipar.MinΔρ;
}

// One direct-minimization step: build the Fock, compute each accelerator's geodesic step,
// then a backtracking line search along the geodesic (the first fraction from 1 that lowers
// the total energy is accepted) -- no density mixing.  Falls back to a diagonalizing
// iteration in the seed step (before the accelerators have orbitals).
template <class T> typename tSCFIterator<T>::cd_t tSCFIterator<T>::DirectMinStep(double Ecur, double mergeTol)
{
    // The STABLE degradation, shared by every non-geodesic exit: drive the Fock from the MIXED density and
    // fold the result back, so a step the geodesic cannot take becomes an ordinary mixed fixed-point step
    // rather than an unmixed diagonalize (which runs away when ill-conditioned).  α=1 (LinearMixer
    // passthrough) => molecular direct-min unchanged.
    // SAFE ONLY WHEN THE ACCELERATOR HOLDS NO LIVE STEP: DoSCFIteration calls NextOrbitals(), and a minimizer
    // with a computed step TAKES it there (GDM: OrbitalsAt(1.0,true)).  Both call sites below satisfy that --
    // one because ComputeStep FAILED, the other because an exhausted RejectStep armed a forced diagonalize.
    auto mixedStep=[&]()
    {
        itsWaveFunction->DoSCFIteration(*itsHamiltonian, itsMixer->FockDensity(*itsCD));
        itsWaveFunction->FillOrbitals(mergeTol);
        cd_t fresh(itsWaveFunction->GetChargeDensity());
        itsMixer->Mix(*fresh, *itsCD);                   // fold ρ_out into ρ_in (itsCD drove this Fock)
        return fresh;
    };
    // Seed / bail: no geodesic this step (normally pre-empted by CanLineSearch() -> FixedPointDriver, but
    // reachable if [F,D] jumps above FDMax between that check and the fresh Fock).
    if (!itsWaveFunction->BuildFockAndComputeSteps(*itsHamiltonian,itsCD.get())) return mixedStep();

    const bool trace=(bool)std::getenv("GPW_GDMTRACE");
    double prevBest=1e300;   // best(E_t) of the PREVIOUS rejected attempt -- the scale-invariance probe
    for (;;)
    {
        double t=1.0, Et=0, best=1e300; int k=0; bool found=false;
        for (;k<12;k++)
        {
            itsWaveFunction->MoveOrbitals(t,false,mergeTol,/*holdBlock*/true);   //trial, ON the geodesic
            cd_t cdt(itsWaveFunction->GetChargeDensity());                //std-managed (no freed-address reuse)
            // Minimize the FREE energy A=E−TS under smearing (MoveOrbitals refilled, so GetEntropyTerm is
            // current); GetEntropyTerm()=0 with no smearing => molecular direct-min unchanged.  GPWPlan1 4b.
            Et=itsHamiltonian->GetTotalEnergy(cdt.get()).GetTotalEnergy()+itsWaveFunction->GetEntropyTerm();
            best=std::min(best,Et);
            if (Et<Ecur) { found=true; break; }
            t*=0.5;
        }
        // DIAGNOSTIC (env GPW_GDMTRACE): did the geodesic line search find a DESCENT (some t lowers Ecur) or
        // fall through all 12 backtracks?  A variational E with a correct geodesic MUST descend for small
        // enough t; a FALLBACK whose best(Et−Ecur) is at the E noise floor is that floor, while one that is
        // ORDERS above it means E(t) does not even tend to E(0) -- a discontinuity, i.e. the occupation
        // jumped branch between trial points (MnO: +1.45e+01 Ha at t=2.4e-4).  The two are worth telling
        // apart by eye, which is why the number is printed and not just the verdict.
        if (trace)
            std::cout << "[gdm] " << (found ? "DESCENT " : "FALLBACK") << " t=" << std::scientific << setprecision(2)
                      << t << " k=" << k << "  best(Et-Ecur)=" << (best-Ecur) << std::defaultfloat << std::endl;
        if (found)
        {
            itsWaveFunction->MoveOrbitals(t,true,mergeTol,/*holdBlock*/true);    //commit at t
            return cd_t(itsWaveFunction->GetChargeDensity());
        }
        // NO DESCENT AT ANY BACKTRACK -- do NOT take the step.  Committing here (the behaviour through
        // 2026-08-09) made the method non-monotone BY CONSTRUCTION: `found` fed only the trace above, so a
        // failed search still moved the orbitals, and on MnO that committed a +14.5 Ha step whose blast
        // radius was five further iterations.  Instead, hand the rejection back to the accelerator (the
        // bail-out contract, tSCFAccelerator::RejectStep) and let IT decide whether a degraded direction is
        // worth another look.  Note the trials have already MOVED the orbitals -- MoveOrbitals mutates even
        // when commit=false, only the accelerator's history is gated on commit -- so a rewind to t=0 is
        // required before either exit; the geodesic at θ=0 is exactly the pre-step point.
        itsWaveFunction->MoveOrbitals(0.0,false,mergeTol,/*holdBlock*/true); //rewind the trials
        // SCALE-INVARIANCE SHORT-CIRCUIT.  A step that merely overshoots gets BETTER as the trust radius
        // shrinks; a failure whose best(E_t) does not move is not about step LENGTH at all, and no further
        // backoff can help (measured on MnO before the held fill: seven rejections, best(E_t−E_cur) = +1.45e+01
        // every time, identical to three figures across six shrinks).  Each retry costs a full 12-point
        // backtrack -- ~12 Fock builds -- so stop asking as soon as the answer stops changing.
        const bool improved = best < prevBest - 1e-10*std::fabs(prevBest);
        const bool stalled  = (prevBest<1e299) && !improved;
        prevBest = std::min(prevBest,best);
        if (stalled)
        {
            if (trace) std::cout << "[gdm] STALLED (best unchanged under trust backoff -- not a step-length "
                                    "failure); abandoning the geodesic" << std::endl;
            itsAccelerator->RejectStep();   // still tell the accelerator, so it arms its forced diagonalize
            while (itsAccelerator->RejectStep()) {}   // ...and drive it to exhaustion, cheaply (no line search)
            return mixedStep();
        }
        if (!itsAccelerator->RejectStep())
        {
            // EXHAUSTED.  By contract the accelerator has armed a diagonalizing step, so NextOrbitals()
            // will NOT re-take the rejected geodesic and the mixed fallback is now safe.
            if (trace) std::cout << "[gdm] EXHAUSTED -> mixed fallback step" << std::endl;
            return mixedStep();
        }
        if (trace) std::cout << "[gdm] RETRY (direction rebuilt, trust shrunk)" << std::endl;
    }
}


template <class T> void tSCFIterator<T>::DisplayEigen() const
{
    itsWaveFunction->DisplayEigen();
}

// (KerkerSetup/KerkerUpdate/FockDensity moved into qchem.ChargeDensity.DensityMixer -- the density-face
//  seam.  The iterator now builds a tDensityMixer at the top of Iterate and delegates via Mix/FockDensity.)

template <class T> EnergyBreakdown tSCFIterator<T>::GetEnergy() const
{
    return TotalEnergy(itsCD.get());
}

// The Hamiltonian's total energy for cd, plus the wavefunction's Mermin −TS stamped in (0 with no
// smearing) so GetTotalEnergy() reads the free energy A=E−TS.  doc/GPWPlan1.md 4b.
template <class T> EnergyBreakdown tSCFIterator<T>::TotalEnergy(const tDM_CD<T>* cd) const
{
    EnergyBreakdown eb=itsHamiltonian->GetTotalEnergy(cd);
    eb.MinusTS = itsWaveFunction->GetEntropyTerm();
    return eb;
}



//====================================================================================================
//  PER-SYSTEM ITERATION DISPLAY (doc/GPWPlan1.md item 2)
//
//  The base layout is MOLECULAR ({#, Etotal, [F,D], Δ[F,D], Δρ, 2+V/K, ρ_mix, accel, cfg} + the optional
//  ReportBandGap() frontier block); SolidSCFIterator overrides both virtuals for the solid/PP layout.  The
//  shared column WRITERS below (WriteRowPrefix / WriteMixAccelCfg / WriteGapColumn) are reused by both, so
//  the two layouts differ only in the physics columns they insert (virial vs ΔE+gap).
//====================================================================================================

// --- Column-width constants: the header labels and the value rows share these, so they line up exactly
//     regardless of the UTF-8 label glyphs (Δ, ρ, ε are multi-byte but one display column). --------------
namespace { enum : int { W_ITER=3, W_E=20, W_FD=10, W_DELTA=10, W_RHO=9, W_VIR=10, W_LOST=10, W_MIX=8, W_ACC=7, W_CFG=3, W_GAP=9, W_ORD=12 }; }

// Display width of a UTF-8 string (counts leading bytes only, so Δ/ρ/ε each count as one column).
static size_t VisWidth(const std::string& s)
{ size_t n=0; for (unsigned char c : s) if ((c&0xC0)!=0x80) ++n; return n; }
// Right-/left-justify a label to a given DISPLAY width (setw counts bytes, which mis-pads UTF-8 labels).
static std::string PadR(const std::string& s, int w){ int v=(int)VisWidth(s); return v>=w ? s : std::string(w-v,' ')+s; }
static std::string PadL(const std::string& s, int w){ int v=(int)VisWidth(s); return v>=w ? s : s+std::string(w-v,' '); }

// One convergence threshold as it appears in the header sub-line, e.g. "(1e-07)".
static std::string Thresh(double t)
{ std::ostringstream os; os << "(" << std::scientific << setprecision(0) << t << ")"; return os.str(); }

// Header cells matching the row layouts (same widths + separators as WriteRowPrefix / WriteMixAccelCfg).
static void WriteHeadPrefix(std::ostream& os)          // " #   Etotal   [F,D]"
{ os << PadR("#",W_ITER) << " " << PadR("Etotal",W_E) << " " << PadR("[F,D]",W_FD) << " "; }
static void WriteHeadMixAccelCfg(std::ostream& os)     // "ρ_mix  accel  cfg"
{ os << " " << PadL("ρ_mix",W_MIX) << " " << PadL("accel",W_ACC) << " " << PadR("cfg",W_CFG); }
// The threshold sub-line's leading pad (under #, Etotal), so each (thresh) sits under its column.
static void WriteThreshLead(std::ostream& os)
{ os << PadR("",W_ITER) << " " << PadR("",W_E) << " "; }

// Row prefix shared by every system: # , Etotal (fixed, 12 dp), [F,D] (the accelerator's orbital gradient;
// "----" for a Null accelerator, which has no [F,D] -- e.g. a pure density-mixing run).
template <class T> void tSCFIterator<T>::WriteRowPrefix(std::ostream& os, const IterationTrace& tr) const
{
    os.setf(ios::fixed,ios::floatfield);
    os << setw(W_ITER) << tr.iteration << " ";
    os << setw(W_E) << setprecision(12) << tr.eb.GetTotalEnergy() << " ";
    if (std::string(tr.accelTag)=="Null") os << PadR("----",W_FD) << " ";
    else os << std::scientific << setw(W_FD) << setprecision(2) << tr.FD << " ";
}

// The three columns common to the END of every row: ρ_mix (Lin/Ker/Pul + α), accel (tag[:count]), cfg ('*'
// when the occupied configuration changed).  Compact, fixed-width, self-identifying (item 2).
template <class T> void tSCFIterator<T>::WriteMixAccelCfg(std::ostream& os, const IterationTrace& tr) const
{
    std::ostringstream mix;                                     // ρ_mix: "Lin 1.00" -- or "----" under direct-min
    if (tr.lineSearch) mix << "----";                           //   (GDM/OT own the density update; NO mixing)
    else               mix << tr.mixTag << " " << std::fixed << setprecision(2) << tr.relax;
    std::ostringstream acc; acc << tr.accelTag; if (tr.accelCount>0) acc << ":" << tr.accelCount;    // "DIIS:3"
    os << " " << PadL(mix.str(),W_MIX) << " " << PadL(acc.str(),W_ACC)
       << " " << PadR(std::string(1, tr.configChanged ? '*' : ' '), W_CFG);
}

// The frontier gap column: ε_LUMO−ε_HOMO (honestly negative for a hole), with a compact flag -- 'h' a
// non-aufbau HOLE, 'm' a partial-occupancy (near-degenerate/metallic) frontier -- else blank.
template <class T> void tSCFIterator<T>::WriteGapColumn(std::ostream& os, const IterationTrace& tr) const
{
    os << " ";
    if (tr.haveHomo && tr.haveLumo) os << std::scientific << setw(W_GAP) << setprecision(2) << tr.gap;
    else                            os << PadR("----",W_GAP);
    os << ' ' << (tr.hole ? 'h' : tr.metallic ? 'm' : ' ');
}

// The ORDER-PARAMETER column (SetOrderParameter; §9): the caller's named scalar on this iteration's density,
// signed and in FIXED notation -- the point of the column is to watch a value DECAY toward zero (or flip
// sign), which scientific notation hides behind a shrinking exponent.  Absent probe => absent column.
template <class T> void tSCFIterator<T>::WriteOrderColumn(std::ostream& os, const IterationTrace& tr) const
{
    if (!tr.orderName) return;
    std::ostringstream v; v << std::fixed << setprecision(6) << std::showpos << tr.order;
    os << " " << PadR(v.str(),W_ORD);
}
template <class T> void tSCFIterator<T>::WriteHeadOrder(std::ostream& os) const
{
    if (itsOrderName.empty() || !itsOrderProbe) return;
    os << " " << PadR(itsOrderName,W_ORD);
}

// --- Base default = the MOLECULAR layout (atoms/molecules; MolecularSCFIterator inherits it unchanged) ---
template <class T>
void tSCFIterator<T>::DisplayColumnHeaders(std::ostream& os, const SCFParams& ipar, size_t idealVirial) const
{
    os << endl << endl;
    WriteHeadPrefix(os);
    os << PadR("Δ[F,D]",W_DELTA) << " " << PadR("Δρ",W_RHO) << " "
       << PadR(std::to_string(idealVirial)+"+V/K",W_VIR) << " ";
    WriteHeadMixAccelCfg(os);
    WriteHeadOrder(os);
    os << endl;
    WriteThreshLead(os);
    os << PadR(Thresh(ipar.MinFD),W_FD) << " " << PadR(Thresh(ipar.MinΔFD),W_DELTA) << " "
       << PadR(Thresh(ipar.MinΔρ),W_RHO) << " " << PadR(Thresh(ipar.MinVirial),W_VIR) << endl;
    os << std::string(100,'-') << endl;
}

template <class T> void tSCFIterator<T>::DisplayColumns(std::ostream& os, const IterationTrace& tr) const
{
    WriteRowPrefix(os, tr);
    os << std::scientific << setw(W_DELTA) << setprecision(2) << tr.dFD  << " ";
    os << std::scientific << setw(W_RHO)   << setprecision(2) << tr.dRho << " ";
    os << std::scientific << setw(W_VIR)   << setprecision(2) << (tr.eb.GetVirial()+(double)tr.idealVirial) << " ";
    WriteMixAccelCfg(os, tr);
    WriteOrderColumn(os, tr);
    // Optional frontier diagnostic (unchanged instrument; solids show it as a permanent column instead).
    if (ReportBandGap())
    {
        WriteGapColumn(os, tr);
        os << endl << "        frontier ε(occ): " << FrontierWindow(itsWaveFunction, 2, 4);
    }
    os << endl;
}

//-------------------------------------------------------------------------------------------------------------------------

// --- Density mixing, per system type (doc/CleanupCandidates.md V1.10b) ---
// The base answer is the structure-neutral LINEAR D-mixer.  It ignores the basis/cell/seed: a molecular
// run has no G-space rho to mix.  SolidSCFIterator overrides this with the periodic mixer.
template <class T> std::unique_ptr<qchem::ChargeDensity::tDensityMixer<T>>
tSCFIterator<T>::CreateMixer(const SCFParams& ipar, const tbs_t<T>*, const Structure*, const tDM_CD<T>*) const
{
    return qchem::ChargeDensity::MakeLinearMixer<T>(ipar.StartingRelaxRo);
}

template class tSCFIterator<double>;
template class tSCFIterator<dcmplx>;

// --- The SOLID / PP layout: ΔE gates (non-variational collocation), NO virial, gap is a permanent column ---
void SolidSCFIterator::DisplayColumnHeaders(std::ostream& os, const SCFParams& ipar, size_t /*idealVirial*/) const
{
    os << endl << endl;
    WriteHeadPrefix(os);
    os << PadR("ΔE/E",W_DELTA) << " " << PadR("Δρ",W_RHO) << " "     // ΔE column is RELATIVE (dE/|E|)
       << PadR("ρ_lost/N",W_LOST) << " ";                            // grid-charge leak per electron (health)
    WriteHeadMixAccelCfg(os);
    os << " " << PadR("gap",W_GAP) << "  ";   // 2 trailing: WriteGapColumn's row cell carries the h/m flag
    WriteHeadOrder(os);
    os << endl;
    WriteThreshLead(os);
    os << PadR(Thresh(ipar.MinFD),W_FD) << " " << PadR(Thresh(ipar.MinΔE),W_DELTA) << " "
       << PadR(Thresh(ipar.MinΔρ),W_RHO) << endl;
    os << std::string(100,'-') << endl;
}

void SolidSCFIterator::DisplayColumns(std::ostream& os, const IterationTrace& tr) const
{
    WriteRowPrefix(os, tr);
    os << std::scientific << setw(W_DELTA) << setprecision(2) << tr.dE   << " ";
    os << std::scientific << setw(W_RHO)   << setprecision(2) << tr.dRho << " ";
    os << std::scientific << setw(W_LOST)  << setprecision(2) << tr.gridLostRel << " ";  // ρ_lost/N
    WriteMixAccelCfg(os, tr);
    WriteGapColumn(os, tr);   // solids always show the gap (folds the former ReportBandGap() instrument)
    WriteOrderColumn(os, tr); // and the caller's order parameter, when one is probed
    os << endl;
}

// The solid mixer.  The branch here is on the KNOB -- did the caller ASK for G-space mixing? -- not on the
// geometry: being periodic is what this class IS, so MakePeriodicMixer takes the Band_FT_IBS basis, the
// UnitCell and the FourierDensity seed as preconditions rather than probing for them.
std::unique_ptr<qchem::ChargeDensity::tDensityMixer<dcmplx>>
SolidSCFIterator::CreateMixer(const SCFParams& ipar, const tbs_t<dcmplx>* bs, const Structure* cell,
                              const tDM_CD<dcmplx>* seed) const
{
    if (ipar.KerkerG0>0.0 || ipar.PulayDepth>0)
        return qchem::ChargeDensity::MakePeriodicMixer(ipar.StartingRelaxRo, ipar.KerkerG0, ipar.PulayDepth,
                                                       ipar.PulayStart, bs, cell, seed);
    return qchem::ChargeDensity::MakeLinearMixer<dcmplx>(ipar.StartingRelaxRo);
}

} //namespace