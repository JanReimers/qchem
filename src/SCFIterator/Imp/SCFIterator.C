// File: SCFIterator/Imp/SCFIterator.C  Partial common implementation for an object that manages SCF convergence.
module;
#include <iostream>
#include <iomanip>
#include <complex>
#include <cassert>
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
import qchem.ChargeDensity.FourierDensity;  // FourierDensity (rho-tilde extraction from the working density)
import qchem.ChargeDensity.FourierMixCD;    // FourierMixCD / KerkerMix (the periodic rho-mixing vehicle)
import qchem.BasisSet.Band_FT_IBS;          // Band_FT_IBS::CreateVxcFitBasisSet (the G-space fit basis)
import qchem.ReciprocalLattice;             // ReciprocalLattice (the Kerker |G| metric)
import qchem.UnitCell;                       // UnitCell (reciprocal cell + volume)

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
    : tSCFIterator(bs,ec,H,acc, ChargeDensity::MakeSeedDensity<T>(seed,bs,st,ec), st, basisOrtho, basisOrthoTol)
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
        if (auto* cell = dynamic_cast<const UnitCell*>(st)) itsKerkerCell = std::make_shared<const UnitCell>(*cell);
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
    size_t idealVirial=itsHamiltonian->IsRelativistic() ? 1 : 2;
    if (ipar.Verbose)
    {
        cout << endl << endl;
        cout << " #           Etotal       " << idealVirial << "+V/K    Δ[F,D]    Δρ    ";
        itsAccelerator->ShowLabels(cout);
        cout << "   relax   Configuration" << endl;
        cout << "                         ";
        cout << "(" << setw(8) << std::scientific << setw(5) << setprecision(0) << ipar.MinVirial  << ")  ";
        cout << "(" << setw(8) << std::scientific << setw(5) << setprecision(0) << ipar.MinΔFD  << ")  ";
        cout << "(" << setw(8) << std::scientific << setw(5) << setprecision(0) << ipar.MinΔρ  << ") ";
        cout << "(" << setw(8) << std::scientific << setw(5) << setprecision(0) << ipar.MinFD  << ") ";
        cout << endl;
        cout << "----------------------------------------------------------------------------------------------------" << endl;
    }

    double ChargeDensityChange=1;
    double  E=0, Eold=0, dE=1e10;
    double FD=0,FDold=1e10,dFD=1e10;
    // double Eoldold=0;
    EnergyBreakdown eb;
    itsConverged=false;
    // Density-face mixer for this run (the density-mixing policy + state; doc/SCFStrategyPlan.md):
    // Kerker ρ̃-mixing when KerkerG0>0 AND the basis/cell/seed are periodic, else linear D-mixing.
    // α=StartingRelaxRo (1.0 default = passthrough -- there is no NullMixer).
    itsMixer = qchem::ChargeDensity::MakeDensityMixer<T>(ipar.StartingRelaxRo, ipar.KerkerG0, ipar.PulayDepth,
                                                         ipar.PulayStart, itsBS, itsKerkerCell.get(), itsCD.get());

    int holeRun=0, momReleases=0;   // 0h MOM-guard state (per run): consecutive hole iterations + releases
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
        const tLoopDriver<T>& driver = itsAccelerator->WantsLineSearch()
                                     ? static_cast<const tLoopDriver<T>&>(itsDirectDriver)
                                     : static_cast<const tLoopDriver<T>&>(itsFixedDriver);
        ChargeDensityChange = driver.Step(lc);
        // cout << "Total charge=" << itsCD->GetTotalCharge() << endl;

        eb=itsHamiltonian->GetTotalEnergy(itsCD.get());
        E=eb.GetTotalEnergy();
        dE=(E-Eold)/fabs(E);
        itsAccelerator->SetEnergy(E); //the ladder gates its hand-off on the energy change
        FD=itsAccelerator->GetError(); //i.e. [F,D]
        dFD=(FD-FDold);
        if (ipar.Verbose) DisplayEnergies(itsIterationCount,eb,itsMixer->GetRelax(),dFD,ChargeDensityChange,idealVirial);
        if (itsObserver) itsObserver({itsIterationCount, E, fabs(E-Eold), FD, ChargeDensityChange});
        // Adaptive [F,D]-keyed policy (LinearMixer only; Kerker/direct-min take the no-op defaults).  The
        // re-damp re-fetches the fresh density + recomputes the energy -- the density LIFECYCLE stays here.
        if (itsMixer->WantsReDamp({E,FD,FDold}))
        {
            SetWorkingCD(cd_t(itsWaveFunction->GetChargeDensity())); //Get new charge density.
            ChargeDensityChange = itsMixer->ReDampMix(itsCD, itsOldCD);
            eb=itsHamiltonian->GetTotalEnergy(itsCD.get());
        }
        itsMixer->UpdateRelax({E,FD,FDold});

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

    return ChargeDensityChange <= ipar.MinΔρ;
}

// One direct-minimization step: build the Fock, compute each accelerator's geodesic step,
// then a backtracking line search along the geodesic (the first fraction from 1 that lowers
// the total energy is accepted) -- no density mixing.  Falls back to a diagonalizing
// iteration in the seed step (before the accelerators have orbitals).
template <class T> typename tSCFIterator<T>::cd_t tSCFIterator<T>::DirectMinStep(double Ecur, double mergeTol)
{
    if (!itsWaveFunction->BuildFockAndComputeSteps(*itsHamiltonian,itsCD.get()))
    {
        itsWaveFunction->DoSCFIteration(*itsHamiltonian,itsCD.get());
        itsWaveFunction->FillOrbitals(mergeTol);
        return cd_t(itsWaveFunction->GetChargeDensity());
    }
    double t=1.0;
    for (int k=0;k<12;k++)
    {
        itsWaveFunction->MoveOrbitals(t,false,mergeTol);                 //trial
        cd_t cdt(itsWaveFunction->GetChargeDensity());                   //std-managed (no freed-address reuse)
        double Et=itsHamiltonian->GetTotalEnergy(cdt.get()).GetTotalEnergy();
        if (Et<Ecur) break;
        t*=0.5;
    }
    itsWaveFunction->MoveOrbitals(t,true,mergeTol);                      //commit at t
    return cd_t(itsWaveFunction->GetChargeDensity());
}


template <class T> void tSCFIterator<T>::DisplayEigen() const
{
    itsWaveFunction->DisplayEigen();
}

// (KerkerSetup/KerkerUpdate/FockDensity moved into qchem.ChargeDensity.DensityMixer -- the density-face
//  seam.  The iterator now builds a tDensityMixer at the top of Iterate and delegates via Mix/FockDensity.)

template <class T> EnergyBreakdown tSCFIterator<T>::GetEnergy() const
{
    return itsHamiltonian->GetTotalEnergy(itsCD.get());
}



template <class T> void tSCFIterator<T>::DisplayEnergies(int i, const EnergyBreakdown& eb, double relax, double dE, double dCD, size_t idealVirial) const
{
    cout.setf(ios::fixed,ios::floatfield);
    cout << setw(3)  << i << " ";
    cout << setw(2+6+12) << setprecision(12) << eb.GetTotalEnergy() << " ";
    cout << setw(8) << std::scientific << setw(8) << setprecision(1) << eb.GetVirial()+idealVirial << " ";
        
    cout << setw(8) << std::scientific << setw(8) << setprecision(1) << dE  << " ";
    cout << setw(8) << std::scientific << setw(7) << setprecision(1) << dCD << " ";
    itsAccelerator->ShowConvergence(cout);
    cout << setw(4) << std::fixed << setw(4) << setprecision(2) << relax << "  ";
    // cout << ConfigString(itsWaveFunction);
    if (ReportBandGap())
    {
        GapInfo g=HomoLumo(itsWaveFunction);
        cout << " εH=";
        if (g.haveHomo) cout << std::fixed << setw(10) << setprecision(5) << g.eHomo; else cout << "     ----";
        cout << " εL=";
        if (g.haveLumo) cout << std::fixed << setw(10) << setprecision(5) << g.eLumo; else cout << "     ----";
        cout << " gap=";
        if (g.haveHomo && g.haveLumo) cout << std::scientific << setw(9) << setprecision(2) << g.gap; else cout << "     ----";
        if (g.metallic) cout << " [partial-occ HOMO]";
        if (g.hole)     cout << " [HOLE: non-aufbau]";
        cout << endl << "        frontier ε(occ): " << FrontierWindow(itsWaveFunction, 2, 4);
    }
    cout << endl;
}

//-------------------------------------------------------------------------------------------------------------------------

template class tSCFIterator<double>;
template class tSCFIterator<dcmplx>;

} //namespace