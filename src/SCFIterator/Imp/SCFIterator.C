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
// ε_HOMO is the highest-energy level still carrying occupation and ε_LUMO the first empty level above it.
// \c metallic flags a partially-occupied frontier level (occ not near a full shell) -- a genuine
// (near-)degenerate crossing where "the gap" is ill-defined; that is itself the finding for NaF at Γ.
struct GapInfo { double eHomo=0, eLumo=0, gap=0; bool haveHomo=false, haveLumo=false, metallic=false; };

template <class T> static GapInfo HomoLumo(const qchem::WaveFunction::tWaveFunction<T>* wf)
{
    const double occTol=1e-6;
    GapInfo g;
    auto els=wf->GetEnergyLevels();
    // First pass: ε_HOMO = last (highest) level with occ>occTol; note any partial occupation on it.
    for (auto it=els.begin(); it!=els.end(); ++it)
    {
        const auto& lvl=it->second;
        if (lvl.occ>occTol) { g.eHomo=lvl.e; g.haveHomo=true;
                              g.metallic = lvl.occ < (double)lvl.degen - occTol; }   // shell not full => crossing
    }
    // Second pass: ε_LUMO = first empty level strictly above ε_HOMO (fall through to the very first empty
    // level if there is no occupation at all, e.g. a bare seed).
    for (auto it=els.begin(); it!=els.end(); ++it)
    {
        const auto& lvl=it->second;
        if (lvl.occ<=occTol && (!g.haveHomo || lvl.e>g.eHomo))
            { g.eLumo=lvl.e; g.haveLumo=true; break; }
    }
    if (g.haveHomo && g.haveLumo) g.gap=g.eLumo-g.eHomo;
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


template <class T> tSCFIterator<T>::tSCFIterator(const tbs_t<T>* bs, const ElectronConfiguration* ec,ham_t* H,acc_t* acc,ChargeDensity::SeedStrategy seed,const Structure* st,qchem::Ortho basisOrtho,double basisOrthoTol)
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
    // Resolve the seed strategy into a concrete (heap, owned) density -- nullptr for CoreGuess.  \a st
    // (the structure) is consumed only by the SAD seeds; bs/st are also forwarded for the HF/DHF bootstrap.
    Initialize(ChargeDensity::MakeSeedDensity<T>(seed,bs,st,ec), bs, st);
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
    double relax=ipar.StartingRelaxRo;
    double relMax=1.0;
    EnergyBreakdown eb;
    itsConverged=false;
    if (ipar.KerkerG0>0.0) KerkerSetup(ipar.KerkerG0);   // periodic rho-mixing (else linear D-mixing, unchanged)

    for (itsIterationCount=1;
        itsIterationCount   <= ipar.NMaxIter && !itsConverged;
        itsIterationCount++)
    {
        // The accelerator decides the loop mode: a direct minimizer (GDM, or a ladder that has
        // handed off to one near convergence) runs the geodesic line search; everything else
        // runs the classic diagonalize+mix fixed-point step.  Queried every iteration so a
        // ladder tail hand-off flips the loop the moment it switches rungs.
        itsDirectMin = itsAccelerator->WantsLineSearch();
        if (itsDirectMin)
        {
            // GDM owns the loop: a geodesic line search drives the energy down directly, with
            // NO density mixing (the line search guarantees descent).
            itsOldCD=itsCD;
            SetWorkingCD(DirectMinStep(Eold,ipar.MergeTol));
            ChargeDensityChange = itsCD->GetChangeFrom(*itsOldCD)/itsCD->GetTotalCharge();
        }
        else
        {
            // Fock from the Kerker-mixed rho-tilde when active (FockDensity()), else from the working D.
            itsWaveFunction->DoSCFIteration(*itsHamiltonian,FockDensity()); //eigen orbitals from the Hamiltonian
            itsWaveFunction->FillOrbitals(ipar.MergeTol);

            itsOldCD=itsCD;
            SetWorkingCD(cd_t(itsWaveFunction->GetChargeDensity())); //Get new charge density.
            if (itsMixedRho)
                ChargeDensityChange = KerkerUpdate(relax);              //rho-mixing: gate on ‖ρ̃_out−ρ̃_in‖, fold in
            else
            {
                ChargeDensityChange = itsCD->GetChangeFrom(*itsOldCD)/itsCD->GetTotalCharge(); //relative MaxAbs change
                if (ChargeDensityChange<1e-5) relMax=0.5;
                itsCD->MixIn(*itsOldCD,1.0-relax);                       //relaxation (linear D-mixing).
            }
        }
        // cout << "Total charge=" << itsCD->GetTotalCharge() << endl;

        eb=itsHamiltonian->GetTotalEnergy(itsCD.get());
        E=eb.GetTotalEnergy();
        dE=(E-Eold)/fabs(E);
        itsAccelerator->SetEnergy(E); //the ladder gates its hand-off on the energy change
        FD=itsAccelerator->GetError(); //i.e. [F,D]
        dFD=(FD-FDold);
        if (ipar.Verbose) DisplayEnergies(itsIterationCount,eb,relax,dFD,ChargeDensityChange,idealVirial);
        if (itsObserver) itsObserver({itsIterationCount, E, fabs(E-Eold), FD, ChargeDensityChange});
        if (!itsMixedRho && FD>FDold && fabs(dFD)>1e-9)   // [F,D]-keyed re-damping is a linear-D-mixing tweak; Kerker skips it
        {
            SetWorkingCD(cd_t(itsWaveFunction->GetChargeDensity())); //Get new charge density.
            ChargeDensityChange = itsCD->GetChangeFrom(*itsOldCD); //Get MaxAbs of change.
            itsCD->MixIn(*itsOldCD,1.0-relax/4.0);
            eb=itsHamiltonian->GetTotalEnergy(itsCD.get());
            relax*=0.8;
        }
        // if (E<Eold && Eold>Eoldold) relax*=0.5;
        if (!itsMixedRho && FD<FDold ) relax*=1.5;   // [F,D]-keyed relax growth is a D-mixing tweak; Kerker holds α fixed
        // if (E>Eold && Eold>Eoldold) relax*=1.5;
        if (relax>relMax) relax=relMax;

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
        // DisplayEigen();
    }
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

// KERKER rho-mixing setup (dcmplx periodic path only).  Build the G-space fit basis from the orbital basis's
// first block, the reciprocal lattice + cell volume from the periodic Structure, and the INITIAL mixed density
// = the seed density's rho-tilde.  A no-op (leaves itsMixedRho null) for the real/molecular path or when the
// basis/structure are not periodic -- so linear D-mixing (the default) is used unchanged.
template <class T> void tSCFIterator<T>::KerkerSetup(double G0)
{
    if constexpr (std::is_same_v<T,dcmplx>)
    {
        namespace CD = qchem::ChargeDensity;
        auto* ftb = itsBS ? dynamic_cast<const BasisSet::Band_FT_IBS*>((*itsBS)[0]) : nullptr;
        auto* cell= dynamic_cast<const UnitCell*>(itsKerkerCell.get());
        auto* fd  = dynamic_cast<const CD::FourierDensity*>(itsCD.get());
        if (!ftb || !cell || !fd)   // not a periodic G-space SCF -> report LOUDLY (Release has no asserts) + bail
        {
            std::cerr << "[Kerker] DISABLED: KerkerG0>0 needs a periodic Band_FT_IBS basis + UnitCell + "
                      << "FourierDensity (got ftb=" << (void*)ftb << " cell=" << (void*)cell
                      << " fd=" << (void*)fd << ") -- falling back to linear D-mixing." << std::endl;
            return;
        }
        itsKerkerG0 = G0;
        itsKerkerFit.reset(ftb->CreateVxcFitBasisSet(cell, qcMesh::MeshParams{}));  // matches the Ham's grid
        ReciprocalLattice recip(cell->MakeReciprocalCell());
        auto rho0 = fd->GetFourierDensity(*itsKerkerFit);
        itsMixedRho = std::make_shared<CD::FourierMixCD>(rho0, recip, itsCD->GetTotalCharge());
        std::cerr << "[Kerker] ENABLED: G0=" << G0 << ", rho-mixing on " << itsMixedRho->GetTotalCharge()
                  << " electrons (" << rho0.size() << " G-vectors)." << std::endl;
    }
}

// One Kerker update: re-collocate the freshly diagonalized density's rho-tilde and fold it into the mixed
// density by the preconditioner rho_mix = rho_in + relax*G^2/(G^2+G0^2)*(rho_out - rho_in) (charge-conserving:
// G=0 is never mixed).  The result drives the NEXT Fock via FockDensity().
template <class T> double tSCFIterator<T>::KerkerUpdate(double relax)
{
    if constexpr (std::is_same_v<T,dcmplx>)
    {
        namespace CD = qchem::ChargeDensity;
        auto* fd    = dynamic_cast<const CD::FourierDensity*>(itsCD.get());
        auto* mixed = dynamic_cast<const CD::FourierMixCD*>(itsMixedRho.get());
        assert(fd && mixed);
        const ΔG_Map  rho_out = fd->GetFourierDensity(*itsKerkerFit);
        const ΔG_Map& rho_in  = mixed->RhoTilde();
        // SCF residual ‖ρ̃_out − ρ̃_in‖_∞: how far the fed density is from self-consistency (0 at the fixed point).
        // This is the RIGHT convergence gate for ρ-mixing -- NOT the D_out change, which shrinks with the damped
        // Kerker step even when ρ is not yet self-consistent (would converge early to a wrong energy).
        double resid = 0.0;
        for (const auto& [dm, ro] : rho_out)
        {
            auto it = rho_in.find(dm);
            resid = std::max(resid, std::abs(dcmplx(ro) - (it!=rho_in.end() ? dcmplx(it->second) : dcmplx(0.0))));
        }
        for (const auto& [dm, ri] : rho_in)
            if (rho_out.find(dm)==rho_out.end()) resid = std::max(resid, std::abs(dcmplx(ri)));
        itsMixedRho.reset(CD::FourierMixCD::KerkerMix(*mixed, rho_out, relax, itsKerkerG0));
        return resid;
    }
    return 0.0;
}

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
        cout << endl << "        frontier ε(occ): " << FrontierWindow(itsWaveFunction, 2, 4);
    }
    cout << endl;
}

//-------------------------------------------------------------------------------------------------------------------------

template class tSCFIterator<double>;
template class tSCFIterator<dcmplx>;

} //namespace