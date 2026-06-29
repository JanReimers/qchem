// File: SCFIterator/Imp/SCFIterator.C  Partial common implementation for an object that manages SCF convergence.
module;
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <string>
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

import qchem.ElectronConfiguration;
import qchem.Math;

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::ios;


namespace qchem::SCFIterator
{

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


template <class T> tSCFIterator<T>::tSCFIterator(const tbs_t<T>* bs, const ElectronConfiguration* ec,ham_t* H,acc_t* acc,ChargeDensity::SeedStrategy seed,const Structure* st)
    : itsHamiltonian (H )
    , itsAccelerator (acc)
    , itsWaveFunction(qchem::WaveFunction::Factory(itsHamiltonian,bs,ec,itsAccelerator) )
    , itsCD          (nullptr)
    , itsOldCD       (nullptr)
    , itsIterationCount(0)
    , itsConverged(false)
{
    assert(itsHamiltonian);
    assert(itsWaveFunction);
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

    // HF/DHF can't build a Fock from a matrix-FREE seed: their exact-exchange K needs the density MATRIX
    // (Vxc::CalcMatrix asserts the density is a DM_CD).  So if the seed has no matrix (a SAD fit) and this
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
            std::unique_ptr<H::Hamiltonian> dftSibling(
                H::Factory(pol, stView, 2.0/3.0, qcMesh::MeshParams{}, bs));   // Dirac exchange (alpha=2/3)
            itsCD=cd_t(itsWaveFunction->Init(*dftSibling, seed, 0.0001));      // D0: a real (matrix-backed) density
            assert(itsCD);
            return;
        }
    }
    // Iteration-0: the WaveFunction builds the Fock from the seed (or core operator if null), diagonalizes,
    // fills, and returns the first real (matrix-backed) density.
    itsCD=cd_t(itsWaveFunction->Init(*itsHamiltonian, seed, 0.0001)); //first real (matrix-backed) density, std-managed
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
            itsCD=DirectMinStep(Eold,ipar.MergeTol);
            ChargeDensityChange = itsCD->GetChangeFrom(*itsOldCD)/itsCD->GetTotalCharge();
        }
        else
        {
            itsWaveFunction->DoSCFIteration(*itsHamiltonian,itsCD.get()); //Just gets a set of eigen orbitals from the Hamiltonian
            itsWaveFunction->FillOrbitals(ipar.MergeTol);

            itsOldCD=itsCD;
            itsCD=cd_t(itsWaveFunction->GetChargeDensity()); //Get new charge density.
            ChargeDensityChange = itsCD->GetChangeFrom(*itsOldCD)/itsCD->GetTotalCharge(); //Get relative MaxAbs of change.
            if (ChargeDensityChange<1e-5) relMax=0.5;
            itsCD->MixIn(*itsOldCD,1.0-relax);                           //relaxation.
        }
        // cout << "Total charge=" << itsCD->GetTotalCharge() << endl;

        eb=itsHamiltonian->GetTotalEnergy(itsCD.get());
        E=eb.GetTotalEnergy();
        dE=(E-Eold)/fabs(E);
        itsAccelerator->SetEnergy(E); //the ladder gates its hand-off on the energy change
        FD=itsAccelerator->GetError(); //i.e. [F,D]
        dFD=(FD-FDold);
        if (ipar.Verbose) DisplayEnergies(itsIterationCount,eb,relax,dFD,ChargeDensityChange,idealVirial);
        if (FD>FDold && fabs(dFD)>1e-9)
        {
            itsCD=cd_t(itsWaveFunction->GetChargeDensity()); //Get new charge density.
            ChargeDensityChange = itsCD->GetChangeFrom(*itsOldCD); //Get MaxAbs of change.
            itsCD->MixIn(*itsOldCD,1.0-relax/4.0);
            eb=itsHamiltonian->GetTotalEnergy(itsCD.get());
            relax*=0.8;
        }
        // if (E<Eold && Eold>Eoldold) relax*=0.5;
        if (FD<FDold ) relax*=1.5;
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
    cout << endl;
}

//-------------------------------------------------------------------------------------------------------------------------

template class tSCFIterator<double>;
template class tSCFIterator<dcmplx>;

} //namespace