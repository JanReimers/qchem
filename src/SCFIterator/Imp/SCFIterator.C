// File: SCFIterator/Imp/SCFIterator.C  Partial common implementation for an object that manages SCF convergence.
module;
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <string>
#include <cctype>

module qchem.SCFIterator;
import qchem.SCFParams;
import qchem.SCFAccelerator;

import qchem.WaveFunction;
import qchem.WaveFunction.Factory;

import qchem.Hamiltonian;
import qchem.Energy; 
import qchem.ChargeDensity;

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
static std::string ConfigString(const qchem::WaveFunction::WaveFunction* wf)
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


SCFIterator::SCFIterator(const bs_t* bs, const ElectronConfiguration* ec,::Hamiltonian* H,SCFAccelerator* acc,DM_CD* cd)
    : itsHamiltonian (H )
    , itsAccelerator (acc)       
    , itsWaveFunction(qchem::WaveFunction::Factory(itsHamiltonian,bs,ec,itsAccelerator) )
    , itsCD          (0)
    , itsOldCD       (0)
    , itsIterationCount(0)
    , itsConverged(false)
{
    assert(itsHamiltonian);
    assert(itsWaveFunction);
    Initialize(cd);
}


void SCFIterator::Initialize(DM_CD* cd)
{
    itsWaveFunction->DoSCFIteration(*itsHamiltonian,cd);
    itsWaveFunction->FillOrbitals(0.0001);

    itsCD=itsWaveFunction->GetChargeDensity(); //Get new charge density.
    assert(itsCD);
    itsOldCD=cd;
    itsIterationCount=0;
    itsConverged=false;
    // DisplayEnergies(itsIterationCount,itsHamiltonian->GetTotalEnergy(itsCD),1.0,0.0,0.0);
}
//
//  Recall that the wavefunction is not owned buy this.
//

SCFIterator::~SCFIterator()
{
    delete itsHamiltonian;
    delete itsAccelerator;
    delete itsWaveFunction;
    delete itsCD;
    delete itsOldCD;
}

bool SCFIterator::Iterate(const SCFParams& ipar)
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
            delete itsOldCD;
            itsOldCD=itsCD;
            itsCD=DirectMinStep(Eold,ipar.MergeTol);
            ChargeDensityChange = itsCD->GetChangeFrom(*itsOldCD)/itsCD->GetTotalCharge();
        }
        else
        {
            itsWaveFunction->DoSCFIteration(*itsHamiltonian,itsCD); //Just gets a set of eigen orbitals from the Hamiltonian
            itsWaveFunction->FillOrbitals(ipar.MergeTol);

            delete itsOldCD;
            itsOldCD=itsCD;
            itsCD=itsWaveFunction->GetChargeDensity(); //Get new charge density.
            ChargeDensityChange = itsCD->GetChangeFrom(*itsOldCD)/itsCD->GetTotalCharge(); //Get relative MaxAbs of change.
            if (ChargeDensityChange<1e-5) relMax=0.5;
            itsCD->MixIn(*itsOldCD,1.0-relax);                           //relaxation.
        }
        // cout << "Total charge=" << itsCD->GetTotalCharge() << endl;

        eb=itsHamiltonian->GetTotalEnergy(itsCD);
        E=eb.GetTotalEnergy();
        dE=(E-Eold)/fabs(E);
        itsAccelerator->SetEnergy(E); //the ladder gates its hand-off on the energy change
        FD=itsAccelerator->GetError(); //i.e. [F,D]
        dFD=(FD-FDold);
        if (ipar.Verbose) DisplayEnergies(itsIterationCount,eb,relax,dFD,ChargeDensityChange,idealVirial);
        if (FD>FDold && fabs(dFD)>1e-9) 
        {
            delete itsCD;
            itsCD=itsWaveFunction->GetChargeDensity(); //Get new charge density.
            ChargeDensityChange = itsCD->GetChangeFrom(*itsOldCD); //Get MaxAbs of change.
            itsCD->MixIn(*itsOldCD,1.0-relax/4.0); 
            eb=itsHamiltonian->GetTotalEnergy(itsCD);
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
DM_CD* SCFIterator::DirectMinStep(double Ecur, double mergeTol)
{
    if (!itsWaveFunction->BuildFockAndComputeSteps(*itsHamiltonian,itsCD))
    {
        itsWaveFunction->DoSCFIteration(*itsHamiltonian,itsCD);
        itsWaveFunction->FillOrbitals(mergeTol);
        return itsWaveFunction->GetChargeDensity();
    }
    double t=1.0;
    for (int k=0;k<12;k++)
    {
        itsWaveFunction->MoveOrbitals(t,false,mergeTol);                 //trial
        DM_CD* cdt=itsWaveFunction->GetChargeDensity();
        double Et=itsHamiltonian->GetTotalEnergy(cdt).GetTotalEnergy();
        delete cdt;
        if (Et<Ecur) break;
        t*=0.5;
    }
    itsWaveFunction->MoveOrbitals(t,true,mergeTol);                      //commit at t
    return itsWaveFunction->GetChargeDensity();
}


void SCFIterator::DisplayEigen() const
{
    itsWaveFunction->DisplayEigen();
}

EnergyBreakdown SCFIterator::GetEnergy() const
{
    return itsHamiltonian->GetTotalEnergy(itsCD);
}



void SCFIterator::DisplayEnergies(int i, const EnergyBreakdown& eb, double relax, double dE, double dCD, size_t idealVirial) const
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

} //namespace