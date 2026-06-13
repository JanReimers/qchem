// File: SCFIterator/Imp/SCFIterator.C  Partial common implementation for an object that manages SCF convergence.
module;
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>

module qchem.SCFIterator;
import qchem.SCFParams;
import qchem.SCFAccelerator;

import qchem.WaveFunction;
import qchem.WaveFunction.Factory;

import qchem.Hamiltonian;
import qchem.Energy; 
import qchem.ChargeDensity;

import qchem.Symmetry.ElectronConfiguration;
import qchem.Math;

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::ios;


namespace qchem::SCFIterator
{


SCFIterator::SCFIterator(const bs_t* bs, const ElectronConfiguration* ec,class Hamiltonian* H,SCFAccelerator* acc,DM_CD* cd)
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
        cout << "   relax" << endl;
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
        itsWaveFunction->DoSCFIteration(*itsHamiltonian,itsCD); //Just gets a set of eigen orbitals from the Hamiltonian
        itsWaveFunction->FillOrbitals(ipar.MergeTol);

        delete itsOldCD;
        itsOldCD=itsCD;
        itsCD=itsWaveFunction->GetChargeDensity(); //Get new charge density.
        ChargeDensityChange = itsCD->GetChangeFrom(*itsOldCD)/itsCD->GetTotalCharge(); //Get relative MaxAbs of change.
        if (ChargeDensityChange<1e-5) relMax=0.5;
        itsCD->MixIn(*itsOldCD,1.0-relax);                           //relaxation.
        // cout << "Total charge=" << itsCD->GetTotalCharge() << endl;

        eb=itsHamiltonian->GetTotalEnergy(itsCD);
        E=eb.GetTotalEnergy();
        dE=(E-Eold)/fabs(E);
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
    cout << setw(4) << std::fixed << setw(4) << setprecision(2) << relax << " ";
    cout << endl;
}

//-------------------------------------------------------------------------------------------------------------------------

} //namespace