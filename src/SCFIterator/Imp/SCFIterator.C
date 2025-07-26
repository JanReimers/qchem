// File: SCFIterator/Imp/SCFIterator.C  Partial common implementation for an object that manages SCF convergence.
module;
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
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

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::ios;




SCFIterator::SCFIterator(const BasisSet* bs, const ElectronConfiguration* ec,Hamiltonian* H,SCFAccelerator* acc,DM_CD* cd)
    : itsHamiltonian (H )
    , itsAccelerator (acc)       
    , itsWaveFunction(WaveFunctionF::Factory(itsHamiltonian,bs,ec,itsAccelerator) )
    , itsCD          (0)
    , itsOldCD       (0)
    , itsIterationCount(0)
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
    if (ipar.Verbose)
    {
        cout << endl << endl;
        cout << " #             Etotal       2+V/K    Del(E)  Del(Ro)";
        itsAccelerator->ShowLabels(cout);
        cout << "   relax" << endl;
        cout << "                                     ";
        cout << "(" << setw(8) << std::scientific << setw(5) << setprecision(0) << ipar.MinDelE  << ") ";
        cout << "(" << setw(8) << std::scientific << setw(5) << setprecision(0) << ipar.MinDeltaRo  << ") ";
        cout << "(" << setw(8) << std::scientific << setw(5) << setprecision(0) << ipar.MinError  << ") ";
        cout << endl;
        cout << "--------------------------------------------------------------------------------------------------------------------------" << endl;
    }

    double ChargeDensityChange=1;
    double E=0,Eold=0,dE=1e10,Error=1e10;
    // double Eoldold=0;
    double relax=ipar.StartingRelaxRo;
    double relMax=0.8;
    EnergyBreakdown eb;

    for (itsIterationCount=1; 
        itsIterationCount   < ipar.NMaxIter && ( 
        ChargeDensityChange > ipar.MinDeltaRo ||
        fabs(dE)            > ipar.MinDelE ||
        Error               > ipar.MinError );
        itsIterationCount++)
    {
        itsWaveFunction->DoSCFIteration(*itsHamiltonian,itsCD); //Just gets a set of eigen orbitals from the Hamiltonian
        itsWaveFunction->FillOrbitals(ipar.MergeTol);

        delete itsOldCD;
        itsOldCD=itsCD;
        itsCD=itsWaveFunction->GetChargeDensity(); //Get new charge density.
        ChargeDensityChange = itsCD->GetChangeFrom(*itsOldCD); //Get MaxAbs of change.
        if (ChargeDensityChange<1e-2) relMax=1.0;
        itsCD->MixIn(*itsOldCD,1.0-relax);                           //relaxation.
        // cout << "Total charge=" << itsCD->GetTotalCharge() << endl;

        eb=itsHamiltonian->GetTotalEnergy(itsCD);
        E=eb.GetTotalEnergy();
        dE=E-Eold;
        Error=itsAccelerator->GetError(); //i.e. [F,D]
        if (ipar.Verbose) DisplayEnergies(itsIterationCount,eb,relax,dE,ChargeDensityChange);
        if (E>Eold ) 
        {
            itsCD=itsWaveFunction->GetChargeDensity(); //Get new charge density.
            ChargeDensityChange = itsCD->GetChangeFrom(*itsOldCD); //Get MaxAbs of change.
            itsCD->MixIn(*itsOldCD,1.0-relax/4.0); 
            eb=itsHamiltonian->GetTotalEnergy(itsCD);
            relax*=0.8;
        }
        // if (E<Eold && Eold>Eoldold) relax*=0.5;
        if (E<Eold ) relax*=1.5;
        // if (E>Eold && Eold>Eoldold) relax*=1.5;
        if (relax>relMax) relax=relMax;

        // Eoldold=Eold;
        Eold=E;
    }
    itsIterationCount--;
    size_t nprec=14,ndigits=log10(-eb.Een)+1,w=1+ndigits+1+nprec;
    if (ipar.Verbose)
    {
        cout << "--------------------------------------------------------------------------------------------------------------------------" << endl;
        cout << "Energy    Breakdown  "
        << "Total: " << std::fixed << setw(w) << setprecision(nprec) << eb.GetTotalEnergy() << "  "
        << "Kinetic: " << std::fixed << setw(w) << setprecision(nprec) << eb.Kinetic << "  "
        << "Potential: " << std::fixed << setw(w) << setprecision(nprec) << eb.GetPotentialEnergy() << endl;
        cout << "Potential Breakdown  "
        << "Een  : " << std::fixed << setw(w) << setprecision(nprec) << eb.Een << "  "
        << "Eee    : " << std::fixed << setw(w) << setprecision(nprec) << eb.Eee << "  "
        << "Eex      : " << std::fixed << setw(w) << setprecision(nprec) << eb.Exc << endl;
        cout << "Virial               V/K  : " << std::fixed << setw(w) << setprecision(nprec) << eb.GetVirial() << "  ";
        if (eb.Exc!=0.0)
            cout << "Eee/Exc: " << std::fixed << setw(w) << setprecision(nprec) << eb.Eee/eb.Exc << endl;

        DisplayEigen();
    }

    return ChargeDensityChange <= ipar.MinDeltaRo;
}


void SCFIterator::DisplayEigen() const
{
    itsWaveFunction->DisplayEigen();
}

EnergyBreakdown SCFIterator::GetEnergy() const
{
    return itsHamiltonian->GetTotalEnergy(itsCD);
}



void SCFIterator::DisplayEnergies(int i, const EnergyBreakdown& eb, double relax, double dE, double dCD) const
{
    cout.setf(ios::fixed,ios::floatfield);
    cout << setw(3)  << i << " ";
    cout << setw(2+6+14) << setprecision(14) << eb.GetTotalEnergy() << " ";
    cout << setw(8) << std::scientific << setw(8) << setprecision(1) << eb.GetVirial()+2.0 << " ";
        
    cout << setw(8) << std::scientific << setw(8) << setprecision(1) << dE  << " ";
    cout << setw(8) << std::scientific << setw(7) << setprecision(1) << dCD << " ";
    itsAccelerator->ShowConvergence(cout);
    cout << setw(4) << std::fixed << setw(4) << setprecision(2) << relax << " ";
    cout << endl;
}

//-------------------------------------------------------------------------------------------------------------------------

