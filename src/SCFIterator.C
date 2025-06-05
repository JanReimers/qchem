// File: SCFIteratorImplementation.C  Partial common implementation for an object that manages SCF convergence.



#include <SCFIterator.H>
#include "Imp/SCFAccelerator.H"
#include <WaveFunction.H>
#include <IterationParams.H>
#include <Hamiltonian.H>
#include <TotalEnergy.H>
#include <ChargeDensity.H>
#include <iostream>
#include <iomanip>
#include <cassert>

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::ios;


SCFIterator::SCFIterator(const BasisSet* bs, const ElectronConfiguration* ec,Hamiltonian* H,SCFAccelerator* acc,DM_CD* cd)
    : itsHamiltonian (H )
    , itsAccelerator (acc)       
    , itsWaveFunction(itsHamiltonian->CreateWaveFunction(bs,ec,*itsAccelerator) )
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
    itsWaveFunction->FillOrbitals();

    itsCD=itsWaveFunction->GetChargeDensity(); //Get new charge density.
    assert(itsCD);
    itsOldCD=cd;
    itsIterationCount=0;
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

bool SCFIterator::Iterate(const SCFIterationParams& ipar)
{
    assert(itsWaveFunction);
    assert(itsHamiltonian);
    assert(itsCD);
    if (ipar.Verbose)
    {
        cout << endl << endl;
        cout << " #             Etotal       Del(E)  Del(Ro)   2+V/K  relax" << endl;
        cout << "--------------------------------------------------------------------------------------------------------------------------" << endl;
    }

    double ChargeDensityChange=1;
    double Eold=0;
    // double Eoldold=0;
    double relax=ipar.StartingRelaxRo;
    double relMax=1.0;
    EnergyBreakdown eb;

    for (itsIterationCount=1; itsIterationCount<ipar.NMaxIter && ChargeDensityChange > ipar.MinDeltaRo; itsIterationCount++)
    {
        itsWaveFunction->DoSCFIteration(*itsHamiltonian,itsCD); //Just gets a set of eigen orbitals from the Hamiltonian
        itsWaveFunction->FillOrbitals();

        delete itsOldCD;
        itsOldCD=itsCD;
        itsCD=itsWaveFunction->GetChargeDensity(); //Get new charge density.
        ChargeDensityChange = itsCD->GetChangeFrom(*itsOldCD); //Get MaxAbs of change.
        itsCD->MixIn(*itsOldCD,1.0-relax);                           //relaxation.
        // cout << "Total charge=" << itsCD->GetTotalCharge() << endl;

        eb=itsHamiltonian->GetTotalEnergy(itsCD);
        double E=eb.GetTotalEnergy();
        if (ipar.Verbose) DisplayEnergies(itsIterationCount,eb,relax,E-Eold,ChargeDensityChange);
        if (E>Eold ) 
        {
            relax*=0.5;
            itsCD=itsWaveFunction->GetChargeDensity(); //Get new charge density.
            ChargeDensityChange = itsCD->GetChangeFrom(*itsOldCD); //Get MaxAbs of change.
            itsCD->MixIn(*itsOldCD,1.0-relax); 
            eb=itsHamiltonian->GetTotalEnergy(itsCD);
        }
        // if (E<Eold && Eold>Eoldold) relax*=0.5;
        if (E<Eold ) relax*=1.5;
        // if (E>Eold && Eold>Eoldold) relax*=1.5;
        if (relax>relMax) relax=relMax;

        // Eoldold=Eold;
        Eold=E;
    }
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
        << "Eee    : " << std::fixed << setw(w) << setprecision(nprec) << eb.Eee << "   "
        << "Eex     : " << std::fixed << setw(w) << setprecision(nprec) << eb.Exc << endl;
        cout << "Virial               V/K  : " << std::fixed << setw(w) << setprecision(nprec) << eb.GetVirial() << endl;
        

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
        
    cout << setw(8) << std::scientific << setw(8) << setprecision(1) << dE  << " ";
    cout << setw(8) << std::scientific << setw(8) << setprecision(1) << eb.GetVirial()+2.0 << " ";
    cout << setw(8) << std::scientific << setw(8) << setprecision(1) << dCD << " ";
    cout << setw(4) << std::fixed << setw(4) << setprecision(2) << relax << " ";
//    cout << setw(8) << setprecision(2) << fitError << " ";
//    cout << setw(9) << setprecision(2) << lam << " ";
    cout << endl;
}

//-------------------------------------------------------------------------------------------------------------------------

