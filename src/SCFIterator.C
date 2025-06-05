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
        std::cout << std::endl << std::endl;
        std::cout << " #       Etotal        Virial          K             V            Ven          Vee          Vxc      Del(E)  Del(Ro) relax" << std::endl;
        std::cout << "--------------------------------------------------------------------------------------------------------------------------" << std::endl;
    }

    double ChargeDensityChange=1;
    double Eold=0;
    double Eoldold=0;
    double relax=ipar.StartingRelaxRo;
    double relMax=1.0;

    for (itsIterationCount=1; itsIterationCount<ipar.NMaxIter && ChargeDensityChange > ipar.MinDeltaRo; itsIterationCount++)
    {
        itsWaveFunction->DoSCFIteration(*itsHamiltonian,itsCD); //Just gets a set of eigen orbitals from the Hamiltonian
        itsWaveFunction->FillOrbitals();

        delete itsOldCD;
        itsOldCD=itsCD;
        itsCD=itsWaveFunction->GetChargeDensity(); //Get new charge density.
        ChargeDensityChange = itsCD->GetChangeFrom(*itsOldCD); //Get MaxAbs of change.
        itsCD->MixIn(*itsOldCD,1.0-relax);                           //relaxation.
        // std::cout << "Total charge=" << itsCD->GetTotalCharge() << std::endl;

        EnergyBreakdown eb=itsHamiltonian->GetTotalEnergy(itsCD);
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

        Eoldold=Eold;
        Eold=E;
    }

    if (ipar.Verbose)
    {
        std::cout << "--------------------------------------------------------------------------------------------------------------------------" << std::endl;
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


using std::cout;
using std::setw;
using std::setprecision;
using std::ios;

void SCFIterator::DisplayEnergies(int i, const EnergyBreakdown& eb, double relax, double dE, double dCD) const
{
    cout.setf(ios::fixed,ios::floatfield);
    cout << setw(2)  << i << " "
         << setw(16) << setprecision(10) << eb.GetTotalEnergy() << " "
         << setw(12) << setprecision(10) << eb.GetVirial() << " "
         << setw(12)  << setprecision(8) << eb.Kinetic << " "
         << setw(12)  << setprecision(8) << eb.GetPotentialEnergy() << " "
         << setw(12)  << setprecision(8) << eb.Een << " "
         << setw(12)  << setprecision(8) << eb.Eee << " "
         << setw(12)  << setprecision(8) << eb.Exc << " ";
    cout.setf(ios::scientific,ios::floatfield);
    cout << setw(7) << setprecision(1) << dE  << " ";
    cout << setw(7) << setprecision(1) << dCD << " " << relax << " ";
//    cout << setw(8) << setprecision(2) << fitError << " ";
//    cout << setw(9) << setprecision(2) << lam << " ";
    cout << std::endl;
}

//-------------------------------------------------------------------------------------------------------------------------

