// File: SCFIteratorImplementation.C  Partial common implementation for an object that manages SCF convergence.



#include "Imp/SCFIterator/SCFIterator.H"
#include <WaveFunction.H>
#include <IterationParams.H>
#include <Hamiltonian.H>
#include <TotalEnergy.H>
#include <ChargeDensity.H>
#include <iostream>
#include <iomanip>
#include <cassert>

SCFIteratorImp::SCFIteratorImp(WaveFunction* W,Hamiltonian* H)
    : itsWaveFunction         (W )
    , itsHamiltonian          (H )
    , itsExactChargeDensity   (0)
    , itsOldExactChargeDensity(0)
{
    assert(itsWaveFunction);
    assert(itsHamiltonian);
}


void SCFIteratorImp::Initialize(ChargeDensity* cd)
{
    assert(cd);

    itsHamiltonian->UseChargeDensity(cd);
    itsWaveFunction->DoSCFIteration(*itsHamiltonian);
    itsWaveFunction->FillOrbitals(0,Spin::None);

    itsExactChargeDensity=itsWaveFunction->GetChargeDensity(); //Get new charge density.
    assert(itsExactChargeDensity);
    itsHamiltonian->UseChargeDensity(itsExactChargeDensity);
    itsOldExactChargeDensity=cd;
}
//
//  Recall that the wavefunction is not owned buy this.
//

SCFIteratorImp::~SCFIteratorImp()
{
    delete itsHamiltonian;
    delete itsExactChargeDensity;
    delete itsOldExactChargeDensity;
}

bool SCFIteratorImp::Iterate(const SCFIterationParams& ipar)
{
    assert(itsWaveFunction);
    assert(itsHamiltonian);
    assert(itsExactChargeDensity);

    double ChargeDensityChange=1;
    double Eold=0;
    double Eoldold=0;
    double relax;
    double relMax=relax=ipar.StartingRelaxRo;

    for (size_t i=0; i<ipar.NMaxIter && ChargeDensityChange > ipar.MinDeltaRo; i++)
    {
        itsWaveFunction->DoSCFIteration(*itsHamiltonian); //Just gets a set of eigen orbitals from the Hamiltonian
        itsWaveFunction->FillOrbitals(0,Spin::None);

        delete itsOldExactChargeDensity;
        itsOldExactChargeDensity=itsExactChargeDensity;
        itsExactChargeDensity=itsWaveFunction->GetChargeDensity(); //Get new charge density.
        ChargeDensityChange = itsExactChargeDensity->GetChangeFrom(*itsOldExactChargeDensity); //Get MaxAbs of change.
        itsExactChargeDensity->MixIn(*itsOldExactChargeDensity,1.0-relax);                           //relaxation.

        itsHamiltonian->UseChargeDensity(itsExactChargeDensity);      //Set all the potentials for this charge denisty distribution.

        if (ipar.Verbose) DisplayEnergies(i,0.0,ChargeDensityChange,0.0);
        double E=itsHamiltonian->GetTotalEnergy().GetTotalEnergy();
        if (E>Eold && Eold<Eoldold) relax*=0.5;
        if (E<Eold && Eold>Eoldold) relax*=0.5;
        if (E<Eold && Eold<Eoldold) relax*=1.2;
        if (E>Eold && Eold>Eoldold) relax*=1.2;
        if (relax>relMax) relax=relMax;

        Eoldold=Eold;
        Eold=E;
    }

    return ChargeDensityChange <= ipar.MinDeltaRo;
}


ChargeDensity* SCFIteratorImp::GetExactChargeDensity() const
{
    return itsExactChargeDensity;
}

using std::cout;
using std::setw;
using std::setprecision;
using std::ios;

void SCFIteratorImp::DisplayEnergies(int i, double lam, double ChargeDensityChange, double fitError) const
{
    TotalEnergy te = itsHamiltonian->GetTotalEnergy();

    cout.setf(ios::fixed,ios::floatfield);
    cout << setw(2)  << i << " "
         << setw(16) << setprecision(8) << te.GetTotalEnergy()
         << setw(6)  << setprecision(3) << te.GetVirial() << " "
         << setw(6)  << setprecision(3) << te.Kinetic << " "
         << setw(6)  << setprecision(3) << te.GetPotentialEnergy() << " "
         << setw(6)  << setprecision(3) << te.Een << " "
         << setw(6)  << setprecision(3) << te.Eee << " "
         << setw(6)  << setprecision(3) << te.Exc << " ";
    cout.setf(ios::scientific,ios::floatfield);
    cout << setw(8) << setprecision(2) << ChargeDensityChange << " ";
    cout << setw(8) << setprecision(2) << fitError << " ";
    cout << setw(9) << setprecision(2) << lam << " ";
}

