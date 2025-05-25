// File: SCFIteratorImplementation.C  Partial common implementation for an object that manages SCF convergence.



#include <SCFIterator.H>
#include <WaveFunction.H>
#include <IterationParams.H>
#include <Hamiltonian.H>
#include <TotalEnergy.H>
#include <ChargeDensity.H>
#include <iostream>
#include <iomanip>
#include <cassert>

SCFIterator::SCFIterator(WaveFunction* W,Hamiltonian* H,DM_CD* cd)
    : itsWaveFunction         (W )
    , itsHamiltonian          (H )
    , itsCD   (0)
    , itsOldCD(0)
{
    assert(itsWaveFunction);
    assert(itsHamiltonian);
    Initialize(cd);
}


void SCFIterator::Initialize(DM_CD* cd)
{
    itsWaveFunction->DoSCFIteration(*itsHamiltonian,cd);
    itsWaveFunction->FillOrbitals();

    itsCD=itsWaveFunction->GetChargeDensity(); //Get new charge density.
    assert(itsCD);
    itsOldCD=cd;
}
//
//  Recall that the wavefunction is not owned buy this.
//

SCFIterator::~SCFIterator()
{
    delete itsHamiltonian;
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
        std::cout << " #       Etotal        Virial          K             V            Ven          Vee          Vxc      Del(Ro) relax" << std::endl;
        std::cout << "---------------------------------------------------------------------------------------------------------------------" << std::endl;
    }

    double ChargeDensityChange=1;
    double Eold=0;
    double Eoldold=0;
    double relax=ipar.StartingRelaxRo;
    double relMax=1.0;

    for (size_t i=0; i<ipar.NMaxIter && ChargeDensityChange > ipar.MinDeltaRo; i++)
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
        if (ipar.Verbose) DisplayEnergies(i,eb,relax,ChargeDensityChange,0.0);
        double E=eb.GetTotalEnergy();
        if (E>Eold && Eold<Eoldold) relax*=0.5;
        if (E<Eold && Eold>Eoldold) relax*=0.5;
        if (E<Eold && Eold<Eoldold) relax*=1.5;
        if (E>Eold && Eold>Eoldold) relax*=1.5;
        if (relax>relMax) relax=relMax;

        Eoldold=Eold;
        Eold=E;
    }

    if (ipar.Verbose)
    {
        std::cout << "---------------------------------------------------------------------------------------------------------------------" << std::endl;
        DisplayEigen();
    }

    return ChargeDensityChange <= ipar.MinDeltaRo;
}


void SCFIterator::DisplayEigen() const
{
    itsWaveFunction->DisplayEigen();
}

DM_CD* SCFIterator::GetExactChargeDensity() const
{
    return itsCD;
}

using std::cout;
using std::setw;
using std::setprecision;
using std::ios;

void SCFIterator::DisplayEnergies(int i, const EnergyBreakdown& eb, double relax, double ChargeDensityChange, double fitError) const
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
    cout << setw(7) << setprecision(1) << ChargeDensityChange << " " << relax << " ";
//    cout << setw(8) << setprecision(2) << fitError << " ";
//    cout << setw(9) << setprecision(2) << lam << " ";
    cout << std::endl;
}

