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
    assert(cd);

    itsHamiltonian->UseChargeDensity(cd);
    itsWaveFunction->DoSCFIteration(*itsHamiltonian,cd);
    itsWaveFunction->FillOrbitals(0);

    itsCD=itsWaveFunction->GetChargeDensity(); //Get new charge density.
    assert(itsCD);
    itsHamiltonian->UseChargeDensity(itsCD);
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
        std::cout << " #       Etotal        Virial       K        V       Ven       Vee       Vxc    Del(Ro)" << std::endl;
        std::cout << "----------------------------------------------------------------------------------------" << std::endl;
    }

    double ChargeDensityChange=1;
    double Eold=0;
    double Eoldold=0;
    double relax;
    double relMax=relax=ipar.StartingRelaxRo;

    for (size_t i=0; i<ipar.NMaxIter && ChargeDensityChange > ipar.MinDeltaRo; i++)
    {
        itsWaveFunction->DoSCFIteration(*itsHamiltonian,itsCD); //Just gets a set of eigen orbitals from the Hamiltonian
        itsWaveFunction->FillOrbitals(0);

        delete itsOldCD;
        itsOldCD=itsCD;
        itsCD=itsWaveFunction->GetChargeDensity(); //Get new charge density.
        ChargeDensityChange = itsCD->GetChangeFrom(*itsOldCD); //Get MaxAbs of change.
        itsCD->MixIn(*itsOldCD,1.0-relax);                           //relaxation.
        // std::cout << "Total charge=" << itsCD->GetTotalCharge() << std::endl;
        itsHamiltonian->UseChargeDensity(itsCD);      //Set all the potentials for this charge denisty distribution.

        if (ipar.Verbose) DisplayEnergies(i,0.0,ChargeDensityChange,0.0,itsCD);
        double E=itsHamiltonian->GetTotalEnergy(itsCD).GetTotalEnergy();
        if (E>Eold && Eold<Eoldold) relax*=0.5;
        if (E<Eold && Eold>Eoldold) relax*=0.5;
        if (E<Eold && Eold<Eoldold) relax*=1.2;
        if (E>Eold && Eold>Eoldold) relax*=1.2;
        if (relax>relMax) relax=relMax;

        Eoldold=Eold;
        Eold=E;
    }

    if (ipar.Verbose)
    {
        std::cout << "------------------------------------------------------------------------------------------" << std::endl;
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

void SCFIterator::DisplayEnergies(int i, double lam, double ChargeDensityChange, double fitError,const DM_CD* cd) const
{
    TotalEnergy te = itsHamiltonian->GetTotalEnergy(cd);

    cout.setf(ios::fixed,ios::floatfield);
    cout << setw(2)  << i << " "
         << setw(16) << setprecision(8) << te.GetTotalEnergy() << " "
         << setw(10) << setprecision(8) << te.GetVirial() << " "
         << setw(8)  << setprecision(8) << te.Kinetic << " "
         << setw(8)  << setprecision(8) << te.GetPotentialEnergy() << " "
         << setw(8)  << setprecision(8) << te.Een << " "
         << setw(8)  << setprecision(8) << te.Eee << " "
         << setw(8)  << setprecision(8) << te.Exc << " ";
    cout.setf(ios::scientific,ios::floatfield);
    cout << setw(8) << setprecision(2) << ChargeDensityChange << " ";
//    cout << setw(8) << setprecision(2) << fitError << " ";
//    cout << setw(9) << setprecision(2) << lam << " ";
    cout << std::endl;
}

