// File: SCFIteratorImplementation.C  Partial common implementation for an object that manages SCF convergence.



#include "WaveFunction.H"
#include "Imp/SCFIterator/SCFIteratorImplementation.H"
#include "Hamiltonian.H"
#include "TotalEnergy.H"
#include "FittedCD.H"
#include "Orbital/ElectronDumper.H"
#include "FunctionsImp/PlotterImplementation.H"
#include "Mesh/LinearMesh.H"
#include <iostream>
#include <iomanip>
#include <cassert>

LinearMesh mesh(-6,6,RVec(1,0,0),241);

SCFIteratorImplementation::SCFIteratorImplementation(WaveFunction* W,Hamiltonian* H,bool showplot)
    : itsWaveFunction         (W )
    , itsHamiltonian          (H )
    , itsExactChargeDensity   (0)
    , itsOldExactChargeDensity(0)
    , itsPlotter              (0)
{
    assert(itsWaveFunction);
    assert(itsHamiltonian);
    if (showplot) itsPlotter=new PlotterImplementation;
}


void SCFIteratorImplementation::Initialize(ChargeDensity* cd, double kT)
{
    assert(cd);

    itsHamiltonian->UseChargeDensity(cd);
    itsWaveFunction->DoSCFIteration(*itsHamiltonian);
    DumpElectrons(itsWaveFunction,kT);

    itsExactChargeDensity=itsWaveFunction->GetChargeDensity(); //Get new charge density.
    assert(itsExactChargeDensity);
    if(itsPlotter)
    {
        itsPlotter->Lines();
        itsPlotter->Plot(*itsExactChargeDensity,&mesh,RVec(1,0,0));
    }

    itsHamiltonian->UseChargeDensity(itsExactChargeDensity);
//	itsHamiltonian->GetTotalEnergy();
}
//
//  Recall that the wavefunction is not owned buy this.
//

SCFIteratorImplementation::~SCFIteratorImplementation()
{
    delete itsHamiltonian;
    delete itsExactChargeDensity;
    delete itsOldExactChargeDensity;
    delete itsPlotter;
}

bool SCFIteratorImplementation::Iterate(double relax, double epsRo, int Nmax, double Smear)
{
    assert(itsWaveFunction);
    assert(itsHamiltonian);
    assert(itsExactChargeDensity);

    double ChargeDensityChange=1;
    double Eold=0;
    double Eoldold=0;
    double relMax=relax;

    for (int i=0; i<Nmax && ChargeDensityChange > epsRo; i++)
    {
        itsWaveFunction->DoSCFIteration(*itsHamiltonian); //Just gets a set of eigen orbitals from the Hamiltonian
        DumpElectrons(itsWaveFunction,Smear);

        delete itsOldExactChargeDensity;
        itsOldExactChargeDensity=itsExactChargeDensity;
        itsExactChargeDensity=itsWaveFunction->GetChargeDensity(); //Get new charge density.
        ChargeDensityChange = itsExactChargeDensity->GetChangeFrom(*itsOldExactChargeDensity); //Get MaxAbs of change.
        itsExactChargeDensity->MixIn(*itsOldExactChargeDensity,1.0-relax);                           //relaxation.

        if(itsPlotter)
        {
            itsPlotter->Lines();
            itsPlotter->Plot(*itsExactChargeDensity,&mesh,RVec(1,0,0));
            int n;
            std::cin >> n;
        }

        itsHamiltonian->UseChargeDensity(itsExactChargeDensity);      //Set all the potentials for this charge denisty distribution.

        double E=DisplayEnergies(i,0.0,ChargeDensityChange,0.0);
        if (E>Eold && Eold<Eoldold) relax*=0.5;
        if (E<Eold && Eold>Eoldold) relax*=0.5;
        if (E<Eold && Eold<Eoldold) relax*=1.2;
        if (E>Eold && Eold>Eoldold) relax*=1.2;
        if (relax>relMax) relax=relMax;

        Eoldold=Eold;
        Eold=E;
    }

    return ChargeDensityChange <= epsRo;
}


ChargeDensity* SCFIteratorImplementation::GetExactChargeDensity() const
{
    return itsExactChargeDensity;
}

using std::cout;
using std::setw;
using std::setprecision;
using std::ios;

double SCFIteratorImplementation::DisplayEnergies(int i, double lam, double ChargeDensityChange, double fitError) const
{
    TotalEnergy te = itsHamiltonian->GetTotalEnergy();

    cout.setf(ios::fixed,ios::floatfield);
    cout << setw(2)  << i << " "
         << setw(16) << setprecision(8) << te.GetTotalEnergy()
         << setw(6)  << setprecision(3) << te.GetVirial() << " "
         << setw(6)  << setprecision(3) << te.Kinetic << " "
         << setw(6)  << setprecision(3) << te.Eee << " "
         << setw(6)  << setprecision(3) << te.Exc << " ";
    cout.setf(ios::scientific,ios::floatfield);
    cout << setw(8) << setprecision(2) << ChargeDensityChange << " ";
    cout << setw(8) << setprecision(2) << fitError << " ";
    cout << setw(9) << setprecision(2) << lam << " ";
    return te.GetTotalEnergy();
}

