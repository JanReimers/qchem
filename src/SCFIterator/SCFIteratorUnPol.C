// File: UnPolarizedSCFIterator.C  SCF convergence for a UnPolarized wave function.



#include "Imp/SCFIterator/SCFIteratorUnPol.H"
#include "Imp/WaveFunction/ElectronDumper.H"
#include <IterationParams.H>
#include <Hamiltonian.H>
#include <WaveFunction.H>
#include <ChargeDensity.H>
#include <iostream>
#include <iomanip>
#include <cassert>


SCFIteratorUnPol::SCFIteratorUnPol(WaveFunction* W, Hamiltonian* H,ChargeDensity* guess,
                                               double nElectrons, double kT, bool showplot)
    : SCFIteratorImp(W,H,showplot)
    , itsTotalCharge(nElectrons)
    , itsEf(0)
{
    Initialize(guess,kT);
    assert(itsTotalCharge>0);
}

bool SCFIteratorUnPol::Iterate(const SCFIterationParams& ipar)
{
    if (ipar.Verbose)
    {
        std::cout << std::endl << std::endl;
        std::cout << " #        Etotal     Virial  K    V    Ven   Vee    Vxc    Del(Ro) Del(Vee)  Lambda     Ef(up) " << std::endl;
        std::cout << "-------------------------------------------------------------------------------" << std::endl;
    }
    bool ret=SCFIteratorImp::Iterate(ipar);

    if (ipar.Verbose)
    {
        std::cout << "-------------------------------------------------------------------------------" << std::endl;
        DisplayEigen();
    }
    return ret;
}

void SCFIteratorUnPol::DisplayEnergies(int i, double lam, double ChargeDensityChange, double fitError) const
{
    SCFIteratorImp::DisplayEnergies(i,lam,ChargeDensityChange,fitError);
//    std::cout.setf(std::ios::fixed,std::ios::floatfield);
//    std::cout << std::setw(9) << std::setprecision(6) << itsEf << std::endl;
    std::cout << std::endl;
}

void SCFIteratorUnPol::DumpElectrons(WaveFunction* wf, double kT)
{
    assert(wf);
    ElectronDumper ed  (0.000001,kT);
    wf->UpdateElectronDumper(ed);
    ed.DumpInElectrons(itsTotalCharge);  //Define occupations for all orbitals.
    itsEf=ed.GetFermiEnergy();
}

void SCFIteratorUnPol::DisplayEigen() const
{
    ElectronDumper ed  (0.000001,0.0);
    itsWaveFunction->UpdateElectronDumper(ed);
    ed.DumpInElectrons(itsTotalCharge);  //Define occupations for all orbitals.
    std::cout << ed << std::endl;
}


