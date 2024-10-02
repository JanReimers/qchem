// File: UnPolarizedSCFIterator.C  SCF convergence for a UnPolarized wave function.



#include "Imp/SCFIterator/UnPolarizedSCFIterator.H"
#include "Hamiltonian.H"
#include "WaveFunction.H"
#include "Orbital/ElectronDumper.H"
#include "ChargeDensity.H"
#include <iostream>
#include <iomanip>
#include <cassert>


UnPolarizedSCFIterator::UnPolarizedSCFIterator(WaveFunction* W, Hamiltonian* H,ChargeDensity* guess,
                                               double nElectrons, double kT, bool showplot)
    : SCFIteratorImplementation(W,H,showplot)
    , itsTotalCharge(nElectrons)
    , itsEf(0)
{
    Initialize(guess,kT);
    assert(itsTotalCharge>0);
}

bool UnPolarizedSCFIterator::Iterate(const SCFIterationParams& ipar)
{
    std::cout << std::endl << std::endl;
    std::cout << " #        Etotal     Virial  K    Vee    Vxc    Del(Ro) Del(Vee)  Lambda     Ef(up) " << std::endl;
    std::cout << "-------------------------------------------------------------------------------" << std::endl;

    bool ret=SCFIteratorImplementation::Iterate(ipar);

    std::cout << "-------------------------------------------------------------------------------" << std::endl;
    DisplayEigen();
    return ret;
}

double UnPolarizedSCFIterator::DisplayEnergies(int i, double lam, double ChargeDensityChange, double fitError) const
{
    double ret=SCFIteratorImplementation::DisplayEnergies(i,lam,ChargeDensityChange,fitError);
//    std::cout.setf(std::ios::fixed,std::ios::floatfield);
//    std::cout << std::setw(9) << std::setprecision(6) << itsEf << std::endl;
    std::cout << std::endl;
    return ret;
}

void UnPolarizedSCFIterator::DumpElectrons(WaveFunction* wf, double kT)
{
    assert(wf);
    ElectronDumper ed  (0.000001,kT);
    wf->UpdateElectronDumper(ed);
    ed.DumpInElectrons(itsTotalCharge);  //Define occupations for all orbitals.
    itsEf=ed.GetFermiEnergy();
}

void UnPolarizedSCFIterator::DisplayEigen() const
{
    ElectronDumper ed  (0.000001,0.0);
    itsWaveFunction->UpdateElectronDumper(ed);
    ed.DumpInElectrons(itsTotalCharge);  //Define occupations for all orbitals.
    std::cout << ed << std::endl;
}


