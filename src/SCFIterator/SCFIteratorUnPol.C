// File: UnPolarizedSCFIterator.C  SCF convergence for a UnPolarized wave function.



#include "Imp/SCFIterator/SCFIteratorUnPol.H"
#include <IterationParams.H>
#include <Hamiltonian.H>
#include <WaveFunction.H>
#include <ChargeDensity.H>
#include <iostream>
#include <iomanip>
#include <cassert>


SCFIteratorUnPol::SCFIteratorUnPol(WaveFunction* W, Hamiltonian* H,ChargeDensity* guess,
                                               double nElectrons)
    : SCFIteratorImp(W,H)
    , itsTotalCharge(nElectrons)
{
    Initialize(guess);
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


void SCFIteratorUnPol::DisplayEigen() const
{
    itsWaveFunction->DisplayEigen();
}


