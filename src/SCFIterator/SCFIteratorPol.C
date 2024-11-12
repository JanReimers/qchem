// File: PolarizedSCFIterator.C  SCF convergence for a Polarized wave fction.

#include "Imp/SCFIterator/SCFIteratorPol.H"
#include <Hamiltonian.H>
#include <WaveFunction.H>
#include <IterationParams.H>
#include <ChargeDensity.H>
#include <Spin.H>
#include <iostream>
#include <iomanip>
#include <cassert>

SCFIteratorPol::SCFIteratorPol(PolarizedWF* W, Hamiltonian* H, ChargeDensity* guess,
                                           double nElectrons, double spin)
    : SCFIteratorImp(W,H)
{
    Initialize(guess);
}

bool SCFIteratorPol::Iterate(const SCFIterationParams& ipar)
{
    if (ipar.Verbose)
    {
        std::cout << std::endl << std::endl;
        std::cout << " #        Etotal     Virial  K    V    Ven    Vee    Vxc    Del(Ro) Del(Vee)  Lambda     Ef(up)   Ef(down) " << std::endl;
        std::cout << "----------------------------------------------------------------------------" << std::endl;
    }
    bool ret=SCFIteratorImp::Iterate(ipar);
    if (ipar.Verbose)
    {
        std::cout << "-----------------------------------------------------------------" << std::endl;
        DisplayEigen();
    }
    return ret;
}




void SCFIteratorPol::DisplayEnergies(int i, double lam, double ChargeDensityChange, double fitError) const
{
    SCFIteratorImp::DisplayEnergies(i,lam,ChargeDensityChange,fitError);
//    std::cout.setf(std::ios::fixed,std::ios::floatfield);
//    std::cout
//        << std::setw(9) << std::setprecision(6) << itsUpEf << " "
//        << std::setw(9) << std::setprecision(6) << itsDownEf << std::endl;
    std::cout << std::endl;
}

void SCFIteratorPol::DisplayEigen() const
{
    assert(itsWaveFunction);
    itsWaveFunction->DisplayEigen();
}

