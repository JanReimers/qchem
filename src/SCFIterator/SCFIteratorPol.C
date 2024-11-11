// File: PolarizedSCFIterator.C  SCF convergence for a Polarized wave fction.

#include "Imp/SCFIterator/SCFIteratorPol.H"
#include "Imp/WaveFunction/ElectronDumper.H"
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
    , itsTotalUp               (0.0)
    , itsTotalDown             (0.0)
    , itsUpEf                  (0)
    , itsDownEf                (0)
{
    DecideElectronCounts(nElectrons,spin);
    Initialize(guess);
}

void SCFIteratorPol::DecideElectronCounts(double total, double spin)
{
    double a2=total+spin;
    double b2=total-spin;
    itsTotalUp=a2/2.0;
    itsTotalDown=b2/2.0;
//    assert(fabs(itsTotalUp-floor(itsTotalUp))==0.0);
//    assert(fabs(itsTotalDown-floor(itsTotalDown))==0.0);
}

void SCFIteratorPol::DumpElectrons(WaveFunction* wf, double kT)
{
    assert(wf);
    PolarizedWF* pwf=dynamic_cast<PolarizedWF*>(wf);
    assert(pwf);

    ElectronDumper uped  (0.0001,kT);
    ElectronDumper downed(0.0001,kT);
    pwf->GetWaveFunction(Spin::Up  )->UpdateElectronDumper(uped  );
    pwf->GetWaveFunction(Spin::Down)->UpdateElectronDumper(downed);
    uped  .DumpInElectrons(itsTotalUp  );  //Define occupations for all spin up   orbitals.
    downed.DumpInElectrons(itsTotalDown);  //Define occupations for all spin down orbitals.
    itsUpEf  =uped  .GetFermiEnergy();
    itsDownEf=downed.GetFermiEnergy();
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
//    PolarizedWF* pwf=dynamic_cast<PolarizedWF*>(itsWaveFunction);
//    assert(pwf);
//    ElectronDumper uped  (0.0001,0.0);
//    ElectronDumper downed(0.0001,0.0);
//    pwf->GetWaveFunction(Spin::Up  )->UpdateElectronDumper(uped  );
//    pwf->GetWaveFunction(Spin::Down)->UpdateElectronDumper(downed);
//    uped  .DumpInElectrons(itsTotalUp  );  //Define occupations for all spin up   orbitals.
//    downed.DumpInElectrons(itsTotalDown);  //Define occupations for all spin down orbitals.
//    std::cout << "Alpha spin :" << std::endl;
//    std::cout <<uped << std::endl;
//    std::cout << "Beta spin :" << std::endl;
//    std::cout << downed << std::endl;
}

