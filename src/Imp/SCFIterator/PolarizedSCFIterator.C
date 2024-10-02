// File: PolarizedSCFIterator.C  SCF convergence for a Polarized wave fction.



#include "Hamiltonian.H"
#include "WaveFunction.H"
#include "Imp/SCFIterator/PolarizedSCFIterator.H"
#include "ChargeDensity.H"
#include "Orbital/ElectronDumper.H"
#include "Misc/Spin.H"
#include <iostream>
#include <iomanip>
#include <cassert>

PolarizedSCFIterator::PolarizedSCFIterator(PolarizedWF* W, Hamiltonian* H, ChargeDensity* guess,
                                           double nElectrons, double spin, double kT, bool showplot)
    : SCFIteratorImplementation(W,H,showplot)
    , itsTotalUp               (0.0)
    , itsTotalDown             (0.0)
    , itsUpEf                  (0)
    , itsDownEf                (0)
{
    DecideElectronCounts(nElectrons,spin);
    Initialize(guess,kT);
}

void PolarizedSCFIterator::DecideElectronCounts(double total, double spin)
{
    double a2=total+spin;
    double b2=total-spin;
    itsTotalUp=a2/2.0;
    itsTotalDown=b2/2.0;
    assert(fabs(itsTotalUp-floor(itsTotalUp))==0.0);
    assert(fabs(itsTotalDown-floor(itsTotalDown))==0.0);
}

void PolarizedSCFIterator::DumpElectrons(WaveFunction* wf, double kT)
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

bool PolarizedSCFIterator::Iterate(const SCFIterationParams& ipar)
{
    std::cout << std::endl << std::endl;
    std::cout << " #        Etotal     Virial  K    Vee    Vxc    Del(Ro) Del(Vee)  Lambda     Ef(up)   Ef(down) " << std::endl;
    std::cout << "----------------------------------------------------------------------------" << std::endl;
    bool ret=SCFIteratorImplementation::Iterate(ipar);
    std::cout << "-----------------------------------------------------------------" << std::endl;
    DisplayEigen();
    return ret;
}




double PolarizedSCFIterator::DisplayEnergies(int i, double lam, double ChargeDensityChange, double fitError) const
{
    double ret=SCFIteratorImplementation::DisplayEnergies(i,lam,ChargeDensityChange,fitError);
//    std::cout.setf(std::ios::fixed,std::ios::floatfield);
//    std::cout
//        << std::setw(9) << std::setprecision(6) << itsUpEf << " "
//        << std::setw(9) << std::setprecision(6) << itsDownEf << std::endl;
    std::cout << std::endl;
    return ret;
}

void PolarizedSCFIterator::DisplayEigen() const
{
    assert(itsWaveFunction);
    PolarizedWF* pwf=dynamic_cast<PolarizedWF*>(itsWaveFunction);
    assert(pwf);
    ElectronDumper uped  (0.0001,0.0);
    ElectronDumper downed(0.0001,0.0);
    pwf->GetWaveFunction(Spin::Up  )->UpdateElectronDumper(uped  );
    pwf->GetWaveFunction(Spin::Down)->UpdateElectronDumper(downed);
    uped  .DumpInElectrons(itsTotalUp  );  //Define occupations for all spin up   orbitals.
    downed.DumpInElectrons(itsTotalDown);  //Define occupations for all spin down orbitals.
    std::cout << "Alpha spin :" << std::endl;
    std::cout <<uped << std::endl;
    std::cout << "Beta spin :" << std::endl;
    std::cout << downed << std::endl;
}

