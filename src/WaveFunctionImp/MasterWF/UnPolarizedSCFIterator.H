// File: UnPolarizedSCFIterator.H  SCF convergence for a UnPolarized wave function.



#include "WaveFunction/SCFIteratorImplementation.H"

class UnPolarizedSCFIterator
    : virtual public SCFIteratorImplementation
{
public:
    UnPolarizedSCFIterator(WaveFunction*, Hamiltonian* ,ChargeDensity* guess, double nElectrons, double kT, bool showplot);

    virtual bool   Iterate(double relax, double epsRo, int Nmax, double kT);
    virtual void   DumpElectrons  (WaveFunction*, double kT);
    virtual double DisplayEnergies(int i, double lam, double chargeDensityChange, double fitError) const;
    virtual void   DisplayEigen   () const;

private:
    double itsTotalCharge;
    double itsEf;
};
