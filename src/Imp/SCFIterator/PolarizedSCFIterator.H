// File: PolarizedSCFIterator.H  SCF convergence for a UnPolarized wave function.
#ifndef _PolarizedSCFIterator_H_
#define _PolarizedSCFIterator_H_

#include "Imp/SCFIterator/SCFIteratorImplementation.H"


class PolarizedSCFIterator
    : virtual public SCFIteratorImplementation
{
public:
    PolarizedSCFIterator(PolarizedWF* W, Hamiltonian*, ChargeDensity* guess, double nElectrons, double spin, double kT, bool showplot);

    virtual bool Iterate(const SCFIterationParams& ipar);
    virtual void DumpElectrons  (WaveFunction*, double kT);
    virtual void DisplayEnergies(int i, double lam, double chargeDensityChange, double fitError) const;
    virtual void DisplayEigen   () const;

private:
    void DecideElectronCounts(double total, double spin);

    double itsTotalUp;
    double itsTotalDown;
    double itsUpEf;
    double itsDownEf;
};

#endif // _PolarizedSCFIterator_H_
