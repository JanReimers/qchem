// File: PolarizedSCFIterator.H  SCF convergence for a UnPolarized wave function.
#ifndef _PolarizedSCFIterator_H_
#define _PolarizedSCFIterator_H_

#include "Imp/SCFIterator/SCFIterator.H"

class PolarizedWF;

class SCFIteratorPol
    : virtual public SCFIteratorImp
{
public:
    SCFIteratorPol(PolarizedWF* W, Hamiltonian*, ChargeDensity* guess, double nElectrons, double spin);

    virtual bool Iterate(const SCFIterationParams& ipar);
    virtual void DisplayEnergies(int i, double lam, double chargeDensityChange, double fitError) const;
    virtual void DisplayEigen   () const;
};

#endif // _PolarizedSCFIterator_H_
