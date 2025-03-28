// File: SCFIteratorUnPol.H  SCF convergence for a UnPolarized wave function.
#ifndef _SCFIteratorUnPol_H_
#define _SCFIteratorUnPol_H_


#include "Imp/SCFIterator/SCFIterator.H"

class SCFIteratorUnPol
    : virtual public SCFIteratorImp
{
public:
    SCFIteratorUnPol(WaveFunction*, Hamiltonian* ,Exact_CD* guess, double nElectrons);

    virtual bool Iterate(const SCFIterationParams& ipar);
    virtual void DisplayEnergies(int i, double lam, double chargeDensityChange, double fitError) const;
    virtual void DisplayEigen   () const;

private:
    double itsTotalCharge;
};

#endif // _SCFIteratorUnPol_H_
