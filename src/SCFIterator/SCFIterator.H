// File: SCFIterator.H  Interface for an object that manages SCF convergence.


#include "types.H"

class SCFIterator
{
public:
    virtual ~SCFIterator() {};

    virtual bool Iterate        (const SCFIterationParams& ipar)=0;
    virtual void DumpElectrons  (WaveFunction*, double Smear)=0;
    virtual void DisplayEnergies(int i, double lam, double chargeDensityChange, double fitError) const=0;
    virtual void DisplayEigen   () const=0;

    virtual ChargeDensity* GetExactChargeDensity () const=0;
};

