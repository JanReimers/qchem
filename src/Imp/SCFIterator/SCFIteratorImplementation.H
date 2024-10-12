// File: SCFIteratorImplementation.H  Partial common implementation for an object that manages SCF convergence.

#include "SCFIterator/SCFIterator.H"

class SCFIteratorImplementation
    : public virtual SCFIterator
{
public:
    SCFIteratorImplementation(WaveFunction*, Hamiltonian*,bool showplot);
    virtual ~SCFIteratorImplementation();
    virtual bool   Iterate(const SCFIterationParams& ipar);
    virtual double DisplayEnergies(int i, double lam, double chargeDensityChange, double fitError) const;

    virtual ChargeDensity* GetExactChargeDensity () const;

protected:
    void Initialize(ChargeDensity* cd, double kT);  //Does on iteration to set up the exact charge density.

    WaveFunction*        itsWaveFunction;  //We don't own this guy!
private:
    Hamiltonian*         itsHamiltonian;

    ChargeDensity*       itsExactChargeDensity;
    ChargeDensity*       itsOldExactChargeDensity;
};

