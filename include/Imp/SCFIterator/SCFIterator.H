// File: SCFIteratorImplementation.H  Partial common implementation for an object that manages SCF convergence.
#ifndef _SCFIteratorImp_H_
#define _SCFIteratorImp_H_

#include <SCFIterator.H>

class SCFIteratorImp
    : public virtual SCFIterator
{
public:
    SCFIteratorImp(WaveFunction*, Hamiltonian*);
    virtual ~SCFIteratorImp();
    virtual bool Iterate(const SCFIterationParams& ipar);
    virtual void DisplayEnergies(int i, double lam, double chargeDensityChange, double fitError) const;

    virtual Exact_CD* GetExactChargeDensity () const;

protected:
    void Initialize(Exact_CD* cd);  //Does on iteration to set up the exact charge density.

    WaveFunction*        itsWaveFunction;  //We don't own this guy!
private:
    Hamiltonian*         itsHamiltonian;

    Exact_CD*       itsExactChargeDensity;
    Exact_CD*       itsOldExactChargeDensity;
};

#endif // _SCFIteratorImp_H_

