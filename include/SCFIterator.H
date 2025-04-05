// File: SCFIterator.H  Interface for an object that manages SCF convergence.


#include "types.H"

class SCFIterator
{
public:
    SCFIterator(WaveFunction*, Hamiltonian*,Exact_CD*);
    ~SCFIterator();
    bool Iterate(const SCFIterationParams& ipar);
    void DisplayEnergies(int i, double lam, double chargeDensityChange, double fitError,const Exact_CD* cd) const;
    void DisplayEigen   () const;

    Exact_CD* GetExactChargeDensity () const;

private:
    void Initialize(Exact_CD* cd);  //Does on iteration to set up the exact charge density.

    WaveFunction*   itsWaveFunction;  //We don't own this guy!
    Hamiltonian*    itsHamiltonian;
    Exact_CD*       itsExactChargeDensity;
    Exact_CD*       itsOldExactChargeDensity;
};

