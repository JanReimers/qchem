// File: SCFIterator.H  Interface for an object that manages SCF convergence.


#include "types.H"

class SCFIterator
{
public:
    SCFIterator(WaveFunction*, Hamiltonian*,DM_CD* cd=0);
    ~SCFIterator();
    bool Iterate(const SCFIterationParams& ipar);
    void DisplayEnergies(int i, const EnergyBreakdown&,  double relax, double chargeDensityChange, double fitError) const;
    void DisplayEigen   () const;

    DM_CD* GetExactChargeDensity () const;

private:
    void Initialize(DM_CD* cd);  //Does on iteration to set up the exact charge density.

    WaveFunction* itsWaveFunction;  //We don't own this guy!
    Hamiltonian*  itsHamiltonian;
    DM_CD*        itsCD;
    DM_CD*        itsOldCD;
};

