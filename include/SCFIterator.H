// File: SCFIterator.H  Interface for an object that manages SCF convergence.

#include <omlfwd.H>
#include <Symmetry/fwd.H>
#include <BasisSet/fwd.H>
#include <ChargeDensity/fwd.H>
#include <Hamiltonian/fwd.H>
#include <SCFAccelerator/fwd.H>
#include <WaveFunction/fwd.H>
#include <SCFIteratorfwd.H>

class SCFIterator
{
public:
    SCFIterator(const BasisSet*, const ElectronConfiguration*, Hamiltonian*,SCFAccelerator*,DM_CD* cd=0);
    virtual ~SCFIterator();
    virtual bool Iterate(const SCFParams& ipar);

    const WaveFunction* GetWaveFunction() const {return itsWaveFunction;}
    EnergyBreakdown     GetEnergy() const;
    size_t              GetIterationCount() const {return itsIterationCount;}
private:
    void Initialize(DM_CD* cd);  //Does on iteration to set up the exact charge density.
    void DisplayEnergies(int i, const EnergyBreakdown&,  double relax, double dE, double dCD) const;
    void DisplayEigen   () const;

    //All owned, see destructor.
    Hamiltonian*    itsHamiltonian;
    SCFAccelerator* itsAccelerator;
    WaveFunction*   itsWaveFunction;  
    DM_CD*          itsCD;
    DM_CD*          itsOldCD;

    size_t          itsIterationCount;
};




