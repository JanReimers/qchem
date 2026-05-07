// File: SCFIterator/SCFIterator.C  Interface for an object that manages SCF convergence.
export module qchem.SCFIterator;
import qchem.SCFIterator.Types;
export import qchem.SCFAccelerator;
export import qchem.WaveFunction;
export import qchem.SCFParams;
export import qchem.BasisSet;

export using qchem::EnergyBreakdown;
using qchem::Hamiltonian::Hamiltonian;
using qchem::WaveFunction::WaveFunction;
using qchem::SCFAccelerators::SCFAccelerator;
using qchem::ChargeDensity::DM_CD;

export namespace qchem::SCFIterator
{

class SCFIterator
{
public:
    SCFIterator(const bs_t*, const ElectronConfiguration*, class Hamiltonian*,SCFAccelerator*,DM_CD* cd=0);
    virtual ~SCFIterator();
    virtual bool Iterate(const SCFParams& ipar);

    const class WaveFunction* GetWaveFunction() const {return itsWaveFunction;}
    EnergyBreakdown     GetEnergy() const;
    size_t              GetIterationCount() const {return itsIterationCount;}
    bool                Converged() const {return itsConverged;}
private:
    void Initialize(DM_CD* cd);  //Does on iteration to set up the exact charge density.
    void DisplayEnergies(int i, const EnergyBreakdown&,  double relax, double dE, double dCD) const;
    void DisplayEigen   () const;

    //All owned, see destructor.
    class Hamiltonian*    itsHamiltonian;
    SCFAccelerator* itsAccelerator;
    class WaveFunction*   itsWaveFunction;  
    DM_CD*          itsCD;
    DM_CD*          itsOldCD;

    size_t          itsIterationCount;
    bool            itsConverged;
};

} //namespace


