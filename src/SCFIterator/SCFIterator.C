// File: SCFIterator/SCFIterator.C  Interface for an object that manages SCF convergence.
export module qchem.SCFIterator;
export import qchem.SCFAccelerator;
export import qchem.WaveFunction;
export import qchem.SCFParams;
export import qchem.BasisSet;

export class SCFIterator
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




