// File: SCFIterator/SCFIterator.C  Interface for an object that manages SCF convergence.
export module qchem.SCFIterator;
import qchem.SCFIterator.Types;
export import qchem.SCFAccelerator;
export import qchem.WaveFunction;
import qchem.WaveFunction.SCF;
export import qchem.SCFParams;

export using qchem::EnergyBreakdown;
using qchem::Hamiltonian::Hamiltonian;
using qchem::WaveFunction::WaveFunction;
using qchem::WaveFunction::SCFWaveFunction;
using qchem::SCFAccelerators::SCFAccelerator;
using qchem::ChargeDensity::DM_CD;

export namespace qchem::SCFIterator
{

class SCFIterator
{
public:
    SCFIterator(const bs_t*, const ElectronConfiguration*, ::Hamiltonian*,SCFAccelerator*,DM_CD* cd=0);
    virtual ~SCFIterator();
    virtual bool Iterate(const SCFParams& ipar);
    // Direct energy minimization (GDM owns the loop): geodesic line search, no density mixing.
    void SetDirectMin(bool b) {itsDirectMin=b;}

    // SCFIterator drives the mutable SCFWaveFunction, but only ever hands clients the const
    // read view (they can query the converged state, never drive someone else's SCF loop).
    const class WaveFunction* GetWaveFunction() const {return itsWaveFunction;}
    EnergyBreakdown     GetEnergy() const;
    size_t              GetIterationCount() const {return itsIterationCount;}
    bool                Converged() const {return itsConverged;}
private:
    void Initialize(DM_CD* cd);  //Does on iteration to set up the exact charge density.
    DM_CD* DirectMinStep(double Ecur, double mergeTol); //one direct-min step (returns new density)
    bool itsDirectMin=false;
    void DisplayEnergies(int i, const EnergyBreakdown&,  double relax, double dE, double dCD, size_t idealVirial) const;
    void DisplayEigen   () const;

    //All owned, see destructor.
    ::Hamiltonian*    itsHamiltonian;
    SCFAccelerator* itsAccelerator;
    SCFWaveFunction*      itsWaveFunction;
    DM_CD*          itsCD;
    DM_CD*          itsOldCD;

    size_t          itsIterationCount;
    bool            itsConverged;
};

} //namespace


