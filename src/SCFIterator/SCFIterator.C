// File: SCFIterator/SCFIterator.C  Interface for an object that manages SCF convergence.
export module qchem.SCFIterator;
import qchem.SCFIterator.Types;
export import qchem.SCFAccelerator;
export import qchem.WaveFunction;
import qchem.WaveFunction.SCF;
export import qchem.SCFParams;

export using qchem::EnergyBreakdown;
using qchem::ChargeDensity::tDM_CD;

export namespace qchem::SCFIterator
{

// Templated on the matrix element type T (rX/cX); SCFIterator is the <double> alias (atoms/
// molecules), cSCFIterator the <dcmplx> instantiation that drives single-k plane-wave DFT through
// the same assemble->diagonalize->fill->build-density loop.
template <class T> class tSCFIterator
{
    typedef qchem::Hamiltonian::tHamiltonian<T>  ham_t;
    typedef qchem::WaveFunction::tWaveFunction<T>    wf_t;
    typedef qchem::WaveFunction::tSCFWaveFunction<T> scfwf_t;
    typedef qchem::SCFAccelerators::tSCFAccelerator<T> acc_t;
public:
    tSCFIterator(const tbs_t<T>*, const ElectronConfiguration*, ham_t*,acc_t*,tDM_CD<T>* cd=0);
    virtual ~tSCFIterator();
    virtual bool Iterate(const SCFParams& ipar);
    // Direct energy minimization (GDM owns the loop): geodesic line search, no density mixing.
    void SetDirectMin(bool b) {itsDirectMin=b;}

    // SCFIterator drives the mutable SCFWaveFunction, but only ever hands clients the const
    // read view (they can query the converged state, never drive someone else's SCF loop).
    const wf_t* GetWaveFunction() const {return itsWaveFunction;}
    EnergyBreakdown     GetEnergy() const;
    size_t              GetIterationCount() const {return itsIterationCount;}
    bool                Converged() const {return itsConverged;}
private:
    void Initialize(tDM_CD<T>* cd);  //Does on iteration to set up the exact charge density.
    tDM_CD<T>* DirectMinStep(double Ecur, double mergeTol); //one direct-min step (returns new density)
    bool itsDirectMin=false;
    void DisplayEnergies(int i, const EnergyBreakdown&,  double relax, double dE, double dCD, size_t idealVirial) const;
    void DisplayEigen   () const;

    //All owned, see destructor.
    ham_t*          itsHamiltonian;
    acc_t*          itsAccelerator;
    scfwf_t*        itsWaveFunction;
    tDM_CD<T>*      itsCD;
    tDM_CD<T>*      itsOldCD;

    size_t          itsIterationCount;
    bool            itsConverged;
};

using SCFIterator  = tSCFIterator<double>;
using cSCFIterator = tSCFIterator<dcmplx>;

} //namespace


