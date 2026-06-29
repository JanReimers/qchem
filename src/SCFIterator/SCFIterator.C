// File: SCFIterator/SCFIterator.C  Interface for an object that manages SCF convergence.
module;
#include <memory>
export module qchem.SCFIterator;
import qchem.SCFIterator.Types;
export import qchem.SCFAccelerator;
export import qchem.WaveFunction;
import qchem.WaveFunction.SCF;
export import qchem.SCFParams;
export import qchem.ChargeDensity.Seed;   // SeedStrategy / MakeSeedDensity

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
    // The seed density is chosen by strategy (see ChargeDensity::SeedStrategy): Default resolves to
    // each path's present-day behaviour -- molecular -> CoreGuess, plane-wave -> Uniform.  \a st (the
    // structure) is needed only by the SAD seeds (atom Z + positions); null is fine otherwise.
    tSCFIterator(const tbs_t<T>*, const ElectronConfiguration*, ham_t*,acc_t*,
                 ChargeDensity::SeedStrategy seed=ChargeDensity::SeedStrategy::Default,
                 const Structure* st=nullptr);
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
    typedef std::shared_ptr<tDM_CD<T>> cd_t;   //!< std-managed density: no manual delete, no freed-address reuse
    void Initialize(tDM_CD<T>* cd);  //Does on iteration to set up the exact charge density.
    cd_t DirectMinStep(double Ecur, double mergeTol); //one direct-min step (returns new density)
    bool itsDirectMin=false;
    void DisplayEnergies(int i, const EnergyBreakdown&,  double relax, double dE, double dCD, size_t idealVirial) const;
    void DisplayEigen   () const;

    //Raw ptrs owned, see destructor; the charge densities are std-managed (cd_t).
    ham_t*          itsHamiltonian;
    acc_t*          itsAccelerator;
    scfwf_t*        itsWaveFunction;
    cd_t            itsCD;       //!< current charge density (shared_ptr: lifetime by std, no reuse)
    cd_t            itsOldCD;    //!< previous charge density

    size_t          itsIterationCount;
    bool            itsConverged;
};

using SCFIterator  = tSCFIterator<double>;
using cSCFIterator = tSCFIterator<dcmplx>;

} //namespace


