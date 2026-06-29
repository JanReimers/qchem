// File: SCFIterator/SCFIterator.C  Interface for an object that manages SCF convergence.
module;
#include <memory>
#include <functional>
export module qchem.SCFIterator;
import qchem.SCFIterator.Types;
export import qchem.SCFAccelerator;
export import qchem.WaveFunction;
import qchem.WaveFunction.SCF;
export import qchem.SCFParams;
export import qchem.ChargeDensity.Seed;   // SeedStrategy / MakeSeedDensity

export using qchem::EnergyBreakdown;
using qchem::ChargeDensity::tDM_CD;
using qchem::ChargeDensity::tChargeDensity;

export namespace qchem::SCFIterator
{

//! Per-iteration progress, handed to an Observer each SCF step so a client (a GUI,
//! a logger) can watch convergence live without owning the loop.  All real-valued
//! (energies/residuals are real for both the rX and cX paths).
struct SCFProgress
{
    size_t iteration;
    double energy;       //!< total energy E (hartree)
    double dE;           //!< |E_n - E_{n-1}|
    double commutator;   //!< [F,D] (the accelerator/DIIS error)
    double drho;         //!< relative charge-density change
};

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

    // Watch convergence live: the observer (if set) fires once per SCF iteration with
    // the current SCFProgress.  Read-only telemetry -- the observer must not drive the loop.
    using Observer = std::function<void(const SCFProgress&)>;
    void SetObserver(Observer obs) {itsObserver=std::move(obs);}

    // SCFIterator drives the mutable SCFWaveFunction, but only ever hands clients the const
    // read view (they can query the converged state, never drive someone else's SCF loop).
    const wf_t* GetWaveFunction() const {return itsWaveFunction;}
    EnergyBreakdown     GetEnergy() const;
    size_t              GetIterationCount() const {return itsIterationCount;}
    bool                Converged() const {return itsConverged;}
private:
    typedef std::shared_ptr<tDM_CD<T>> cd_t;   //!< std-managed WORKING density (matrix-backed); no manual delete
    //! Seed the SCF: build the iteration-0 Fock from \a seed (a DFT-face tChargeDensity -- may be a fit, e.g.
    //! the SAD seed; or null for the core guess), diagonalize, and take the first real (matrix) density.
    //! \a bs / \a st are forwarded for the HF/DHF bootstrap (build a DFT sibling when the seed has no matrix
    //! but the Hamiltonian needs one -- see project_numericcd_refactor); null is fine for the core guess.
    void Initialize(tChargeDensity<T>* seed, const tbs_t<T>* bs, const Structure* st);
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
    Observer        itsObserver;   //!< optional live-progress sink (default empty)
};

using SCFIterator  = tSCFIterator<double>;
using cSCFIterator = tSCFIterator<dcmplx>;

} //namespace


