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
import qchem.LASolver;   // qchem::Ortho (the basis-overlap orthogonalisation knob, forwarded to the WF)
import qchem.BasisSet.Fit_IBS;   // BasisSet::FIT_SF_ABS<T> (the G-space fit basis for Kerker rho-tilde extraction)
export import qchem.ChargeDensity.DensityMixer;   // tDensityMixer<T> (the density-face of SCF convergence)
import qchem.SCFIterator.LoopDriver;   // tLoopDriver<T> + Fixed/DirectMin concretes (the loop-face seam)

export using qchem::EnergyBreakdown;
using qchem::ChargeDensity::tDM_CD;
using qchem::ChargeDensity::tChargeDensity;

export namespace qchem::SCFIterator
{

//! Process-wide diagnostic toggle (default OFF), mirroring \c qchem::ReportOverlapConditioning and
//! \c qchem::Hamiltonian::ReportGridCharge.  When true, the verbose per-iteration SCF line appends the
//! frontier spectrum -- ε_HOMO, ε_LUMO and the gap ε_LUMO−ε_HOMO (Ha) -- from the current orbital energy
//! levels.  The band-gap instrument for the NaF Γ-instability (doc/GPWPlan §0b″): a near-degenerate
//! HOMO/LUMO (gap → 0) is the giant-response hypothesis; watching the gap per iteration lets the spurious
//! level be seen diving across the Fermi edge one step before each energy spike.  Flip it around
//! \c Iterate and reset it (it is a static, so it leaks between tests otherwise).
bool& ReportBandGap();

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
    // \a basisOrtho / \a basisOrthoTol select how the orbital-overlap S is orthogonalised: Cholesky
    // (default; needs S positive-definite) or Eigen/SVD with a cutoff that drops near-null eigen/singular
    // values -- canonical orthogonalisation for a linearly-dependent basis (e.g. GPW's diffuse Gaussians
    // on a dense lattice, whose Bloch overlap goes singular).  \a basisOrthoTol<=0 keeps all.
    tSCFIterator(const tbs_t<T>*, const ElectronConfiguration*, ham_t*,acc_t*,
                 ChargeDensity::SeedStrategy seed=ChargeDensity::SeedStrategy::Default,
                 const Structure* st=nullptr,
                 qchem::Ortho basisOrtho=qchem::Cholesky, double basisOrthoTol=0.0);
    virtual ~tSCFIterator();
    virtual bool Iterate(const SCFParams& ipar);

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
    cd_t DirectMinStep(double Ecur, double mergeTol); //one direct-min step (returns new density; used by DirectMinDriver)

    void DisplayEnergies(int i, const EnergyBreakdown&,  double relax, double dE, double dCD, size_t idealVirial) const;
    void DisplayEigen   () const;

    //Raw ptrs owned, see destructor; the charge densities are std-managed (cd_t).
    ham_t*          itsHamiltonian;
    acc_t*          itsAccelerator;
    scfwf_t*        itsWaveFunction;
    cd_t            itsCD;       //!< current charge density (shared_ptr: lifetime by std, no reuse)
    cd_t            itsOldCD;    //!< previous charge density
    //! The SCF density lineage (ChargeDensity::Lineage).  SetWorkingCD makes each new itsCD the head, so a
    //! superseded density (itsOldCD, a stale copy) reports isActive()==false and trips the Hamiltonian's
    //! assert if reused.  Created in Init; one per SCF run.
    qchem::ChargeDensity::LineagePtr itsLineage;
    //! Assign the new working density AND make it the head of itsLineage (so the previous one goes inactive).
    //! Every itsCD (re)assignment in the working loop goes through here.
    void SetWorkingCD(cd_t cd) {itsCD=std::move(cd); if (itsCD && itsLineage) itsCD->JoinLineage(itsLineage);}

    size_t          itsIterationCount;
    bool            itsConverged;
    Observer        itsObserver;   //!< optional live-progress sink (default empty)

    // Density-face state.  The mixer (Linear / Kerker; see qchem.ChargeDensity.DensityMixer) is built per-run
    // from SCFParams at the top of Iterate and owns the mixing policy + state (relax, the Kerker ρ̃, ...).
    const tbs_t<T>*  itsBS = nullptr;         //!< orbital basis (for the Kerker G-space fit basis) -- from the ctor
    //! A PERSISTENT copy of the periodic cell (reciprocal lattice + volume for Kerker).  The ctor's raw \c st
    //! comes from a temporary (\c Lattice_3D::GetStructure returns a fresh \c make_shared), so it dangles by the
    //! time \c Iterate runs -- we deep-copy it here (periodic path only) so the Kerker mixer has a live cell.
    std::shared_ptr<const Structure> itsKerkerCell;
    std::unique_ptr<qchem::ChargeDensity::tDensityMixer<T>> itsMixer;  //!< the density-face concrete for this run
    // The two loop-face concretes (stateless).  Iterate selects one per macro-iteration by the accelerator's
    // WantsLineSearch() and dispatches Step() -- virtual dispatch in place of the old mode `if`.
    FixedPointDriver<T> itsFixedDriver;
    DirectMinDriver<T>  itsDirectDriver;
};

using SCFIterator  = tSCFIterator<double>;
using cSCFIterator = tSCFIterator<dcmplx>;

} //namespace


