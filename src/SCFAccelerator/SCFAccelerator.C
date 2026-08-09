// FIle: SCFAccelerator.C  Interface for an SCF accelerator alogrithm
module;
#include <iosfwd>
export module qchem.SCFAccelerator;
export import qchem.Symmetry.Irrep;
export import qchem.LASolver;

export namespace qchem::SCFAccelerators
{

// Templated on the matrix element type T (rX/cX convention).  For T=double, hmat_t<double> IS
// smat_t<double>=rsmat_t, so the existing real accelerators (DIIS/GDM/Ladder/Null) -- which bind the
// <double> aliases below -- are unchanged.  The plane-wave path uses tSCFAcceleratorNull<dcmplx>.
template <class T> class tSCFIrrepAccelerator
{
public:
    virtual ~tSCFIrrepAccelerator() {};
    // Feed the current (AO) Fock matrix and the (orthonormal-basis) density matrix.
    virtual void UseFD(const hmat_t<T>& F, const hmat_t<T>& DPrime)=0;
    // Produce the next set of orbital coefficients (U, U', e), as SolveOrtho returns.
    // A Fock-extrapolator (DIIS) extrapolates F then diagonalizes; a direct minimizer
    // (GDM) rotates the current orbitals along the Grassmann manifold instead.
    virtual typename LASolver<T>::UUd_t NextOrbitals()=0;

    // Direct-minimization line-search hooks (only GDM-style minimizers implement these).
    //   ComputeStep(): compute the search direction/geodesic for the current Fock without
    //     moving the orbitals.  Returns false if the accelerator does not support a line
    //     search, or wants the caller to fall back to NextOrbitals() (e.g. its seed step).
    //   OrbitalsAt(t,commit): orbitals at geodesic fraction t; commit=false is a pure trial
    //     (for evaluating the energy), commit=true takes the step.
    virtual bool ComputeStep() {return false;}
    virtual typename LASolver<T>::UUd_t OrbitalsAt(double t, bool commit) {return NextOrbitals();}
    //! STEP REJECTION (the bail-out callback; see the aggregate's twin for the full contract).
    virtual bool RejectStep() {return false;}
};

template <class T> class tSCFAccelerator
{
public:
    virtual ~tSCFAccelerator() {};
    virtual tSCFIrrepAccelerator<T>* Create(const LASolver<T>*,const Irrep&, int occ)=0;
    virtual bool CalculateProjections()=0;
    virtual void ShowLabels     (std::ostream&) const=0;   // DEPRECATED (doc/GPWPlan1.md item 2): the free-form
    virtual void ShowConvergence(std::ostream&) const=0;   //   label/value blob, superseded by the compact
                                                           //   Tag()/Count() accel column; retire once it proves out.
    virtual double GetError() const=0;
    //! Compact self-identifier for the per-iteration \c accel column (doc/GPWPlan1.md item 2): a short tag
    //! ("Null", "DIIS", "GDM", ...) plus a \c Count() (DIIS projection depth; 0 when not meaningful).  The
    //! ladder delegates both to its active rung.  Distinct from \c GetError() (the [F,D] the column also shows).
    virtual const char* Tag  () const=0;
    virtual int         Count() const {return 0;}
    // Has this accelerator run out of steam (ladder hand-off signal)?  Default: never.
    virtual bool Exhausted() const {return false;}
    // The SCF iterator reports the current total energy each macro-iteration.  The ladder
    // uses the energy change to decide hand-offs (see SCFAcceleratorLadder); others ignore it.
    virtual void SetEnergy(double E) {}
    // Should the SCF iterator run its direct-minimization loop (geodesic line search, no
    // density mixing)?  True for a GDM-style minimizer; the ladder reports its active rung's.
    virtual bool WantsLineSearch() const {return false;}
    // Is a geodesic step actually READY this iteration (seeded + residual under its engage threshold)?  The
    // iterator runs the direct-min driver only when WantsLineSearch() AND CanLineSearch(); otherwise it does a
    // MIXED fixed-point step -- so a not-yet-ready minimizer degrades to a stable mixed diagonalize, never an
    // UNMIXED one (which runs away on ill-conditioned systems).  Default false: non-minimizers never line-search.
    virtual bool CanLineSearch() const {return false;}

    //! STEP REJECTION -- the BAIL-OUT callback.  An accelerator PROPOSES a step; only the CALLER can judge
    //! it, because only the caller can evaluate the energy (the accelerator sits BELOW the wavefunction /
    //! Hamiltonian in the library DAG, which is exactly why the step body lives in the iterator at all).
    //! So acceptance is the caller's decision, and this is how the decision is reported back: the iterator
    //! calls \c RejectStep() when the proposed step was NOT acceptable -- for a direct minimizer, a line
    //! search that found no descent at ANY backtrack.
    //!
    //! The accelerator responds by DEGRADING its model (a CG minimizer drops its conjugacy history and
    //! shrinks its trust radius) and answers ONE question: is a RETRY worth attempting?
    //!   true  -- a fresh, more conservative direction is ready; the caller may line-search again.
    //!   false -- EXHAUSTED.  The accelerator MUST then leave itself in a state where \c ComputeStep()
    //!            returns false AND \c CanLineSearch() is false, so the caller's fallback (diagonalize,
    //!            then density-mix) is SAFE.  That guarantee is the whole point of the contract: without
    //!            it the caller cannot fall back at all, because \c NextOrbitals() on a minimizer holding
    //!            a live step TAKES that step (GDM: \c OrbitalsAt(1.0,true)) -- the very step just
    //!            rejected, at full length.
    //!
    //! Default false = "I have no step to reject and no retry to offer".  That is correct for every
    //! non-minimizer: their \c ComputeStep() is false to begin with, so the caller's fallback was already
    //! safe and a rejection simply degrades to it.
    virtual bool RejectStep() {return false;}
};

using SCFIrrepAccelerator  = tSCFIrrepAccelerator<double>;  using cSCFIrrepAccelerator = tSCFIrrepAccelerator<dcmplx>;
using SCFAccelerator       = tSCFAccelerator<double>;       using cSCFAccelerator      = tSCFAccelerator<dcmplx>;

} //namespace

