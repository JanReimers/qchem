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

//=== TUNING KNOBS AND POLICY -- PUBLIC (lifted out of the .Internal. modules 2026-08-09) ==========
// These are what a USER sets, so they belong on the public face: a caller choosing DIIS wants Nproj and
// the residual gates, and a caller building a ladder must be able to name its scheduling signal.  They
// lived beside their implementations only because the molecular path reached them through a json blob;
// the moment a FACADE had to name them (Step 4) that placement became the thing forcing an Internal
// import across a library boundary.  The concrete accelerators still own the ALGORITHMS -- only the
// vocabulary moved.

struct DIISParams
{
    size_t Nproj;  //Number of terms to keep for proections.
    double FDMax;  //DIIS STARTS extrapolating once the residual [F,D] < FDMax.  (Named for the residual it
                   //  gates on -- the DIIS error IS [F,D], NOT the energy; the old "EMax" name was a trap.)
    double FDMin;  //DIIS STOPS once [F,D] < FDMin (converged; hand back to plain diagonalization).
    double SVTol;  //DIIS bails out when the minimum singular value of B matrix is < SVTol;
};

struct GDMParams
{
    double FDMax;         //Engage the geodesic step once the residual [F,D] < FDMax (below that, run a
                          //  diagonalizing step).  Named for what it gates on -- [F,D], NOT the energy.
    double Trust=0.1;     //Trust radius: cap the largest geodesic rotation angle (radians) per step.
    //! Step-REJECTION policy (the \c RejectStep bail-out above).  On a rejected step the minimizer drops its
    //! CG history and multiplies the LIVE trust radius by \c TrustBackoff; once that falls below
    //! \c TrustMin it declares itself EXHAUSTED and forces a diagonalizing step, which is what makes the
    //! caller's fallback safe.  The budget is log(TrustMin/Trust)/log(TrustBackoff) retries -- 6 at these
    //! defaults from the GPW factory's Trust -- and an ACCEPTED step restores the full radius, so the
    //! backoff is per-step recovery rather than a ratchet.  Retries are not free (a full backtrack each), but
    //! a failure that is SCALE-INVARIANT is short-circuited by the caller long before the budget runs out.
    double TrustBackoff=0.25;
    double TrustMin=1e-4;
    //! PRECONDITIONER FLOOR on the diagonal orbital-Hessian denominator (eps_a - eps_i).
    //! The preconditioner divides the o-v gradient by that gap, which assumes it is POSITIVE -- true only at
    //! an aufbau point.  Two things go wrong without a floor, and the second is a correctness bug rather than
    //! a conditioning one:
    //!   * NEAR-DEGENERATE (gap -> 0): the step blows up along exactly the soft modes, which is the measured
    //!     GDM wander on diffuse bases (V1.24).
    //!   * NON-AUFBAU (gap < 0, i.e. an empty level below an occupied one -- a HOLE): the quotient FLIPS
    //!     SIGN, so that component of "steepest descent" points UPHILL, and the Polak-Ribiere denominator
    //!     Re<g,PG> can go negative and blow up beta with it.  A minimizer holding a hole -- which is the
    //!     normal state of a MOM or symmetry-broken run -- was therefore building an ascent direction from a
    //!     preconditioner that is not positive-definite.
    //! Clamping the denominator from below restores positive-definiteness, which is all a preconditioner has
    //! to be: it is a METRIC, not the Hessian, so a floored gap changes the step LENGTH along stiff modes and
    //! never the descent property.  0.1 Ha is well below a real HOMO-LUMO gap and well above the ties.
    double PCFloor=0.1;
};

//! Which convergence signal drives the TAIL hand-off (see the design notes above + doc/SCFStrategyPlan.md).
//! The scheduling signal is a LOOP-FACE policy, NOT a per-rung param -- the ladder already consumes both
//! signals (SetEnergy + GetError), so the choice lives here, once, rather than smeared across the rungs.
//!   Error       -- switch once the residual [F,D] < switchat.  Right for well-conditioned MOLECULES, where
//!                  [F,D] (the true orbital gradient) drives cleanly to 0.
//!   EnergyChange -- switch once |ΔE/E| < switchat.  Right for SOLIDS, where [F,D] is contaminated by
//!                  charge-transfer slosh / diffuse-mode ill-conditioning / giant-response spikes and may
//!                  never settle, while the total energy (a variational scalar) does.
enum class ScheduleSignal { Error, EnergyChange };

using SCFIrrepAccelerator  = tSCFIrrepAccelerator<double>;  using cSCFIrrepAccelerator = tSCFIrrepAccelerator<dcmplx>;
using SCFAccelerator       = tSCFAccelerator<double>;       using cSCFAccelerator      = tSCFAccelerator<dcmplx>;

} //namespace

