// File: Calculation/SolidCalculation.C
//
// qchem::SolidCalculation -- the production front door for a PERIODIC (GPW / plane-wave) SCF.  The
// solid peer of qchem::MoleculeCalculation (molecules) and qchem::AtomCalculation (atoms): it owns the whole
// object graph -- basis, electron configuration, Hamiltonian, accelerator, iterator -- runs the
// canonical assemble->converge recipe, and exposes high-level questions rather than the graph.
//
// WHY IT EXISTS, beyond symmetry with the other two (doc/CleanupCandidates.md Step 4).  The user's
// reason is the operative one: there needs to be something SPECIFIC AND NAMED to point at when a
// decision belongs ABOVE the SCF iterator.  Those decisions were previously spread through a driver
// function in an integration test, which meant that in practice they had no home at all -- D6 moved
// the Becke recipe out of that file, V1.26 moved the grid selector, and each time the question "where
// does this belong?" had no better answer than "not there".  This class is that answer.
//
// DECISIONS THIS CLASS OWNS (the point of it):
//   - WHICH XC QUADRATURE the run uses.  It gathers the run's sharpness (basis alpha_max + local-PP
//     alpha_pp, both off ABSTRACT capability faces) and hands it to qcMesh::ResolveXCMesh, which costs
//     the two grids and picks (V1.26/V2.4).  A caller may still pin cellKind and be obeyed.
//   - WHICH ACCELERATOR, by policy name rather than by constructing internals.
//   - The SPIN BOOKKEEPING: multiplicity -> (nUp,nDown) -> Crystal_EC, with the parity check that
//     catches "singlet with an odd electron count" before it silently empties a channel.
//
// WHY IT COULD NOT BE WRITTEN UNTIL NOW.  Both the solid Hamiltonian and the complex accelerator were
// reachable ONLY through .Internal. modules, which per CLAUDE.md may be imported across a library
// boundary by unit tests and nothing else.  That -- not carelessness -- is why the driver lived in a
// test file: it was the one place the cheat was legal.  Step 4 opened the two public doors first
// (qchem.Hamiltonian.Factory's cHamiltonian overload; qchem.SCFAccelerator.Factory's typed-options
// overload), so this file imports ZERO internals, exactly as the molecular facades already manage.
module;
#include "forward.H"   // RunDiagnosticsTests, ResponseFacadeTests -- the unit-test friends (CLAUDE.md: tests may cheat)
#include <functional>  // the order-parameter probe
#include <memory>
#include <optional>
#include <vector>     // the anneal schedule
#include <vector>
#include <string>
#include <utility>
export module qchem.SolidCalculation;

import qchem.Lattice_3D;                      // Lattice_3D (crystal + BZ mesh)
import qchem.Structure;                       // Structure
import qchem.ScalarFunction;                  // ScalarFunction<double> (the rho(r) face)
import qchem.BasisSet;                        // Real_BS (orbital source), Complex_BS (the Bloch basis)
import qchem.BasisSet.Lattice.BasisSet;    // GPWFactory, GPWParams, MakeRasterLadder, CellImages
export import qchem.Hamiltonian.HubbardEstimator;   // ACBN0Estimate (EstimateHubbardU), HubbardEstimate (the loop)
export import qchem.Response;                  // ChannelResponse / ResponseFailure (IndependentResponse, LinearResponsePlan R0)
export import qchem.Response.Solver;           // SelfConsistentResponse / QMeshResponse + KrylovParams (HubbardLinearResponse)
import qchem.SCFAccelerator.Factory;          // Type, SolidAcceleratorOptions, the typed solid door
import qchem.SCFIterator;                     // SolidSCFIterator, SCFParams, IterationRow, EnergyBreakdown
import qchem.ChargeDensity;                   // cDM_CD
import qchem.ChargeDensity.Seed;              // SeedStrategy
import qchem.Mesh;                            // qcMesh::MeshParams / UnitCellKind
import qchem.LASolver;                        // qchem::Ortho
import qchem.Outcome;                         // Outcome<T,E> -- the fallible-call vocabulary (N1/T1)
export import qchem.SolidState;               // RestartRefusal (Restart's failure), the saved-state layer (CK-1)
export import qchem.SolidSetup;               // SolidCalcOptions + the resolved SolidSetup a run is built from (step 6)
export import qchem.PotentialEnergySurface;   // the face SolidSession implements (step 8)
export import qchem.SCFFailure;               // SCFFailure: why an attempt has no answer (shared with the molecular facade)
export import qchem.Calculation;              // the abstract base of the three facades

export namespace qchem
{

//! \brief WHAT THE RUN DID, iteration by iteration, in the two quantities the OUTCOME DETECTORS judge
//! (doc/OpenWork.md N1/T3-T4).  Measured by the facade itself, so a library user gets the detectors the
//! MnO campaign had to write into its own test file.
//!
//! THE ORGANISING PRINCIPLE (N1): you cannot enumerate bad COMBINATIONS of knobs, but you can detect bad
//! OUTCOMES.  These two instruments catch combinations nobody has thought of yet, which is why they run
//! ALWAYS rather than behind a diagnostic flag.
//!
//! THE TWO INSTRUMENTS, and why exactly these:
//!   - **The INTEGRATED site moment** \f$\max_A|\mu_A|\f$ over the Becke basins.  INTEGRATED, never a
//!     point sample of \f$m(\mathbf r)\f$: a point probe cannot tell an AFM pair \f$(+m,-m)\f$ from a
//!     disproportionated one \f$(0,-2m)\f$ -- measured 2026-08-08, a run reporting "order SURVIVED,
//!     m_stag=0.29" actually sat at sites \f$(+0.0006,-0.579)\f$ -- and its value depends on WHICH
//!     direction you step away from the nucleus.  It is a spin DENSITY, not a moment.
//!   - **The Hartree term** \f$E_{ee}\f$.  A low-G charge runaway roughly DOUBLES it while the total
//!     stays a plausible negative number: measured 2026-08-25 across four MnO collapse arms it read
//!     13.48 (healthy) / 29.0 / 35.1 / 13.6, ranking them correctly WITHOUT having been designed for
//!     any of them -- and the 13.6 row is the one that matters, a run whose moment died with NO slosh.
//!     Two failures with the same symptom and different mechanisms, which is exactly why a detector on
//!     the moment alone is not enough.
//!
//! BOTH ARE FREE.  The moment is a block sum over a raster the XC term has already built for this
//! iteration's density serial; \f$E_{ee}\f$ is already in every iteration's breakdown.  Neither adds a
//! matrix build, so there is nothing to gate behind a flag.
class RunDiagnostics
{
public:
    //! \name The ORDER (magnetic) channel
    //!@{
    //! Did this run have any order to lose?  False on an unpolarized run, on a run whose XC quadrature
    //! owns no site basins (a uniform mesh -- \c SiteMoments correctly returns empty and the detectors
    //! SKIP rather than guess), and on a genuinely closed-shell one.
    bool   HasOrder() const;
    //! \brief The YARDSTICK: the largest integrated site moment this run ever carried -- the RAW SEED's,
    //! or any iteration's, whichever is bigger.
    //!
    //! WHY THE HIGH-WATER MARK AND NOT THE SEED ALONE.  Both ends of the run mislead.  Iteration 1 is a
    //! terrible baseline for a quantity that GROWS before it dies (MnO reads 0.0046 at iteration 1,
    //! peaks at 0.106 by iteration 7 as the exchange splitting builds it back, then decays to 7e-5 --
    //! judged against iteration 1 that is "survived"; judged against the peak it is a 1400x collapse).
    //! And the seed alone is not enough either, for the mirror-image reason: the density available at
    //! construction is ALREADY post-iteration-0 (the wave function builds a Fock from the seed and fills
    //! it), so a run whose order dies in the very first fill would show no baseline at all -- measured
    //! on Na2, a seed staggered at +/-1 e reads +/-0.07 e one fill later.  Hence: the raw seed, measured
    //! before it is consumed, OR any later peak.
    double OrderPeak() const;
    double OrderFinal() const;   //!< the last iteration's integrated moment
    //! \brief Did the order RISE AND THEN DIE -- fall below 1% of \c OrderPeak and stay there?
    //! The same rule the MnO test's post-mortem used, now where every caller gets it.
    bool   OrderCollapsed() const;
    //! \brief The first STEP from which the order stayed dead; 0 when it never died.
    //! A STEP, not an SCF iteration number: the trajectory accumulates across an annealed schedule's
    //! stages (and across a re-\c Converge), because an order that died in stage 1 must not become
    //! invisible just because stage 2 got a fresh iterator and restarted its counter.
    size_t OrderDiedAt() const;
    //!@}

    //! \name The CHARGE channel
    //!@{
    bool   HasHartree() const;    //!< at least two iterations were observed (one point is not a trend)
    double HartreeFloor() const;  //!< the run's own healthy level: the MINIMUM \f$E_{ee}\f$ it visited
    double HartreePeak()  const;  //!< the largest \f$E_{ee}\f$ it visited
    //! \brief Did the Hartree term RUN AWAY -- did the run END above \c kSloshFactor times its own floor?
    //!
    //! SELF-CALIBRATING, and it has to be: \f$E_{ee}\f$ is extensive and has no universal scale (13.5 Ha
    //! is healthy for the MnO cell and would be a catastrophe for Si), so an absolute threshold would be
    //! an MnO constant wearing a library's clothes.
    //!
    //! \note IT IS A REFINEMENT OF A NON-CONVERGENT OUTCOME, NOT A VERDICT ON A CONVERGED ONE, and it
    //! says so ITSELF -- a converged run answers \c false whatever the ratio, rather than leaving the
    //! caller to remember the rule.  A converged density is stationary by definition, so there is
    //! nothing left sloshing; whereas a run that legitimately RESTRUCTURES on its way to the answer can
    //! raise \c Eee a long way and be perfectly healthy (measured: Na2 in a box, 0.114 -> 0.195 Ha =
    //! 1.71x, converged and correct).  Convicting such a run would trade a missed diagnosis for a WRONG
    //! one, which is much the worse error.
    bool   ChargeSloshed() const;
    double SloshRatio()    const; //!< the LAST \f$E_{ee}\f$ over \c HartreeFloor (1 if it never moved)
    //!@}

    //! \brief The one-line post-mortem -- what the MnO test used to print for itself.  It always says
    //! something: a run with no basins reports that it could not measure, which is a different fact from
    //! "the order was zero" and the two must not be allowed to read alike.
    std::string Summary() const;

    //! \name The rules, named and in one place so they can be argued with
    //!@{
    //! A quantity below 1% of its own high-water mark, and staying there, is DEAD.  Shared by the
    //! trajectory detector and the postcondition, so "collapsed" means one thing in this class.
    static constexpr double kCollapseFraction = 0.01;
    //! Below this many electrons the run never carried order, so there is nothing for a collapse to
    //! contradict.  It is the guard that keeps a closed-shell imposed run -- the common case -- silent.
    static constexpr double kOrderFloor = 0.05;
    //! \brief Ending above this multiple of the run's OWN \f$E_{ee}\f$ floor is a charge runaway.
    //!
    //! CALIBRATED, not guessed (doc/OpenWork.md T3).  Healthy converged runs measured 2026-08-26:
    //! MnO AFM-II imposed **1.05** (Eee 12.51 -> 13.18 Ha, after an early transient to 15.23 -- which
    //! is why the ratio is taken at the END and not at the peak), Na2-in-a-box **1.71** (a genuine
    //! restructuring from diffuse atoms into a bond).  Against that, the two MnO collapse arms measure
    //! **2.48** (flat Kerker, `MNO_KERKER_G0=0.01`: Eee 14.17 -> 35.10 Ha) and **2.01** (the banked
    //! aufbau/shared-mu recipe with `QCHEM_XC_DM_SOURCE=1`: 14.42 -> 29.00 Ha).  The flat-Kerker arm is
    //! worth looking at: after iteration 9 its Eee is a clean PERIOD-2 LIMIT CYCLE (35.210, 35.098,
    //! 35.210, ...) for seventy iterations, which is charge sloshing in the most literal sense
    //! available.  1.5 sits above every healthy run measured and below both collapses -- and because the
    //! detector is only ever consulted on a run that ALREADY failed, being wrong here costs a label,
    //! never a good answer.
    static constexpr double kSloshFactor = 1.5;
    //!@}

private:
    friend class SolidCalculation;
    //! The detector RULES are pure functions of the two trajectories below, so they are unit-tested on
    //! SYNTHETIC ones (src/Calculation/tests/RunDiagnostics.C).  Before that they were reachable only by
    //! coaxing a full SCF into the outcome under test -- which made a LOGIC test hostage to whether a
    //! marginal run converged inside its iteration cap, and it was silently disarmed twice that way.
    friend ::RunDiagnosticsTests;
    bool                itsConverged = false; //!< the last attempt's verdict -- ChargeSloshed() needs it
    double              itsSeedOrder = 0.0;   //!< the RAW seed's max_A |mu_A|, taken before Init consumes it
    bool                itsHasBasins = false; //!< the XC quadrature owns an atom-centred partition at all
    std::vector<double> itsOrder;             //!< per iteration, max_A |mu_A| (integrated)
    std::vector<double> itsEee;               //!< per iteration, the Hartree term
};

//! \brief ONE STAGE of an ANNEALED convergence: its own SCF parameters and its own accelerator.
//!
//! A descending-kT schedule with density + MOM continuation between stages is a DECISION ABOVE THE SCF
//! ITERATOR, which is what this class exists to own -- it lived in an integration-test driver instead,
//! which is why an MnO recipe was spread across four thousand lines of that file.
//!
//! Each stage gets a FRESH Hamiltonian and accelerator: a kT change must not carry stale DIIS history
//! across the re-seed.  What DOES carry across is the density, and -- under \c momFromSeed -- the occupied
//! SUBSPACE, so the configuration a hot stage settled on survives into the colder one.
struct SCFStage
{
    SCFParams             params;        //!< this stage's parameters; \c SmearingkT is its temperature
    SCFAccelerators::Type accelerator;   //!< per stage, so a recipe can anneal on DIIS and finish on GDM
};

//! \brief A periodic SCF: attempt it, and get back either the ANSWERS or the reason there are none.
class SolidCalculation : public virtual Calculation
{
    struct Imp;   //!< fwd: \c Converged holds one (the definition lives in the .C, as always)
public:
    //! \brief PROOF that the SCF converged, and the only place the answers live.
    //!
    //! WHY A SEPARATE TYPE rather than \c Fallible<double> on each accessor (user, 2026-08-25).
    //! "Did it converge?" is ONE fact governing a whole FAMILY of queries; wrapping each accessor repeats
    //! it and asks the caller to check N times.  Worse, a Fallible-style accessor ASSERTS on misuse, and
    //! an assert is compiled out under NDEBUG -- exactly how the imposed-mesh site-blocks defect stayed
    //! invisible in Release.  Here there is no check to forget and nothing to compile out: \b the
    //! \b accessors \b do \b not \b exist \b without \b the \b proof.  (CLAUDE.md: give capabilities only
    //! to types that have them.)
    //!
    //! \warning NON-OWNING, and invalidated by the next \c Converge() -- the same contract \c Density()
    //! has always carried.  Take it fresh from each attempt; do not store one across a re-run.
    class Converged
    {
    public:
        double                 Energy()         const;   //!< total energy (hartree)
        qchem::EnergyBreakdown EnergyTerms()    const;
        double                 TotalCharge()    const;   //!< the electron count (a charge-conservation check)
        size_t                 IterationCount() const;
        //! rho(r) of the converged state (see the class \warning for its lifetime).
        const ScalarFunction<double>& Density() const;
        //! \brief \f$m(r)=\rho_\uparrow-\rho_\downarrow\f$, or NULL on an unpolarized run.
        //!
        //! WHY (ρ, m) AND NOT (ρ↑, ρ↓).  They are the same information -- \f$\rho_\uparrow=(\rho+m)/2\f$ --
        //! but they are not equally good to hand out.  Charge and spin are the channels that BEHAVE
        //! differently: the Hartree restoring force acts on ρ alone, which is why Kerker preconditions the
        //! charge channel and has no business on the spin one (doc/OpenWork.md N3).  Exposing (↑,↓) invites
        //! exactly the channel-blind handling that costs magnetic runs their basin; exposing (ρ, m) makes
        //! the physical split the obvious one.  It is also what a magnetic diagnostic actually wants --
        //! \f$m(r)\f$ near a site, in one call rather than a cross-cast and a subtraction.
        const ScalarFunction<double>* SpinDensity() const;
        //! \brief The converged density OBJECT, for a caller that must CROSS-CAST it for a capability --
        //! \c FourierDensity to take a G-space component, \c cSpinResolved_CD to reach the channels.
        //!
        //! \c Density()/\c SpinDensity() are the convenient views and cover most consumers; this is the
        //! escape hatch for a campaign diagnostic that needs to ASK THE DENSITY WHAT IT CAN DO, which is
        //! how capabilities are reached everywhere else in this tree.  Same lifetime as the rest.
        const qchem::ChargeDensity::cDM_CD& DensityMatrix() const;
    private:
        friend class SolidCalculation;
        explicit Converged(const Imp* imp) : itsImp(imp) {}
        const Imp* itsImp;
    };

    using sf_t     = ScalarFunction<double>;
    using Observer = qchem::SCFIterator::SolidSCFIterator::Observer;

    //! \brief Build the whole graph over a RESOLVED setup (doc/CalcNetworkReview.md §7, step 6) -- the geometry point at the
    //! setup's reference lattice -- and converge once with \a params and the setup's accelerator.  The setup is
    //! \c SolidSetup::Resolve(lattice, molecular basis | recipe, options, accelerator options); the only door (user ruling
    //! 2026-10-30: clean replace, no (lattice, basis, options) constructors).
    //! \a onIteration, when given, observes THIS constructor's SCF from its first iteration -- the constructor converges, so an
    //! observer attached afterwards (\c OnIteration) has already missed stage 0, the stage an annealed recipe's trace most wants;
    //! it stays attached for later \c Converge calls.  A constructor argument, not a setup field: an observer is not a resolved
    //! decision (the molecular \c MoleculeCalculation takes it the same way).
    explicit SolidCalculation(std::shared_ptr<const SolidSetup> setup, const SCFParams& params = {}, Observer onIteration = nullptr);
    //! \brief The ANNEALED constructor: stage 0 on its own accelerator, then the rest of \a schedule, each stage seeded from the one
    //! before it.  An empty schedule throws.  \a onIteration observes every stage, from stage 0's first iteration.
    SolidCalculation(std::shared_ptr<const SolidSetup> setup, const std::vector<SCFStage>& schedule, Observer onIteration = nullptr);
    ~SolidCalculation();

    //! \brief START FROM A SAVED STATE (CK-1; doc/OpenWork.md §2 "SCF checkpoint/restart"): build the graph for
    //! \a setup exactly as the constructor does, check the file's fingerprint against the BUILT basis, seed the
    //! SCF with the saved density matrices instead of \c opts.seed, and converge with \a params.
    //!
    //! EXACT RESUME when nothing in the fingerprint differs; WARM START when only the warm class does (+U, the
    //! grids, the multiplicity, kT, the recipe -- the SCU loop and grid continuation); both say so on the
    //! console, with every difference named.  FAILS -- REFUSED, the graph discarded -- on an absent/foreign
    //! file or any difference that makes the saved D meaningless here (structure, species, basis, k-mesh,
    //! electron count, spin group, functional).  A refusal is a value because the caller's fallback (a fresh
    //! seed) is a decision only the caller can make.
    //!
    //! No mixer or accelerator history is carried: a restart is a fresh stage.  With \c opts.momFromSeed the
    //! MOM reference is the restored density's own first fill, as for any seed.  A \c Ladder accelerator STARTS
    //! ON ITS GDM RUNG for the restart's stage when \a params.SmearingkT is 0 (\c SolidAcceleratorOptions::ladderStartOnGDM):
    //! the DIIS rung's mixing would step a converged density away before the first energy is read.  A smeared restart
    //! starts on DIIS (GDM would engage on a restored integer D' and never apply the smearing).
    //! The setup's accelerator (\c opts.accelerator) drives the restored stage.  \a onIteration observes the restored stage from
    //! its first iteration (a restart's first iterate is what an exact resume is judged by).
    static Outcome<std::unique_ptr<SolidCalculation>, RestartRefusal>
        Restart(const std::string& path, std::shared_ptr<const SolidSetup> setup, const SCFParams& params = {},
                Observer onIteration = nullptr);

    //! \name BUILT, NOT CONVERGED (2026-10-31, the deck's \c deck::Build: lane G binds the GUI to the deck)
    //! For a caller that must HOLD the calculation before its SCF starts -- a GUI hands it to a worker thread and watches it.
    //! The graph is built exactly as the constructors build it; the caller runs the SCF.  ONE path: the annealed constructor
    //! is \c Prepare + \c Converge(schedule), \c Restart is \c PrepareRestart + \c Converge(params).
    //!@{
    //! \a schedule's stage 0 standing (its parameters, accelerator and seed) and no SCF; run it with \c Converge(schedule).
    //! \a onIteration observes from that first iteration.  An empty schedule throws.
    static std::unique_ptr<SolidCalculation> Prepare(std::shared_ptr<const SolidSetup> setup, const std::vector<SCFStage>& schedule,
                                                     Observer onIteration = nullptr);
    //! \c Restart up to its SCF: the fingerprint judged, the saved density restored, the restored stage standing for \a params;
    //! run it with \c Converge(params).  Refused exactly when \c Restart is; an exact resume's first convergence prints the resumed
    //! energy beside the saved one, whoever runs it.
    static Outcome<std::unique_ptr<SolidCalculation>, RestartRefusal>
        PrepareRestart(const std::string& path, std::shared_ptr<const SolidSetup> setup, const SCFParams& params = {},
                       Observer onIteration = nullptr);
    //!@}

    //! \brief Write the LAST ITERATE's state -- every (k,σ) block's D, C, ε, f, the fingerprint and the run
    //! summary -- to \a path (HDF5, layout in qchem.SolidState; atomic: written aside, renamed into place).
    //! Converged or not: the file records which, and a restart reports it.  THROWS on an I/O failure.
    void SaveState(const std::string& path) const;

    SolidCalculation(const SolidCalculation&)            = delete;   // owns raw resources
    SolidCalculation& operator=(const SolidCalculation&) = delete;

    //! \brief Re-run the SCF (e.g. with tighter tolerances) and hand back the ANSWERS or the REASON.
    //! Invalidates any \c Converged handed out by an earlier attempt.
    //! \c [[nodiscard]] via \c Outcome: dropping the result is a compile-time warning.
    Outcome<Converged, SCFFailure> Converge(const SCFParams& params = {});
    //! \brief Run an ANNEALED schedule: each stage in turn, seeded from the one before it.
    //!
    //! The outcome is the FINAL stage's -- which is the honest one, since the earlier stages exist only to
    //! deliver a good starting density to it.  Every stage's own summary is still printed.
    //! An empty schedule is a caller error (there is no sensible "no stages" run) and throws.
    Outcome<Converged, SCFFailure> Converge(const std::vector<SCFStage>& schedule);
    //! \brief The outcome of the attempt the CONSTRUCTOR made, so a caller that never re-converges still
    //! has to face it.  (The ctor converges on construction; without this, its result would be silently
    //! unreachable and the old \c Energy() would happily serve a failed run.)
    Outcome<Converged, SCFFailure> Result() const;
    //! Live per-iteration telemetry, REPLACING the constructor's.  Set BEFORE a Converge() call to observe it (the ctor's
    //! initial convergence has run by then: to observe it, pass the observer to the constructor).
    void OnIteration(Observer obs);

    //! Did it converge?  A cheap STATUS question, kept because reporting and tests ask it without wanting
    //! the answers.  It is NOT the gate on the results any more -- \c Converge()/\c Result() are.
    bool                   DidConverge()    const;
    size_t                 IterationCount() const;

    //! \brief What the OUTCOME DETECTORS saw (doc/OpenWork.md N1/T3-T4).  Available whether the run
    //! succeeded or failed -- on a failure it is the evidence behind \c SCFFailure::details, and on a
    //! success it is the reassurance that the order survived and the charge channel never ran away.
    const RunDiagnostics& Diagnostics() const;

    //! \name The LAST ITERATE, converged or not -- DIAGNOSTICS, deliberately named so
    //! A fixed-iteration A/B (run both arms for exactly one step and compare the arithmetic) is a real and
    //! legitimate need, and the point of N1/T1 is NOT to make unconverged numbers unreachable -- it is to
    //! make them IMPOSSIBLE TO MISTAKE FOR AN ANSWER.  These say what they are in their names; \c Energy()
    //! on the proof says what it is by existing at all.
    //!@{
    qchem::EnergyBreakdown LastIterateTerms()  const;
    double                 LastIterateCharge() const;
    //! The last iterate's DENSITY OBJECT (null before the first Converge) -- for a symmetry or order probe
    //! that deliberately runs a BOUNDED number of iterations and asks the density what it can do
    //! (\c cSpinResolved_CD for the channels), exactly as \c Converged::DensityMatrix does for an answer.
    const qchem::ChargeDensity::cDM_CD* LastIterateDensity() const;
    //!@}

    //! \brief ACBN0 (DFT+U increment 3): \f$(\bar U,\bar J)\f$ per Hubbard manifold FROM THE LAST ITERATE's
    //! orbitals -- the paper's outer-loop estimate (SCF at \f$U^{(n)}\f$, estimate, run again at
    //! \f$U^{(n+1)}=\bar U-\bar J\f$ until \f$|\Delta U|<10^{-4}\f$ eV).  Fed block by block from the wave
    //! function's occupied orbitals; the mechanism (Löwdin renormalisation, bare on-site integrals) is the
    //! Hamiltonian's.  Throws when the run carries no +U manifold.  A DIAGNOSTIC of the last iterate like its
    //! neighbours above: run it on a \c Converged run, and say which iterate it came from when you quote it.
    std::vector<qchem::Hamiltonian::ACBN0Estimate> EstimateHubbardU() const;

    //! \brief THE ACBN0 OUTER LOOP (Agapito et al.'s practice; the true variational functional is their future
    //! work too): estimate \f$(\bar U,\bar J)\f$ from the current orbitals, set every manifold's \f$U\f$ to
    //! \f$U_{\rm eff}=\bar U-\bar J\f$ on the SAME Hamiltonian, re-converge from the current density with
    //! \a params, repeat until every manifold moves less than \c tolU_eV or \c maxOuter is spent.  The result
    //! carries the whole trajectory, so a slow drift and an oscillation are visible, not hidden in a final
    //! number.  The run's +U banner value is the STARTING U; the final one is in the trajectory.
    //! ★ R4 (doc/LinearResponsePlan.md §4c, ruling D4): the SAME loop drives either estimator -- the facade picks the strategy.
    //! \c LinearResponse: each step's U is hp.x's \f$U=(\chi_0^{-1}-\chi^{-1})_{II}\f$ over the \c nq q-mesh's supercell, with
    //! +U FROZEN at that step's state (\c HubbardLinearResponse), written onto the PERTURBED manifolds only (the
    //! spectators keep theirs); the fixed point \f$U_{\rm out}(U_{\rm in})=U_{\rm in}\f$ is Timrov et al. PRB 103,
    //! 045141 (2021)'s self-consistent U.  Started at \f$U_{\rm in}=0\f$, step 1 IS \f$U_0\f$ and the last is
    //! \f$U_{\rm SC}\f$ -- ONE run fills the \c [U table] line (user, 2026-09-27).
    struct HubbardLoop
    {
        enum class Estimator {ACBN0, LinearResponse};
        size_t    maxOuter=20;
        double    tolU_eV=1e-4;
        Estimator estimator=Estimator::ACBN0;
        //! \name LinearResponse only
        //!@{
        ivec3_t             nq{1,1,1};             //!< the q-mesh (1 = q = 0 of this cell)
        KrylovParams        krylov{.tol=1e-8};
        //! The manifolds that are Hubbard SITES (perturbed, and given the estimated U); EMPTY = those carrying a U at the
        //! start, else all.  ⚠ Starting at U_in = 0 nothing carries a U, so a deck with projector spectators (orthofull)
        //! must LIST its sites, or every spectator is perturbed too -- a different U (§3d Q10) at 4x the cost.
        std::vector<size_t> perturbed;
        //!@}
    };
    struct HubbardLoopResult
    {
        std::vector<std::vector<double>> U_eV;     //!< per outer step, per manifold: the U estimated FROM that step's SCF (manifolds the
                                                   //!< estimator does not set keep the U they ran with)
        std::vector<size_t> estimated;             //!< the manifolds the loop sets (ACBN0: all; LinearResponse: the perturbed set)
        std::vector<double> Uin_eV;                //!< the U each manifold STARTED with (the U_0 end of the table runs at these)
        //! Per outer step, the state the estimate was made FROM: its band gap (eV; NaN if the basis has no empty level) and each
        //! manifold's integrated moment \f$n_\uparrow-n_\downarrow\f$ through the +U projectors (0 on an unpolarized run).
        std::vector<double>              gap_eV;
        std::vector<std::vector<double>> moment;
        std::vector<qchem::Hamiltonian::HubbardEstimate> last;   //!< the final step's estimates (either strategy; ACBN0's detail is EstimateHubbardU's)
        size_t outer=0;                            //!< outer steps run (SCF re-convergences = outer-1)
        bool   converged=false;                    //!< |ΔU| < tolU_eV on every estimated manifold
        bool   scfConverged=false;                 //!< the LAST SCF converged (an unconverged SCF poisons the estimate)
        std::string failure;                       //!< LinearResponse: the response that FAILED and ended the loop (empty = none)
    };
    //! \brief Each +U manifold's integrated spin moment \f$\mathrm{Tr}\,n^\uparrow-\mathrm{Tr}\,n^\downarrow\f$ (electrons, k-weighted) of
    //! the CURRENT orbitals, through the +U term's own projectors -- the same number the [U table] reports per outer step, and
    //! valid on ANY XC mesh (\c RunDiagnostics' order needs Becke basins).  THROWS when the run carries no Hubbard manifold.
    //! PUBLIC because a moment the panel asserts on (doc/SCFRobustnessPlan.md §4 #10: non-magnetic Si must keep m = 0 under a
    //! one-spin alpha*P) is an observable a caller judges a run by, and a uniform-mesh run has no other integrated moment.
    std::vector<double> MakeHubbardMoments() const;
    HubbardLoopResult ConvergeHubbardU(const SCFParams& params, const HubbardLoop& loop);
    HubbardLoopResult ConvergeHubbardU(const SCFParams& params) {return ConvergeHubbardU(params, HubbardLoop{});}

    //! \brief \f$\chi_0\f$ -- the INDEPENDENT-PARTICLE channel response of the last iterate over its Hubbard
    //! manifolds, on the \a Nq q-mesh (doc/LinearResponsePlan.md stage R0; Timrov et al. PRB 98, 085127's
    //! \f$\chi_0\f$).  Sum over states on the stored orbitals (every virtual included), the run's OWN occupancy
    //! rule (the final stage's configuration) and the +U term's OWN projectors; gated on the response gap
    //! against the final iteration's [F,D] commutator as the eigenvalue noise.  Reports as it computes
    //! (per q, then the real-space home-cell block, in eV\f$^{-1}\f$ -- hp.x's unit).
    //! FAILS on an incommensurate q-mesh or an inverted/unresolved coupled pair.  THROWS when the run carries
    //! no +U manifold (list the channels at U=0 to probe a run without +U) or reduced its k-mesh (D5: run
    //! without \c imposeSymmetry).  A diagnostic of the last iterate, like \c EstimateHubbardU.
    Outcome<Response::ChannelResponse,Response::ResponseFailure> IndependentResponse(ivec3_t Nq) const;
    //! \brief THE GRID AUDIT (doc/GridEpsilonPlan.md §3; C2 M2c/M3): every raster term's DELIVERED error at the converged density.
    //! Rebuilds the GPW basis and the Hamiltonian on REFERENCE rasters \a referenceMargin beyond what each of the run's rasters delivers
    //! (x_ref = x_run + margin; the default 8 puts the reference's own tail ~1e-3 of the run's: (x+8)^3 e^{-x-8} / x^3 e^{-x}), carries
    //! the converged D onto it (\c SnapshotState / \c RestoreDensity: the orbital functions are the same Gaussians), and compares term by
    //! term: max|ΔH_ij| over every block, and ΔE.  A FIXED-D evaluator: no SCF step, the run's result is untouched (its own Hamiltonian
    //! is re-refreshed for the restored D; every later caller refreshes for its own).  Beside each row it states the run's own grid-audit
    //! verdict (grids.auditDensity / auditVlocLong) AND the reference's, so the predicted alias of whichever side is the coarse one can be
    //! read against the delivered difference.  Reports to the
    //! console and to the run report (section `gridAudit`).  THROWS on a plane-wave run (no raster) or a margin <= 0.
    //! label/entries: \c TermEnergy's.  \c bound / \c eta: the term's RIGOROUS alias bound at the COARSE side's x (GPWTolerances.C
    //! \c AliasTail) and the measured eta = delivered / bound -- on max|dH_ij| for V_loc long (what forces see), on |dE|/|E| for Hartree;
    //! \c etaOn names which ("H", "E"), or "measured" for the XC raster (no rigorous bound: v_xc[rho] is nonlinear), or "" (analytic).
    //! WHERE a term's max|ΔH_ij| sits: its (spin, block), its two functions as "atom, shell, l", and D_ij there (does the density even see it?).
    struct GridAuditElement { std::string spin, functionI, functionJ; size_t block=0, i=0, j=0; double absD=0.0; };
    //! One term's row.  Beside max|ΔH_ij|, never instead of it, and never a screen (doc/Pins.md pins 19/20; OpenWork ε row step 0b):
    //! \c trDdH = the BZ-weighted Re Tr(D ΔH) of each (spin, block) channel, in canonical block order, at the converged D.  SIGNED, so it
    //! cancels and never carries a verdict alone; and not ΔE term by term (Hartree and XC are not linear in D).
    //! \c sNormDH = the largest channel's S-metric norm ‖V†ΔHV‖_F in the run's orthonormal frame (= sqrt(Tr(S⁻¹ΔH S⁻¹ΔH)) when VV† = S⁻¹, the
    //! dual of \c tDensityMatrixDistance's), which bounds the trace: |Tr(D ΔH)| ≤ sqrt(Σf²)·sNormDH.  Report only: no gate reads either yet.
    struct GridAuditTerm
    {
        std::string label, entries; double maxAbsΔH=0.0, ΔE=0.0, E=0.0, x=0.0, bound=0.0, eta=0.0; std::string etaOn;
        std::vector<double> trDdH; double sNormDH=0.0; GridAuditElement maxElement;
    };
    struct GridAuditVerdict { std::string audit, status; double alias=0.0; };   //!< one grids.audit* entry (audit = its key)
    struct GridAuditResult
    {
        std::string reference;                 //!< the reference: its rasters' x, or the egg-box shift
        std::vector<GridAuditVerdict> runVerdicts, referenceVerdicts;   //!< each side's own grid-audit verdicts
        std::vector<GridAuditTerm> terms;      //!< every term, in the Hamiltonian's order (analytic terms read 0)
        double maxAbsΔH=0.0, ΔE=0.0;           //!< the whole Hamiltonian's (ΔE = E_run - E_reference)
    };
    static constexpr double kReferenceMargin=8.0;   //!< EMPIRICAL: the audit's default x_ref - x_run, chosen for the tail ratio \c GridAudit derives
    static constexpr int    kPilotRounds=4;         //!< the pilot's cap on audit + re-convergence rounds (\c CalibrateRasters)
    GridAuditResult GridAudit(double referenceMargin = kReferenceMargin) const;
    //! \brief WHAT THIS RUN'S CLIENT CONSUMES (user ruling 2026-10-31, doc/GridEpsilonPlan.md §1): the measure the one accuracy target
    //! is judged on.  A calculation built on its own answers an \c Energy; one a \c SolidSession builds for a geometry point (a
    //! surface, forces, a relaxation) answers \c MatrixElements.  The run answers it from who built it -- the user never picks a metric.
    ClientConsumes Consumes() const override;
    //! \brief The setup with THIS run's raster sizing FROZEN into it (C2 step 3): the tolerances the pilot audit settled on, marked
    //! resolved -- what a \c SolidSession steps the rest of its surface on, so every point samples on the same rasters (a re-sized raster
    //! per point would make the surface discontinuous).  The setup itself when its rasters were resolved already.
    std::shared_ptr<const SolidSetup> SizedSetup() const;
    //! \brief What the pilot COST (C2 step 3; the controller's condition: a user pays it, so it is reported and gated, never hidden): its
    //! audit rounds, and the SCF iterations of its re-convergences (on top of the run's own).  All zero on a resolved setup.
    struct RasterPilotCost { int rounds=0; size_t reconvergeIterations=0; bool settled=false; };   //!< settled: every term met its target (not the round cap)
    RasterPilotCost PilotCost() const;
    //! \brief THE EGG-BOX GATE (doc/GridEpsilonPlan.md §5 verdict 5): every atom translated rigidly by \a shift (bohr), basis and Hamiltonian
    //! rebuilt there (the molecular basis remade from the setup's frozen recipe), the SAME AO density matrix restored (a rigid
    //! translation moves every function with its atom, so D is unchanged), and every term compared -- with the same GridAuditResult.
    //! An exact method is translation-invariant: every ANALYTIC term must read max|dH_ij| = dE = 0, and the raster terms' dE / dH are
    //! the egg-box.  FIXED-D, and sufficient as a gate: at a converged (stationary) D the self-consistent egg-box differs from this only
    //! at second order.  THROWS on a plane-wave run, on an imposed group (the SALC phases follow the origin, so blocks would not match
    //! function for function), and on a setup resolved from a pre-built basis (no recipe to remake it at the shifted atoms).
    GridAuditResult EggBox(const rvec3_t& shift) const;
    //! \brief THE XC-MESH AUDIT: the run against the SAME rasters with the XC quadrature replaced by \a reference, at the same D, term by
    //! term -- the gate of an atom-centred mesh, which the egg-box cannot see (the mesh moves with the atoms).  The pilot's nested Becke
    //! ladder is this comparison against \c qcMesh::BeckeAngularRungAbove and \c BeckeRadialRungAbove; a test gates a settled mesh against a rung of its own choosing.
    GridAuditResult XCMeshAudit(const qcMesh::MeshParams& reference) const;
    //! \brief The SELF-CONSISTENT channel response at q = 0 over the run's Hubbard manifolds
    //! (doc/LinearResponsePlan.md stage R2): \f$(1-\mathcal R_0\mathcal K)\,\delta D=\mathcal R_0V_J\f$ with the
    //! Hamiltonian's ANALYTIC kernel (Hartree + ALDA f_xc; +U frozen or at U = 0 answers zero, H4), linearised
    //! about the converged orbitals' own density.  Reports \f$\chi_0\f$, \f$\chi\f$ and
    //! \f$U_I=(\chi_0^{-1}-\chi^{-1})_{II}\f$ (eV) as it computes.  In a cell of \f$N\f$ primitive cells that
    //! U is the LR-cDFT U of that supercell (Timrov §III) -- the finite-size limit is R3's q-mesh.
    //! ★ +U IS FROZEN FOR THE WHOLE SOLVE, and restored after (doc/LinearResponsePlan.md §3d Q6): linear-response
    //! U is DEFINED with V_Hub held at its ground-state value (Timrov eq 20), so the +U term contributes no kernel.
    //! \a perturbed: the channels J to perturb (indices into the term's manifolds); EMPTY = the manifolds that
    //! carry a +U (hp.x's "Hubbard sites"), or every manifold when none does.  Every channel is MEASURED; U is read
    //! from the J x J inverse, because an inverse-response U depends on the set it spans (§3d Q10).
    //! FAILS on an inverted/unresolved pair (E1) or a Krylov solve short of \a krylov.tol.  THROWS when the run
    //! carries no +U manifold, reduced its k-mesh (D5), or has a REAL TRIM block (run with \c forceComplex: the
    //! real-block response face is not built).
    //! \note const, although it freezes and unfreezes the +U term: the Hamiltonian comes back exactly as it was.
    Outcome<Response::SelfConsistentResponse,Response::ResponseFailure> HubbardLinearResponse(
        const KrylovParams& krylov={.tol=1e-8}, std::vector<size_t> perturbed={}) const;
    //! \brief The same over the \a Nq q-MESH (doc/LinearResponsePlan.md §3d, R3 step 4): one ground state and one
    //! kernel, a Krylov solve per (q, perturbed channel), each q reported as it finishes; then \f$\chi_0\f$, \f$\chi\f$
    //! and \f$U=(\chi_0^{-1}-\chi^{-1})_{II}\f$ over the q-mesh's SUPERCELL (hp.x's definition) -- the U of an
    //! ISOLATED perturbation in that supercell, which converges with \a Nq.  \a Nq = 1 is the call above.
    //! FAILS as the q = 0 call does, and on a q-mesh that does not divide the k-mesh.
    //! \a res: how finely each manifold is split into channels (OpenWork §2 "Response: REPORT orbital-resolved channels"):
    //! \a perturbed still names MANIFOLDS, and every channel of a named manifold is perturbed; the result's \c group maps each
    //! channel back to its manifold, and the U lines add the 1x1 inversion (Macke 2024's pointwise) to hp.x's full one.
    Outcome<Response::QMeshResponse,Response::ResponseFailure> HubbardLinearResponse(
        ivec3_t Nq, const KrylovParams& krylov={.tol=1e-8}, std::vector<size_t> perturbed={},
        Response::HubbardResolution res=Response::HubbardResolution::Manifold,
        //! \a qStars (doc/KAxisPlan.md P3, Q9): solve ONE q per star of the ground state's group and generate the rest; asserts
        //! the state has that group, so one member is solved directly as an audit (a disagreement FAILS).  Manifold only.
        bool qStars=false,
        //! \a littleGroupK (doc/KAxisPlan.md P5, D5-a): solve each (q, perturbed channel) on one k-block per orbit of its little
        //! group -- the ops of the ground state's group that keep q and the perturbation, each MEASURED on the orbitals first
        //! (§5 ruling 1: a state that broke an op loses it; none kept = the full mesh).  Exact: χ equals the full-mesh χ.
        bool littleGroupK=false) const;
    //! \brief The INTER-SITE \f$V_{IJ}\f$ of a q-mesh response (doc/LinearResponsePlan.md §5 R4 row): every home-cell row of
    //! \c QMeshResponse::Interaction off its diagonal, with the PAIR's geometry -- what a reader of a V needs is "between which two
    //! atoms, how far apart".  Sorted by distance.
    //! ⚠ A FINITE q-MESH ALIASES V: each entry is the SUM over every image of \c to that differs by a q-mesh SUPERCELL translation;
    //! \c R and \c distance name its NEAREST member.  So an entry is one bond only when the supercell is large enough to separate the
    //! neighbours -- measured on Si 2x1x1 (2026-10-05): three of site 0's four equivalent nearest bonds fold onto ONE entry (1.28 eV)
    //! and the fourth is alone (-0.016 eV), both at d = 4.44 bohr.  Read V on a mesh whose supercell separates the shell you care about.
    struct InterSiteV
    {
        size_t  from=0, to=0;   //!< channel indices of the response's labels (= manifold indices when one channel per manifold; from = home cell)
        ivec3_t R;              //!< the image of \c to, in primitive lattice vectors
        double  distance=0;     //!< bohr
        double  V=0;            //!< Hartree
    };
    std::vector<InterSiteV> InterSiteHubbardV(const Response::QMeshResponse&) const;

    //! \brief The FINITE-DIFFERENCE cross-check of \c HubbardLinearResponse -- LR-cDFT's own definition
    //! (Cococcioni & de Gironcoli 2005; Timrov §III): perturb manifold \a J with the static \f$\pm\alpha\hat P_J\f$
    //! (QE's \c Hubbard_alpha), re-converge each with +U FROZEN (as the linear response holds it), and difference
    //! every channel's occupation: \f$\chi_{IJ}=[n_I(+\alpha)-n_I(-\alpha)]/2\alpha\f$.  No kernel, no solver: the
    //! independent route to suspect a DFPT number with (user, 2026-09-29).  At q = 0 of THIS cell only -- a
    //! q-mesh number is checked in the matching supercell.
    struct FiniteDifferenceChi
    {
        std::vector<std::string> labels;   //!< every channel (the measured set)
        size_t perturbed=0;                //!< J
        double alpha=0;                    //!< the step (Hartree)
        rvec_t nPlus, nMinus;              //!< the channel occupations at +alpha and -alpha
        rvec_t chi;                        //!< \f$\chi_{IJ}\f$ per I, 1/Hartree
        bool   restored=false;             //!< the closing unperturbed, unfrozen SCF converged
    };
    //! Three SCFs (+alpha, -alpha, then the unperturbed state RESTORED, unfrozen), each a fresh stage with
    //! \a params.  FAILS with the first perturbed SCF that does not converge (the state is still restored).
    //! THROWS when the run carries no +U manifold or \a J is out of range.
    //! \a reseed: false (default) starts each SCF from the CURRENT state; true starts each from the run's SEED
    //! STRATEGY (R2 gate (b)'s route) and restores from the ground-state density kept aside.  ⚠ A restart sees the
    //! MIXING undressed (measured 2026-09-29, Si Γ): its error vectors are one mode, DIIS keeps 2 of them, so a
    //! relax of 1.0 diverged where it converges from the seed; relax 0.2 restarts cleanly.  Tune \a params, not
    //! this flag.  ⚠ The re-seed is the iterator's PLAIN seed path (\c opts.seed only -- no IonicSAD site
    //! targets), so a MAGNETIC run keeps reseed=false.
    Outcome<FiniteDifferenceChi,SCFFailure> HubbardFiniteDifferenceChi(size_t J, double alpha, const SCFParams& params,
                                                                        bool reseed=false);

    // ⛔ Energy() / EnergyTerms() / TotalCharge() / Density() DELIBERATELY DO NOT LIVE HERE any more
    // (doc/OpenWork.md N1/T1).  They are on Converged, reachable only through Converge()/Result(), because
    // every one of them used to return a plausible number from a run that had failed.

    //! \brief WHICH XC quadrature this run actually used, after \c Auto was resolved.
    //!
    //! Exposed because the resolution is a real decision made here on the run's behalf, and a caller that
    //! wants to reproduce, pin, or report the run needs to see what was chosen rather than re-deriving it.
    const qcMesh::MeshParams& ResolvedXCMesh() const;

    //! \brief The GPW (Bloch) basis this run built, for a caller that must ASK IT SOMETHING -- e.g. build a
    //! fit basis to take a Fourier component of the converged density.
    //!
    //! Exposed for the same reason \c ResolvedXCMesh is: the facade MADE this object on the run's behalf,
    //! and a campaign diagnostic that needs to interrogate it should not have to rebuild a second one and
    //! hope the two agree.  It is a const view -- the calculation keeps ownership.
    const BasisSet::Complex_BS& Basis() const;

private:
    //! \name THE FD RESPONSE ORACLE'S DOOR (src/Response/tests, ruling D6) -- the converged run's Hamiltonian,
    //! wave function and composite basis, which the facade deliberately does not expose.  Friend only.
    //!@{
    friend class ::ResponseFacadeTests;
    qchem::Hamiltonian::cHamiltonian&       ResponseHamiltonian () const;
    const WaveFunction::cWaveFunction&      ResponseWaveFunction() const;
    //!@}
    //! \name THE BUILD, in the two halves Restart needs a gap between (the fingerprint is judged against
    //! the BUILT basis, before anything expensive that depends on the seed)
    //!@{
    struct BuildOnly {};
    explicit SolidCalculation(BuildOnly);
    //! THE GEOMETRY POINT over \a setup (doc/CalcNetworkReview.md §2): the Bloch basis at \a lat's geometry, its vetting report,
    //! and the electron configuration over its k-blocks.  Every DECISION was made by \c SolidSetup::Resolve; this only builds.
    //! \a mol: the molecular basis AT \a lat's geometry (the setup's own at the reference; a geometry step's, made from its recipe).
    void BuildPoint(std::shared_ptr<const SolidSetup> setup, const Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol);
    //! The +U manifolds' groups and radials, the Hamiltonian, the accelerator, the seed (\a seed, or
    //! \c opts.seed when null), the iterator and the run banner.  \a exactResume: \a seed is a restored state with no
    //! WARM difference (SeedRole::ExactResume -- it becomes the first working density, not just the first Fock's).
    //! \a accel: the first stage's accelerator, as asked.  \a restored: \a seed is a SAVED state of this same Hamiltonian (Restart);
    //! false for a carried one (a geometry step: the Hamiltonian moved with the atoms) -- the Auto accelerator policy reads it.
    void BuildRun(std::unique_ptr<qchem::ChargeDensity::cChargeDensity> seed, const SCFParams& params,
                  SCFAccelerators::Type accel, bool restored, bool exactResume);
    //! \name THE GEOMETRY POINT a SolidSession steps to (step 8)
    //!@{
    friend class SolidSession;
    //! \a lat: the setup's cell with the atoms at new positions; \a mol: the molecular basis made there from the setup's frozen
    //! recipe.  \a carry null: the setup's seed strategy and the whole \a schedule (the first point); else seeded from \a carry and
    //! converged on the schedule's FINAL stage only (the earlier stages exist to deliver a density the carry already is).
    //! \a lowdin: move \a carry onto this point's basis through the Löwdin frame (\c LowdinCarry) before it seeds the run.
    //! \a onIteration: the session's observer, attached before the point's first iteration.
    static std::unique_ptr<SolidCalculation> AtGeometry(std::shared_ptr<const SolidSetup>, const Lattice_3D& lat,
                                                        std::shared_ptr<const BasisSet::Real_BS> mol, const SavedState* carry,
                                                        const std::vector<SCFStage>& schedule, bool lowdin, Observer onIteration);
    SavedState Carry() const;   //!< the last iterate's per-block D, for the next point
    //!@}
    //! Point + run + the first Converge (the constructors' common body), and the rest of an annealed schedule.
    void Start(std::shared_ptr<const SolidSetup> setup, const SCFParams& params, SCFAccelerators::Type accel);
    //! Point + run, no SCF: \c Start's build half (and \c Prepare's whole body).
    void BuildFirstStage(std::shared_ptr<const SolidSetup> setup, const SCFParams& params, SCFAccelerators::Type accel);
    friend class ::DeckTests;   //!< tests/Deck.C: an inline deck and its named twin build the same fingerprint
    void RunRestOf(const std::vector<SCFStage>& schedule);
    //! This run's saved-state fingerprint, with \a kT the recipe's temperature.
    StateFingerprint Fingerprint(double kT) const;
    //!@}
    //! Stand up one stage's Hamiltonian + accelerator + iterator, seeded from \a carried when given.  \a params are the
    //! stage's (an \c Auto accelerator is resolved against them).
    void BuildStage(const std::string& why, SCFAccelerators::Type, std::unique_ptr<qchem::ChargeDensity::cDM_CD> carried, const SCFParams& params);
    //! Build this stage's accelerator: resolve \a asked (\c Auto, see \c SolidCalcOptions::accelerator) for \a params and
    //! a \a restored (Restart: same Hamiltonian) or seeded / carried start, and record both for the banner and for the stages
    //! that re-resolve it.
    void MakeAccelerator(SCFAccelerators::Type asked, const SCFParams& params, bool restored);
    //! \c Converge(params) after the stage's own adjustment (a restored Auto stage's linear step capped at 0.3).
    Outcome<Converged, SCFFailure> ConvergeStage(const SCFParams& params);
    //! The one place a Converged/SCFFailure is minted, so \c Converge() and \c Result() cannot drift
    //! apart in what they call a success.
    Outcome<Converged, SCFFailure> Outcome_() const;
    //! Install the composed order probe + observer on the current iterator (see the .C for why the
    //! facade takes both hooks and hands the caller's back through them).
    void AttachProbes();
    //! The ONE comparison under both audits: the run against a reference built on \a raster at the atoms shifted by \a shift (zero
    //! for the raster audit), at the same D, term by term.  \a reference names the reference in the report.  \a xcMesh: the reference's XC
    //! quadrature (null = the run's own; the pilot's nested Becke audit passes the next ladder rung).  \a ballScale: the reference's density
    //! G-ball cutoff over the run's (1 = the run's own; the pilot's nested ball audit passes 2 -- \c cutoffFactor when the ball is auto,
    //! \c densityEcut when it is stated).  \a xcΔH: when given, also hands back the XC term(s)' H_run - H_reference, one matrix per (channel,
    //! Exc term) in the fold's order -- the Becke climb's cache (\c CalibrateRasters), which compares two REFERENCE meshes through the run's.
    GridAuditResult CompareAtFixedD(const BasisSet::Gaussian::GPWTolerances& raster, const rvec3_t& shift, const std::string& reference,
                                    const qcMesh::MeshParams* xcMesh=nullptr, double ballScale=1.0,
                                    std::vector<hmat_t<dcmplx>>* xcΔH=nullptr) const;
    //! The reference tolerances of a grid audit: every raster demand at what the run's raster DELIVERS plus \a margin.
    BasisSet::Gaussian::GPWTolerances ReferenceRasters(double margin) const;
    //! \brief THE PILOT (C2 step 3, doc/GridEpsilonPlan.md §1): on a setup whose rasters are not resolved, after the FINAL stage converged
    //! with \a params, audit every raster term against references \c kReferenceMargin finer, fit each over-target term's measured tail
    //! (\c GridAccuracy::MeasuredDemand) at the target for what the client consumes (\c Consumes), and -- when any term demands more --
    //! rebuild the point on the raised rasters and re-converge from this D (a restored stage: the Ladder on its GDM rung at kT = 0).
    //! Repeats until every term meets its target (at most \c kPilotRounds); idempotent once settled.  No-op on a resolved setup.
    //! A run whose client consumes H on a Becke XC mesh also audits that mesh against the next rung on each axis of the nested ladder
    //! (\c qcMesh::BeckeAngularRungAbove, \c BeckeRadialRungAbove).  When the mesh was the facade's choice (\c xcMesh Auto) it CLIMBS at the
    //! audited D (\c qcMesh::ClimbBeckeLadder): round 0 walks the ladder to settled with no re-convergence between rungs, and every later round
    //! is the confirming audit (one raise at most).  A STATED Becke mesh is audited and reported, never raised (the grid is a user decision).
    void CalibrateRasters(const SCFParams& params);
    std::unique_ptr<Imp> itsImp;
};

//! \brief THE SOLID POTENTIAL ENERGY SURFACE (doc/CalcNetworkReview.md §3.1, §7; step 8): one fixed \c SolidSetup, sampled at
//! positions the caller chooses.  Each \c ConvergeAt builds a new geometry point -- the setup's cell with the atoms moved, its
//! molecular basis remade from the setup's frozen recipe (the trim is NOT re-vetted) -- seeded from the previous point's density
//! carried onto it.
//!
//! WHAT A POINT REFUSES (each a \c PESFailure, never a crash): positions that break the setup's IMPOSED group (ruling Q4); a
//! basis that, at these positions, would lose functions in the ortho step the trim was decided to avoid (ruling Q3: the frozen
//! basis left its conditioning window -- checked when the setup states an \c orthoTol); a moved cell (the basis and grids
//! are built for one); forces or stress (no term has its gradient yet: step 9).  +U is held at the setup's values along the surface (ruling Q5).
class SolidSession : public virtual PotentialEnergySurface
{
public:
    //! HOW the previous point's density reaches the next point's basis (the same recipe, its centres moved).  \c AO hands D over
    //! verbatim, as AO coefficients; \c Lowdin holds the orbitals' Löwdin coefficients \f$S^{1/2}C\f$ fixed (\c LowdinCarry), so the
    //! carried D stays idempotent with the right electron count in the new metric.
    enum class Carry { AO, Lowdin };
    //! \a schedule: the SCF recipe -- the first point runs all of it from the setup's seed, every later point its FINAL stage
    //! from the carried density, moved onto the new basis as \a carry says.  THROWS on a free Gaussian setup resolved with
    //! \c SolidCalcOptions::xcMeshSymmetric (the XC mesh would follow the group detected at each geometry and switch along a path).
    //! \a onIteration observes every point's SCF (each \c ConvergeAt), from its first iteration -- as \c MoleculeSession.
    SolidSession(std::shared_ptr<const SolidSetup> setup, std::vector<SCFStage> schedule, SolidCalculation::Observer onIteration = nullptr,
                 Carry carry=Carry::AO);
    ~SolidSession();
    PESGeometry Reference() const override;
    //! THROWS when \a cell is empty, when \a R has not one position per atom, or when the atoms MOVED on a setup resolved from a
    //! pre-built basis (no recipe to remake it with -- resolve the setup from a \c GaussianBasisRecipe).  A cell that differs
    //! from the setup's is \c PESFailure::CellChange (the basis and grids are built for one cell).
    Outcome<PESPoint, PESFailure> ConvergeAt(const std::optional<rmat3d_t>& cell, const std::vector<rvec3_t>& R,
                                             const PESRequest& want = {}) override;
    //! The current point (the last \c ConvergeAt's), for a caller that wants its rich answers; null before the first.  Replaced
    //! by the next \c ConvergeAt.
    const SolidCalculation* Current() const;
private:
    struct Point;
    std::shared_ptr<const SolidSetup> itsSetup;
    std::vector<SCFStage>             itsSchedule;
    SolidCalculation::Observer        itsObserver;
    Carry                             itsCarry;
    std::unique_ptr<Point>            itsPoint;
};

} //export namespace qchem
