// File: Calculation/SolidCalculation.C
//
// qchem::SolidCalculation -- the production front door for a PERIODIC (GPW / plane-wave) SCF.  The
// solid peer of qchem::Calculation (molecules) and qchem::AtomCalculation (atoms): it owns the whole
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
#include "forward.H"   // RunDiagnosticsTests -- the unit-test friend (CLAUDE.md: tests may cheat)
#include <functional>  // the order-parameter probe
#include <memory>
#include <vector>     // the anneal schedule
#include <string>    // SCFFailure::details
#include <vector>
#include <string>
#include <utility>
export module qchem.SolidCalculation;

import qchem.Lattice_3D;                      // Lattice_3D (crystal + BZ mesh)
import qchem.Structure;                       // Structure
import qchem.ScalarFunction;                  // ScalarFunction<double> (the rho(r) face)
import qchem.BasisSet;                        // Real_BS (orbital source), Complex_BS (the Bloch basis)
import qchem.BasisSet.Lattice_3D.BasisSet;    // GPWFactory, GPWParams, RasterPolicy, CellImages
import qchem.Hamiltonian.Factory;             // Pol, VxcFit, the cHamiltonian solid door
import qchem.SCFAccelerator.Factory;          // Type, SolidAcceleratorOptions, the typed solid door
import qchem.SCFIterator;                     // SolidSCFIterator, SCFParams, SCFProgress, EnergyBreakdown
import qchem.ChargeDensity;                   // cDM_CD
import qchem.ChargeDensity.Seed;              // SeedStrategy
import qchem.Mesh;                            // qcMesh::MeshParams / UnitCellKind
import qchem.LASolver;                        // qchem::Ortho
import qchem.Outcome;                         // Outcome<T,E> -- the fallible-call vocabulary (N1/T1)

export namespace qchem
{

//! \brief How to set up a periodic calculation.  Designated-initializer friendly:
//!     SolidCalculation calc(lattice, basis, {.Nelec=8, .species={{"Si",4}}});
struct SolidCalcOptions
{
    //! \name The system
    //!@{
    int Nelec = 0;                                        //!< VALENCE electrons per cell (Zion, not Z).
    //! Spin multiplicity 2S+1.  0 (default) = minimal spin: the ζ=0 collapse for even \c Nelec, a doublet
    //! for odd.  >=1 promotes to the SPIN-NATIVE two-channel run (nUp-nDown = 2S).  NB unlike the molecular
    //! facade, 1 here means the EXPLICIT two-channel singlet (nUp=nDn) -- the ζ=0 cross-check of the
    //! polarized machinery against the unpolarized anchors, not a shortcut to it.
    int multiplicity = 0;
    //! Pseudopotential species as `(element, valence)`.  Multi-species is the norm for a compound; a
    //! per-Z router PP is built so each atom gets its own GTH entry.
    std::vector<std::pair<std::string,int>> species;
    //!@}

    //! \name Grids (the efficiency lever for ionic systems)
    //!@{
    double densityEcut  = -1.0;   //!< <0 AUTO: floor at \c cutoffFactor*alpha_max (recommended).
    double cutoffFactor = 2.0;    //!< C in that floor; 2 = the density's own product exponent.
    double ladderFactor = 4.0;    //!< multigrid ladder gradation (DEPTH is automatic).
    BasisSet::Lattice_3D::RasterPolicy raster = BasisSet::Lattice_3D::RasterPolicy::BallOnly;
    BasisSet::Lattice_3D::CellImages   images = BasisSet::Lattice_3D::CellImages::Periodic;
    rvec3_t kShift = rvec3_t(0,0,0);
    //! \brief The XC real-space quadrature.  DEFAULT \c Auto = "you choose": the class costs a uniform
    //! mesh sized to this run's sharpness against the atom-centred Becke mesh and takes the cheaper, with
    //! a margin (V1.26/V2.4).  Pin \c cellKind to \c Uniform or \c Becke to overrule it and be obeyed --
    //! the grid is a user decision, and the library's job is only to warn about a dominated one.
    qcMesh::MeshParams xcMesh{.cellKind=qcMesh::UnitCellKind::Auto};
    //! WHICH fit basis represents v_xc -- ORTHOGONAL to the grid above.  \c Auto picks Delta whenever the
    //! plane-wave fit cannot serve (Becke grid, or a polarized run).
    Hamiltonian::VxcFit vxcFit = Hamiltonian::VxcFit::Auto;
    //!@}

    //! \name Convergence machinery
    //!@{
    SCFAccelerators::Type accelerator = SCFAccelerators::Type::DIIS;
    bool globalFermi    = false;  //!< metal: one μ across the k-mesh.
    //! Impose the detected space group (IBZ k-fold + per-iteration density star-average + site-adapted
    //! Becke mesh).  OPT-IN (default \c false): it roughly doubles the XC grid, and see the warning below.
    //!
    //! \note DELIBERATELY DIVERGES from the test harness's \c GpwOptions, which defaults this to \c true --
    //! while its own comment says "OPT-IN per the §3 pin".  That comment states the intent and the value
    //! contradicts it (doc/CleanupCandidates.md V1.30).  This facade follows the stated intent, because
    //! default-ON is the more dangerous way to be wrong: a symmetry imposition you did not ask for is
    //! invisible in the result, whereas a missing one merely costs time.
    //! \warning Do NOT combine with a polarized ANTIFERROMAGNETIC density until the fold is spin-aware:
    //! the op set comes from the CHEMICAL space group, which contains sublattice-exchanging operations,
    //! so star-averaging each channel under it averages the magnetic sublattices together and silently
    //! demagnetises the run (doc/CleanupCandidates.md V1.28).
    bool imposeSymmetry = false;
    qchem::ChargeDensity::SeedStrategy seed = qchem::ChargeDensity::SeedStrategy::Uniform;
    qchem::Ortho ortho    = qchem::Cholesky;
    double       orthoTol = 0.0;
    //! ANSATZ POLICY (doc/RealComplexPlan.md §6): force every Bloch block COMPLEX even where the run
    //! admits real ones (a TRIM k-point under a realness-preserving Hamiltonian -- with the default
    //! LDA stack every Γ / zone-boundary block runs REAL: real S/T/V/KB, real eigensolve, real
    //! quadrature GEMMs).  The physics is identical either way (the 3c-3 acceptance gate pins ON==OFF
    //! to machine precision); complex is the DOWNGRADE direction (§1), kept as the escape hatch for
    //! complex-instability experiments and as the gate's A/B door.
    //!
    //! NOT needed for MOM any more (R2.21, 2026-08-17).  It used to be: a real block's MOM reference had
    //! nowhere to live in a run-typed reference home, so switching \c SCFParams::UseMOM on threw mid-SCF
    //! and the documented workaround was to come here and pay for an all-complex run.  The occupation
    //! state now keys each block's reference by the BLOCK's scalar, so MOM runs REAL on a mixed mesh --
    //! measured identical to the forced-complex twin (Si (3,1,1), ΔE = 1.8e-15) at the real run's cost.
    //! No knob, and no surprise: turning MOM on no longer changes what the run is made of.
    bool forceComplex = false;

    //! \name Knobs the GPW test driver carried that this struct did not (2026-08-25)
    //! Added so a periodic recipe can be stated in ONE place instead of being assembled across a test's
    //! options struct, its run driver and its accelerator helper (user: "some parameters are set in TESTF,
    //! some in RunGpw, some in MakeGpwAccelerator, all spread out by 100's of lines").
    //!@{
    //! One Fermi level shared by BOTH spin channels (the moment becomes an OUTPUT) rather than a fixed
    //! per-channel electron count.
    bool spinsShareFermi = false;
    //! Impose the DETECTED grey group with the magnetic decoration suppressed -- the erasure control.
    bool greyImposition  = false;
    //! Pin the MOM reference to the SEED's own occupied subspace before iteration 1, so a hot stage's
    //! CHARACTER continues, not just its density.
    bool momFromSeed     = false;
    //! \brief THE MAGNETIC DECORATION: one entry per cell site, +1 / -1 / 0, naming which sublattice each
    //! site belongs to.  It is what makes an imposition SHUBNIKOV rather than merely spatial -- without it
    //! an "imposed AFM" run star-averages under the spatial group and the order is ERASED.
    //!
    //! EMPTY (the default) = DERIVE it from the seed, which reproduces the historical driver behaviour and
    //! is right for the common case.  SET IT to state a specific ordering: a cell usually admits several
    //! (AFM-I vs AFM-II vs ferrimagnetic), the derived one is only the seed's guess, and comparing orderings
    //! is the whole point of a magnetic campaign.  Ignored unless \c imposeSymmetry.
    //! \note \c greyImposition forces the empty decoration regardless -- that arm is the erasure control.
    std::vector<int> siteSpins = {};
    // NB: the test driver's `realTRIMBlocks` is NOT added here -- it is \c forceComplex with the opposite
    // polarity (both feed GPWParams::hamPreservesReal).  One decision, one field.
    //! Names this run in its own console output -- the fingerprint, the stage banners, the order trace.
    std::string label = "gpw";
    //! \brief Per-iteration telemetry, live FROM CONSTRUCTION.
    //!
    //! It belongs in the options block rather than only on \c OnIteration because the constructor
    //! CONVERGES: an observer attached afterwards has already missed stage 0 -- which, on an annealed
    //! recipe, is the stage whose trajectory the fingerprint most wants.  Put it here and the run is
    //! observed from its first iteration, with the rest of the recipe, in one place.
    qchem::SCFIterator::SolidSCFIterator::Observer onIteration = nullptr;
    //! \name The run's ORDER PARAMETER -- named, and measured every iteration
    //! A campaign watches a scalar the library cannot know about (a staggered moment, a charge disproportion,
    //! a distortion amplitude).  It belongs beside \c onIteration for the same reason: the ctor CONVERGES,
    //! so a probe attached afterwards has already missed stage 0.  Empty probe = no order column.
    //!@{
    std::string orderName;
    std::function<double(const qchem::ChargeDensity::cDM_CD&)> orderProbe;
    //!@}
    //!@}
    //!@}
};

//! \brief WHY an SCF attempt produced no usable answer (doc/OpenWork.md N1).
//!
//! Carried as a VALUE, not thrown: a non-convergent SCF is a legitimate outcome the caller must act on,
//! not a broken invariant.  \c why exists because the detectors distinguish cases that look identical from
//! the energy alone -- a run that simply ran out of iterations, one whose Hartree term ran away
//! (measured 2026-08-25: Eee 13.5 -> 29.0 Ha), and one that converged to the WRONG STATE, e.g. an imposed
//! AFM cell that relaxed to zero moment, which contradicts its own constraint.
struct SCFFailure
{
    enum class Why
    {
        NotConverged,   //!< ran to the iteration limit with the residual still above tolerance
        ChargeSlosh,    //!< the Hartree term ran away (the low-G charge mode un-damped)
        //! \brief the run did NOT keep the magnetic order it was set up to carry -- either it converged
        //! to zero moment, or the moment rose and then collapsed mid-run.  An imposed Shubnikov group is
        //! an ASSERTION about the answer, so contradicting it is a failed postcondition, not a diagnostic.
        OrderLost
    };
    Why         why        = Why::NotConverged;
    size_t      iterations = 0;
    double      residual   = 0.0;
    //! The last iterate's energy.  DIAGNOSTIC ONLY -- it is deliberately not reachable as "the energy",
    //! because handing back a plausible number from a failed run is the defect this whole type exists to
    //! prevent (four such numbers were produced in one session: -45.5, -46.3, -56.4, -38.5).
    double      lastEnergy = 0.0;
    std::string details;
};

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
    //! aufbau/shared-mu recipe with `GPW_XC_DM_SOURCE=1`: 14.42 -> 29.00 Ha).  The flat-Kerker arm is
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
class SolidCalculation
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
        //! \c FourierDensity to take a G-space component, \c cPolarized_CD to reach the channels.
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

    //! Build the whole graph and converge once with \a params.  \a mol is the MOLECULAR orbital basis the
    //! GPW basis is built over (shared: one per k-block).  \a acc tunes the chosen accelerator; its
    //! defaults are the GPW production recipe.
    SolidCalculation(const Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                     const SolidCalcOptions& opts = {},
                     const SCFParams& params = {},
                     const SCFAccelerators::SolidAcceleratorOptions& acc = {});
    //! \brief The ANNEALED constructor: build the graph and run \a schedule end to end.  The peer of the
    //! single-stage ctor above, and the shape an annealed recipe wants -- one options block, one list of
    //! stages, nothing assembled by the caller in between.
    SolidCalculation(const Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                     const SolidCalcOptions& opts,
                     const std::vector<SCFStage>& schedule,
                     const SCFAccelerators::SolidAcceleratorOptions& acc = {});
    ~SolidCalculation();

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
    //! Live per-iteration telemetry.  Set BEFORE a Converge() call to observe it (the ctor's initial
    //! convergence runs before any observer can be attached).
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
    //!@}

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
    //! Stand up one stage's Hamiltonian + accelerator + iterator, seeded from \a carried when given.
    void BuildStage(SCFAccelerators::Type, std::unique_ptr<qchem::ChargeDensity::cDM_CD> carried);
    //! The one place a Converged/SCFFailure is minted, so \c Converge() and \c Result() cannot drift
    //! apart in what they call a success.
    Outcome<Converged, SCFFailure> Outcome_() const;
    //! Install the composed order probe + observer on the current iterator (see the .C for why the
    //! facade takes both hooks and hands the caller's back through them).
    void AttachProbes();
    std::unique_ptr<Imp> itsImp;
};

} //export namespace qchem
