// File: Calculation/Calculation.C  The abstract base of the three SCF facades (atom, molecule, solid).
//
// USER RULING (2026-11-02, V-FACADE-PARITY): the abstract base takes the name \c Calculation, the concrete molecular facade became
// \c MoleculeCalculation.  The base is where the TEMPLATE METHOD lives -- the high-level policy setup -> Prepare -> Converge ->
// Restart as a sequence of virtual calls, each step pure virtual or defaulted here and overridden only where the GEOMETRY differs --
// and where the shared state contract goes.  The steps move up one at a time, each flipping cells of the facade matrix
// (UTCalculation Facades/FacadeMatrix.Cell/<Kind>_<Capability>), so nothing lands here that two facades do not already share.
//
// THE STATE CONTRACT (V-FACADE-STATE, executed 2026-11-09): \c SaveState writes what every kind shares -- the wave-function blocks
// (D, C, ε, f, S per Irrep), the run summary, the structure / electron / block identity, the framing (\c qchem.SCFState) -- and each
// facade's \c RecordState override adds its own fingerprint entries.  A restart is each facade's static factory (it must build
// from its own setup), judged here (\c JudgeRestart): read, refuse, compare, carry across moved atoms.
module;
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>
export module qchem.Calculation;
export import qchem.Accuracy;   // ClientConsumes (Consumes' answer)
export import qchem.Energy;     // EnergyBreakdown (ReportResult's argument)
export import qchem.SCFState;   // the shared saved-state half (RestartRefusal: every Restart's failure)
import qchem.Structure;
import qchem.Types;             // rvec3_t (a level's k)
import qchem.Orbitals;          // Orbital (the level walk's handle)
import qchem.Symmetry.Spin;     // SpinGroup, Spin
export import qchem.SCFFailure; // SCFFailure: why an attempt has no answer (Result's failure half)
export import qchem.Outcome;    // Outcome<Converged, SCFFailure>
import qchem.ScalarFunction;    // ScalarFunction<double> (the spin density's face)
import qchem.ChargeDensity;     // tDM_CD<T> (the converged density)
export import qchem.SCFParams;  // SCFParams (Converge's edit argument)
import qchem.RunOptions;        // RunOptions (the setup's resolved options: Converge's settings, the checkpoint path)

export namespace qchem
{

//! \brief The abstract base of \c AtomCalculation, \c MoleculeCalculation and \c SolidCalculation (the template-method home).
//!
//! Data-free, derived from virtually (the project's interface convention).  A facade never opens a report run: it writes into the
//! run its CALLER opened (user ruling D, 2026-11-02).
class Calculation
{
public:
    virtual ~Calculation() = default;
    //! \brief What this run's client consumes of the Hamiltonian -- the measure the ONE accuracy target eps is a tolerance OF (pin 31;
    //! V-FACADE-ACCURACY, user ruling 2026-11-03: every facade answers it through this one rule).  An energy by default; a facade
    //! whose runs can feed a geometry or a surface (\c SolidCalculation under a SolidSession) answers MatrixElements there.
    virtual ClientConsumes Consumes() const {return ClientConsumes::Energy;}

    //! \name THE SAVED STATE (V-FACADE-STATE)
    //!@{
    //! \brief SAVE this run's state to \a path (HDF5, "qchem-state 2.0", see \c qchem.SCFState): converged or not -- a restart says
    //! which.  Written to a temporary and renamed, so a crash never destroys the previous good state.  THROWS before the first
    //! Converge (no wave function) and on an I/O failure.
    void SaveState(const std::string& path) const;
    //! \brief WHAT A RESTART OF THIS RUN FROM \a path WOULD DO, without doing it (a GUI's "can I start from that file?"): the
    //! REFUSAL (unreadable, another kind, a refuse-class difference -- each named), or the warm differences, one line each with
    //! what it means (EMPTY = an exact resume).  Judged against this run as it stands (its last stage's kT).
    Outcome<std::vector<std::string>, RestartRefusal> JudgeState(const std::string& path) const;
    //! What this run's restart did (exact or warm, every difference, the displacement, the electrons the restored density held on
    //! this basis); EMPTY for a run that did not start from a saved state.
    //! A REFERENCE (into the facade): a range-for over a by-value temporary's ->warm dangled (lane A, 2026-11-09).
    virtual const std::optional<RestartSummary>& RestartedFrom() const = 0;
    //!@}

    //! \name THE LEVELS (structure-neutral; OpenWork §2 "Levels on the Calculation face", lane G 2026-11-12)
    //!@{
    //! One eigenlevel of the converged state: what a Levels table and an MO diagram list, for ANY structure.
    struct Level
    {
        double      e = 0.0;           //!< eigen-energy (hartree)
        double      occupation = 0.0;  //!< electrons in it, PHYSICAL (a solid's BZ weight is \c weight, never folded in), a degenerate
                                       //!< level's partners included (an atom's p shell holds up to 6); 0 for a virtual
        Spin        spin = Spin::None; //!< Up / Down for a polarized run, None for an unpolarized one
        rvec3_t     k{0,0,0};          //!< the Bloch vector it belongs to, FRACTIONAL (reciprocal-lattice units); Γ for a molecule or an atom
        double      weight = 1.0;      //!< its k-point's BZ weight (1 for a molecule or an atom): Σ weight·occupation = N
        std::string label;             //!< the orbital's QN label (its irrep under symmetry)
    };
    //! \brief EVERY level of the converged state, virtuals included, ascending in energy (spin channels and k-points interleaved;
    //! ties keep the saved state's block order: spin outer, basis blocks inner).  Built by THIS base from the facade's state blocks
    //! (\c StateBlocks), so every kind answers it the same way.  An owned list (a \c Make: it walks the wave function on each call).
    //! THROWS before the first Converge (no wave function).
    std::vector<Level> MakeLevels() const;
    //!@}

protected:
    //! A level and the orbital it is, in \c MakeLevels' order: what a facade that also SAMPLES orbitals builds its own list from, so
    //! its indices are the neutral list's (the molecule's \c GetLevels).  The orbital is BORROWED from the wave function.
    struct LevelOrbital
    {
        Level                     level;
        const Orbitals::Orbital*  orbital = nullptr;
    };
    std::vector<LevelOrbital> WalkLevels() const;

    //! \brief THE RESULT LINE: a run reports its own outcome, contemporaneously, as one \c "[label] CONVERGED iters= charge= Eelec=
    //! Etot= (terms)" console line and a \c result section of the run report.  Each facade calls it where its SCF ends, from what
    //! every facade has; a solid's site moments follow as its own lines.  Rendered at Verbose only, since the line already says it.
    //! A grid Hartree route adds what its mesh lost of the charge, and that charge's energy against \a energyTarget when one is
    //! stated (user ruling 2026-11-06: reported against eps, not rescaled away).  \a iterations is the RUN's total; \a pilotIterations
    //! is the share a grid pilot spent re-converging on the grids it raised, printed beside it and filed as \c pilotIterations.
    static void ReportResult(const std::string& label, bool converged, size_t iterations, double charge, const EnergyBreakdown& E,
                             double energyTarget=0.0, size_t pilotIterations=0);

    //! \name THE STATE HOOKS each facade answers (the template method's steps)
    //!@{
    //! What every kind's run is, beyond its structure and its blocks.
    struct StateIdentity
    {
        int         Nelec = 0;
        SpinGroup   spin  = SpinGroup::UnPolarized;
        int         multiplicity = 0;          //!< as stated (0 = the minimal spin)
        std::string functional;                //!< the model / XC functional, named (WARM: an LDA -> PBE restart is a continuation)
        double      kT = 0.0;                  //!< the last stage's smearing (WARM: the occupations are re-filled)
    };
    //! The record a save writes and a restart compares.  \c withOrbitals: the blocks carry their orbitals and the summary is filled
    //! (a save); a restart's record needs neither (the fingerprint and the block walk only).
    struct StateRecord
    {
        bool                    withOrbitals = false;
        StateFingerprint        fp;
        StateSummary            summary;
        std::vector<StateBlock> blocks;
    };
    //! \brief THE BASE fills what every kind shares: the blocks (\c StateBlocks), the summary (a save), and the structure, electron
    //! and block keys of the fingerprint (\c AddStructure, \c AddElectrons, \c AddBlocks, scf.kT).  A facade OVERRIDES it, calls
    //! it, and adds its own keys (its basis, its meshes, its model's knobs).  \a kT, when >= 0, replaces the identity's: the stage
    //! a restart is about to run.
    virtual void RecordState(StateRecord&, double kT) const;
    virtual std::string             StateKind() const = 0;                         //!< "atom" | "molecule" | "solid"
    virtual const Structure&        StateStructure() const = 0;                    //!< the run's atoms (positions: the geometry)
    virtual StateIdentity           MakeStateIdentity() const = 0;
    //! The run's blocks in the canonical order (spin irreps outer, basis blocks inner), with their orbitals when \a withOrbitals.
    virtual std::vector<StateBlock> StateBlocks(bool withOrbitals) const = 0;
    virtual StateSummary            MakeStateSummary() const = 0;                  //!< THROWS before the first Converge
    //! The kind's schema history, applied to a saved fingerprint before it is compared (none by default; the solid's has several).
    virtual FingerprintMigration    StateMigration() const {return nullptr;}
    //!@}

    //! \brief A RESTART'S JUDGEMENT, shared by every facade's Restart once its basis is built: refuse another kind's state, compare
    //! \a saved's fingerprint with this run's at \a kT (after \c StateMigration), carry every block across MOVED atoms as AO coefficients
    //! (\c CarryToMovedAtoms: the measured best first iterate; Tr(D S') then N to first order, a diagnostic), and print the EXACT
    //! RESUME / WARM START lines under \a tag.  On success \a saved holds the densities to restore, and the summary says what the
    //! restart is (the facade keeps it for \c RestartedFrom).  A REFUSAL is printed too.
    Outcome<RestartSummary, RestartRefusal> JudgeRestart(SavedState& saved, double kT, const std::string& tag) const;

    //! \brief A VIEW'S JUDGEMENT (\c StateView, lane G 2026-11-13): may \a saved be DRAWN on this run's basis, as it stands?  The
    //! restart's refusals (another kind, a refuse-class difference) and ONE MORE: ANY moved atom (over \c kViewMoveTol), with no
    //! carry -- a view draws the state where it was computed, and a saved D on edited positions is a silent wrong picture (the
    //! carry exists to SEED an SCF, which a view never runs; controller ruling 2026-11-13).  The refusal names the largest move,
    //! its atom and the remedy.  The warm differences are accepted and listed (\c RestartSummary::warm): drawing
    //! \f$\rho=\sum D\chi\chi\f$ needs none of them.  Needs only the basis and the setup (\c RecordState without orbitals), so
    //! a facade answers it before any Hamiltonian exists.  Writes nothing; the caller reports.
    Outcome<RestartSummary, RestartRefusal> JudgeView(const SavedState& saved) const;
    //! The largest atom move (bohr) a view accepts: the positions must be the SAME numbers, so this is a round-off guard, not a
    //! physics tolerance.  A deck -> json -> deck round trip reproduces them bit-exactly (FacadeMatrix ViewMatchesTheLiveState).
    static constexpr double kViewMoveTol = 1e-10;

    //! \brief THE CHECKPOINT (CK-1), every kind's (V-FACADE-PARITY P3): \c SaveState to \a path, converged or not -- the unattended
    //! run's safety net.  EMPTY \a path = no checkpoint.  A failed WRITE WARNS under \a label and the run goes on: losing a finished
    //! SCF to a full disk would make the net the hazard (\c SaveState called directly throws as usual).
    void Checkpoint(const std::string& path, const std::string& label) const;
};

//! \brief The ZERO FIELD, analytic: \f$m(r)=0\f$, the spin density of an unpolarized run (user ruling R2, 2026-11-17: "Return zero is
//! better than throwing ... as long as we don't accidentally call it 10^6 times for nothing").  Free to sample, never a raster of zeros;
//! a hot loop that only needs to know WHETHER there is spin asks the run's \c SpinGroup, not this field.
struct ZeroField : public virtual ScalarFunction<double>
{
    double         operator()(const rvec3_t&) const override {return 0.0;}
    vec3_t<double> Gradient  (const rvec3_t&) const override {return vec3_t<double>(0.0, 0.0, 0.0);}
};

//! \brief THE ONE PROOF TYPE (V-FACADE-PARITY P2, user ruling R1 2026-11-17: "Non templated base is a well known pattern"): the
//! layer of \c Calculation that knows the run's scalar type \a T (\c double for the atom and the molecule, \c dcmplx for the solid's
//! Bloch blocks), and so can hand back the converged density AS a density, \c tDM_CD<T>, with no cast down from \c ScalarFunction.
//!
//! WHY A SEPARATE PROOF TYPE rather than \c Fallible<double> on each accessor (user, 2026-08-25).  "Did it converge?" is ONE fact
//! governing a whole FAMILY of queries; wrapping each accessor repeats it and asks the caller to check N times.  Worse, a Fallible-style
//! accessor ASSERTS on misuse, and an assert is compiled out under NDEBUG -- exactly how the imposed-mesh site-blocks defect stayed
//! invisible in Release.  Here there is no check to forget and nothing to compile out: \b the \b answers \b do \b not \b exist
//! \b without \b the \b proof.  (CLAUDE.md: give capabilities only to types that have them.)
//!
//! Each kind answers the protected \c Proof* hooks and \c MakeFailure; \c Converged and \c Result are written ONCE, here.
template <class T> class tCalculation : public virtual Calculation
{
public:
    //! \brief PROOF that the SCF converged, and the only place the answers live.
    //! \warning NON-OWNING, and invalidated by the next \c Converge(): take it fresh from each attempt; do not store one across a re-run.
    class Converged
    {
    public:
        double          Energy()         const {return itsCalc->ProofEnergyTerms().GetTotalEnergy();}   //!< total energy (hartree)
        EnergyBreakdown EnergyTerms()    const {return itsCalc->ProofEnergyTerms();}
        double          TotalCharge()    const {return itsCalc->ProofCharge();}       //!< the electron count (a charge-conservation check)
        size_t          IterationCount() const {return itsCalc->ProofIterations();}   //!< the run's SCF iterations, a pilot's included
        //! \brief \f$\rho(r)\f$ of the converged state, AS A DENSITY: sample it, or ask it what it can do (its capabilities are the
        //! density's own faces, never a cast down from \c ScalarFunction).  The same lifetime as the proof.
        const ChargeDensity::tDM_CD<T>& Density() const {return itsCalc->ProofDensity();}
        //! \brief \f$m(r)=\rho_\uparrow-\rho_\downarrow\f$; on an unpolarized run the ZERO FIELD (R2), never a null or a throw.
        //!
        //! WHY (ρ, m) AND NOT (ρ↑, ρ↓).  They are the same information -- \f$\rho_\uparrow=(\rho+m)/2\f$ -- but charge and spin are the
        //! channels that BEHAVE differently: the Hartree restoring force acts on ρ alone, which is why Kerker preconditions the charge
        //! channel and has no business on the spin one.  Exposing (ρ, m) makes the physical split the obvious one, and it is what a
        //! magnetic diagnostic wants: \f$m(r)\f$ near a site, in one call.
        const ScalarFunction<double>&   SpinDensity() const {return itsCalc->ProofSpinDensity();}
    private:
        friend class tCalculation;
        explicit Converged(const tCalculation* c) : itsCalc(c) {}
        const tCalculation* itsCalc;
    };

    //! \brief AN EDIT OF THE SETUP'S SCF SETTINGS (R10 (d), 2026-11-17): it starts from the setup's \c opts.scf and changes only
    //! the fields it names, so the setup's eps, kT and mixer survive a one-field override.  One line at the call site:
    //! \c calc.Converge([](SCFParams& p){ p.NMaxIter=60; }).  (A whole \c SCFParams in its place would reset every OTHER field to its
    //! library default -- the hazard R10 retired.)
    using SCFEdit = std::function<void(SCFParams&)>;
    //! \brief THE TEMPLATE METHOD (V-FACADE-PARITY P3): run the SETUP's SCF settings (R10: they live in \c RunOptions::scf, never in
    //! this call), edited by \a edit when one is given; checkpoint the state; answer the outcome.  Every kind's \c Converge is this one:
    //! the steps it calls are the kind's (\c RunStage) and the base's (\c Checkpoint, \c Result).  The stage's params AS RUN are what
    //! the run report's `scf` header files, so an edit is visible there.  Invalidates any \c Converged from an earlier attempt.
    Outcome<Converged, SCFFailure> Converge(const SCFEdit& edit = {})
    {
        SCFParams params=Options().scf;
        if (edit) edit(params);
        return ConvergeStage(params);
    }

    //! \brief The outcome of the LAST attempt (a Prepare'd run that never converged says "did not converge in 0 iterations"): the
    //! answers, or the reason there are none.  \c [[nodiscard]] via \c Outcome.
    Outcome<Converged, SCFFailure> Result() const
    {
        using O=Outcome<Converged, SCFFailure>;
        if (std::optional<SCFFailure> f=MakeFailure()) return O::Fail(std::move(*f));
        return O::Ok(Converged(this));
    }

protected:
    //! \brief ONE STAGE, the same steps on every kind: the kind runs it (\c RunStage), the base checkpoints it (CK-1: after every
    //! stage, converged or not, when \c saveStateTo is set) and answers it.  \c Converge runs one; a kind's own multi-stage driver (the
    //! solid's annealed schedule, its +U loop, its raster pilot) runs one per stage, so each is checkpointed alike.
    Outcome<Converged, SCFFailure> ConvergeStage(const SCFParams& params)
    {
        if (RunStage(params)) Checkpoint(Options().saveStateTo, Options().label);
        return Result();
    }
    //! \name THE STEPS EACH KIND TAKES (\c Converge's)
    //!@{
    //! The run's options as the setup resolved them: \c scf is what \c Converge runs, \c saveStateTo where it checkpoints.
    virtual const RunOptions& Options() const = 0;
    //! \brief Build, iterate and report ONE SCF stage on \a params (as run: the caller's edit applied).  FALSE when it was refused
    //! before any iteration (a resource refusal): there is no new state to checkpoint, and \c Result says why.
    virtual bool RunStage(const SCFParams& params) = 0;
    //!@}
    //! \name WHAT EACH KIND ANSWERS (read only through a \c Converged, so only after a converged attempt)
    //!@{
    //! Why the last attempt has no answer, its detectors ranked by the kind; EMPTY when it converged.
    virtual std::optional<SCFFailure>       MakeFailure()      const = 0;
    virtual EnergyBreakdown                 ProofEnergyTerms() const = 0;
    virtual double                          ProofCharge()      const = 0;
    virtual size_t                          ProofIterations()  const = 0;
    virtual const ChargeDensity::tDM_CD<T>& ProofDensity()     const = 0;
    //! m(r): the kind's own (built once per attempt, on first use) when polarized, else \c TheZeroField.
    virtual const ScalarFunction<double>&   ProofSpinDensity() const = 0;
    //!@}
    //! The one zero field every unpolarized kind hands out (stateless, so one instance serves every run).
    static const ScalarFunction<double>& TheZeroField() {static const ZeroField z; return z;}
};


} // namespace qchem
