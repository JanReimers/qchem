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
    //! stated (user ruling 2026-11-06: reported against eps, not rescaled away).
    static void ReportResult(const std::string& label, bool converged, size_t iterations, double charge, const EnergyBreakdown& E,
                             double energyTarget=0.0);

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
};

} // namespace qchem
