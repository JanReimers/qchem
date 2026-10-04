// File: Calculation/Deck.C  The INPUT DECK: a run as one serialized, versioned JSON value (D-ENV step 6;
// doc/Records/EnvKnobInventory.md §9).
//
// THE RULES (user, 2026-10-03):
//   * The deck is the ONLY way to set a tier-1/2 value.  The environment is not a second path.
//   * The file the user wrote is never the record of what ran.  Every run writes a numbered REVISION
//     `<name>.r<NNN>.json`: the deck + every `--set` + every default filled in + a header + a provenance block.
//     Re-running that file reproduces the run; re-running it with one `--set` is the intentional change.
//   * A key beginning with `_` (`_doc`, `_comment`) is a COMMENT and is ignored -- the convention of materials.json; a record does not keep it.
//   * A reader REJECTS an unknown key (a typo in a deck must not silently run the default), and a missing key keeps
//     the field's default (so an old deck keeps working when a field is added).  Units are the library's: atomic
//     units where it matters are in the key name: energies a person quotes are eV (`U_eV`, converted at the boundary -- RAM is atomic units).
//
// This unit is the DATA LAYER: the typed option structs <-> JSON, the dotted-path `--set`, the revision claim.
module;
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>
#include <nlohmann/json.hpp>
export module qchem.Deck;

import qchem.SolidCalculation;     // SolidCalcOptions
import qchem.SCFParams;            // SCFParams
import qchem.Mesh;                 // MeshParams
import qchem.RunPolicy;            // RunPolicySpec
import qchem.BasisSet.Gaussian.Point.Factory;   // BasisSetData (the deck names the basis file)
import qchem.BasisSet.Gaussian.Point.ShellTrim;  // ShellTrim (a stated trim)
import qchem.Types;                // ivec3_t
import qchem.SCFAccelerator.Factory;   // SCFAccelerators::Type (a stage names its accelerator)
import qchem.Materials;            // Material (the structure section's name resolves to one)
import qchem.BasisSet.Gaussian.Lattice.GPWTolerances;

export namespace qchem::deck
{
using json = nlohmann::json;

//! The schema version written into every deck.  Bumped when a key is RENAMED or its meaning changes (adding a key does not).
inline constexpr int kSchemaVersion = 1;

//! \name Typed option structs <-> JSON.  \c To* writes EVERY field (a resolved deck is complete); \c From* overwrites only
//! the keys present and THROWS \c std::runtime_error on a key it does not know, naming the object path and the legal keys.
//!@{
json ToJson(const BasisSet::Gaussian::GPWTolerances&);   void FromJson(const json&, BasisSet::Gaussian::GPWTolerances&);
json ToJson(const SCFParams&);                           void FromJson(const json&, SCFParams&);
json ToJson(const qcMesh::MeshParams&);                  void FromJson(const json&, qcMesh::MeshParams&);
//! The solid options.  NOT serialized (not choices): \c onIteration (a callback).  \c Hubbard manifold \c siteOps / \c greyOps
//! (derived by the facade from the structure).  A manifold's U is \c U_eV (and \c Uirrep_eV, \c alpha_eV): eV in the file, converted to Hartree on read; the writer picks the eV value that converts back to the SAME double.
json ToJson(const RunPolicySpec&);                       void FromJson(const json&, RunPolicySpec&);
json ToJson(const SolidCalcOptions&);                    void FromJson(const json&, SolidCalcOptions&);
//!@}

//! \brief One run, as the deck states it.  The STRUCTURE is just a NAME (user, 2026-10-03): the key into `materials.json` (a
//! crystal or box -> a solid run) or `molecules.json` (a molecule -> a molecular run, not yet supported by the deck).  Everything
//! the structure files own -- lattice, atoms, spin decoration, the pseudopotential vocabulary -- is resolved FROM the name, never
//! restated in the deck.
//! \c solid.Nelec / \c solid.species left at their defaults (0 / empty) are DERIVED from the material; stated, they are honoured
//! (a charged cell, a different pseudopotential valence) -- and the resolved deck always records the values actually used.
struct RunSpec
{
    //! The basis the GPW basis is built over: WHICH data file (an exact \c BasisSetData name, e.g. "VALENCE_LOWQ_VA"), and whether
    //! it is wrapped in the spherical lattice view (contaminant-free d/f).  (Shell trims are not in the deck yet.)
    struct Basis
    {
        //! One shell removed from the basis (pin 22): element, the angular momentum dropped, and the shell's (single) exponent that names it.
        struct Trim { int Z=0; int l=0; double alpha=0.0; };
        BasisSet::Gaussian::BasisSetData data = BasisSet::Gaussian::BasisSetData::VALENCE_LOWQ_SR;
        bool spherical=false;
        std::vector<Trim> trim;     //!< a STATED trim (built once, no vet loop)
        bool vet=false;             //!< the vet-stage trim at \c solid.orthoTol: near-dependent diffuse shells removed ONCE on the full k-mesh (exclusive with \c trim)
    };
    //! Saved states (CK-1).  \c save: empty = never, "auto" = `<outDir>/states/<stem>.h5` after every stage, else a path.  \c restartFrom: a REVISION
    //! STEM (`MnO_AFM2.r003` -> `<outDir>/states/MnO_AFM2.r003.h5`, and that revision becomes this run's \c parent) or a path; runs the schedule's FINAL
    //! stage only (a restart is one stage -- the earlier stages exist to deliver the density the file already holds).
    struct State { std::string save, restartFrom; };
    //! One stage of an annealed recipe: its SCF parameters and its accelerator (\c SCFStage, by name).
    struct Stage { SCFParams scf; SCFAccelerators::Type accelerator = SCFAccelerators::Type::DIIS; };

    //! One post-convergence action, run IN ORDER on the converged calculation (the deck's `postSCF` list; each is `{"name":{params}}`).  Each reports
    //! itself to the console; its result summary goes in the revision record under `results`.  Prerequisites are checked by \c Resolve BEFORE the SCF starts.
    struct PostAction
    {
        enum class Kind { EstimateHubbardU, HubbardLoop, IndependentResponse, HubbardLinearResponse, HubbardFiniteDifference };
        Kind                kind = Kind::EstimateHubbardU;
        size_t              maxOuter = 20;        //!< hubbardLoop: outer ACBN0 steps allowed
        double              tolU_eV  = 1e-4;      //!< hubbardLoop: |dU| convergence per manifold (eV)
        int                 nq       = 1;         //!< independentResponse: the q-mesh is nq x nq x nq (must divide the k-mesh)
        std::vector<size_t> perturb;              //!< hubbardLinearResponse (empty = the manifolds carrying +U) / hubbardFiniteDifference (empty = {0}): manifold indices
        double              tol      = 1e-8;      //!< hubbardLinearResponse: Krylov relative residual
        size_t              maxIter  = 200, restart = 40;   //!< hubbardLinearResponse: the GMRES budget and restart length
        double              alpha    = 0.0;       //!< hubbardFiniteDifference: the perturbation step, HARTREE in RAM (`alpha_eV` in the file)
        bool operator==(const PostAction&) const = default;
    };

    std::string      structure;       //!< REQUIRED: a name in materials.json or molecules.json
    ivec3_t          kmesh{1,1,1};    //!< Monkhorst-Pack divisions of the Brillouin zone
    Basis            basis;
    State            state;
    SolidCalcOptions solid;           //!< the periodic-run options (used when \c structure is a cell)
    SCFParams        scf;             //!< the single-stage recipe (with \c solid.accelerator) ...
    std::vector<Stage> schedule;      //!< ... OR an annealed one; a deck gives one or the other (never both)
    std::vector<PostAction> postSCF;  //!< post-convergence actions, in order (empty = none)
    //! The stages the run executes: \c schedule, or the one stage (\c scf, \c solid.accelerator).
    std::vector<Stage> Stages() const
    { return schedule.empty() ? std::vector<Stage>{{scf,solid.accelerator}} : schedule; }
};
json ToJson(const RunSpec&);
//! Strict like the other readers: \c structure is required; an unknown top-level key (or any nested one) throws.
void FromJson(const json&, RunSpec&);
//! \brief Resolve \a spec against the structure files: returns the Material (cell + atoms + species) and FILLS \c spec.solid.Nelec /
//! \c species when they were left to derive, so \c ToJson(spec) afterwards is the complete record.  THROWS (listing the known names)
//! on an unknown structure, and for a molecule name (molecular decks are the next increment).  Also PRE-FLIGHT-VALIDATES the deck's `postSCF` against
//! the run (a response needs +U manifolds, a full k-mesh, the complex ansatz, a q-mesh dividing the k-mesh ...) so a deck that cannot work fails in
//! milliseconds, not after the SCF.
Materials::Material Resolve(RunSpec& spec);



//! \brief Apply `"a.b.c=value"` to \a deck: creates the intermediate objects; \c value is parsed as JSON (`1e-8`, `true`,
//! `[1,2]`, `"x"`) and, failing that, taken as a bare string (`tol.mode=fast`).  An array index is a number (`species.0.1=4`).
//! THROWS on a missing `=` or an empty path.  Type/key validity is checked later, by \c From* on the result.
void ApplySet(json& deck, std::string_view assignment);

//! \brief Atomically claim the next free revision `<dir>/<name>.r<NNN>.json` (exclusive create: parallel jobs cannot collide, a
//! number is never reused or overwritten) and return its path.  The file exists, empty, on return; \c WriteRevision fills it.
std::filesystem::path ClaimRevision(const std::filesystem::path& dir, const std::string& name);

//! \brief Provenance of a run: where its deck came from.  \c parent is the revision file the run started from (empty when it
//! started from a hand-written deck), \c commandLine the invocation, \c overrides the `--set` strings in the order applied.
struct Provenance
{
    std::filesystem::path inputDeck;                 //!< the deck as loaded (empty = none, all defaults)
    std::string           commandLine;
    std::string           restartedFrom;             //!< the revision whose saved state this run continues (set by Run from state.restartFrom)
    std::vector<std::string> overrides;
    std::string           policyResolved;            //!< the RESOLVED CP2K-deviation table (RunPolicy::Banner), for the record
    std::string           codeVersion;               //!< git hash (+"-dirty"), supplied by the caller
    std::vector<std::string> ignoredEnvironment;     //!< retired env variables found set (recorded, never honoured; from step 6a on)
    std::vector<std::string> activeEnvironment;      //!< INTERIM (until 6a retires the hooks): env overrides still HONOURED by the library, so the record is honest
};
//! Write \a resolved + header + provenance to \a path (a claimed revision).  Header: schema version, \c codeVersion, the input
//! deck's checksum.  The file is the COMPLETE record: loading it reproduces the run with no \c --set.
void WriteRevision(const std::filesystem::path& path, const json& resolved, const Provenance& prov);
//! Load a deck file (hand-written or a revision): strips the header/provenance, WARNS on stderr when its \c codeVersion differs
//! from \a currentCodeVersion, and returns the option payload.  THROWS if the schema version is newer than \c kSchemaVersion.
json LoadDeck(const std::filesystem::path& path, const std::string& currentCodeVersion);

//! \brief VET a deck without running it: resolve it (the postSCF pre-flight included), build the lattice, and build the basis (a stated trim applied; the
//! vet-stage trim loop is NOT run -- it is an SCF-adjacent cost), reporting what the run WOULD be.  THROWS on anything \c Run would refuse before its SCF (an
//! unknown structure, a basis file with no block for an element, a bad postSCF, a spherical view the span refuses ...).  Cheap: seconds, no SCF, no revision.
struct VetReport
{
    std::string              summary;   //!< one line: structure, atoms, Nelec, k-mesh, basis + functions, stages, manifolds, postSCF
    std::vector<std::string> notes;     //!< things that are not errors but a run would hit (a restartFrom state that does not exist yet)
};
VetReport Vet(RunSpec spec, const std::filesystem::path& outDir);

//! \brief What a deck run produced.  \c energy is present ONLY for a converged run (a failed run's last iterate is deliberately not
//! offered as "the energy" -- see \c SCFFailure::lastEnergy).
struct RunOutcome
{
    bool                   converged = false;
    std::optional<double>  energy;                //!< total energy (Ha) of the converged state
    std::filesystem::path  revision;              //!< the resolved deck this run wrote: the record
    std::string            summary;               //!< the one-line verdict (or the failure reason)
    struct PostResult { std::string action; bool ok=true; std::string summary; };
    std::vector<PostResult> postSCF;              //!< one per \c RunSpec::postSCF entry, in order (also written to the revision's `results`)
};
//! \brief RUN a deck: resolve it, CLAIM and WRITE its revision `<outDir>/<structure>.rNNN.json` (BEFORE the SCF, so even a crashed
//! run leaves its record), build the lattice, basis and \c SolidCalculation, converge the schedule, and report to the console.
//! The caller supplies the provenance (command line, \c --set strings, input deck); the code version is the caller's too.
RunOutcome Run(RunSpec spec, Provenance prov, const std::filesystem::path& outDir);

} // namespace qchem::deck
