// File: Calculation/Deck.C  The INPUT DECK: a run as one serialized, versioned JSON value (D-ENV step 6;
// doc/Records/EnvKnobInventory.md §9).
//
// THE RULES (user, 2026-10-03):
//   * The deck is the ONLY way to set a tier-1/2 value.  The environment is not a second path.
//   * The file the user wrote is never the record of what ran.  Every run writes a numbered REVISION
//     `<name>.r<NNN>.json`: the deck + every `--set` + every default filled in + a header + a provenance block.
//     Re-running that file reproduces the run; re-running it with one `--set` is the intentional change.
//   * A reader REJECTS an unknown key (a typo in a deck must not silently run the default), and a missing key keeps
//     the field's default (so an old deck keeps working when a field is added).  Units are the library's: atomic
//     units, and every key that carries one says so in its name where it is not obvious (`U_Ha`).
//
// This unit is the DATA LAYER: the typed option structs <-> JSON, the dotted-path `--set`, the revision claim.
module;
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>
#include <nlohmann/json.hpp>
export module qchem.Deck;

import qchem.SolidCalculation;     // SolidCalcOptions
import qchem.SCFParams;            // SCFParams
import qchem.Mesh;                 // MeshParams
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
//! (derived by the facade from the structure).  A manifold's U is \c U_Ha: HARTREE, the library's unit -- no hidden eV conversion.
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
    std::string      structure;       //!< REQUIRED: a name in materials.json or molecules.json
    SolidCalcOptions solid;           //!< the periodic-run options (used when \c structure is a cell)
    SCFParams        scf;
};
json ToJson(const RunSpec&);
//! Strict like the other readers: \c structure is required; an unknown top-level key (or any nested one) throws.
void FromJson(const json&, RunSpec&);
//! \brief Resolve \a spec against the structure files: returns the Material (cell + atoms + species) and FILLS \c spec.solid.Nelec /
//! \c species when they were left to derive, so \c ToJson(spec) afterwards is the complete record.  THROWS (listing the known names)
//! on an unknown structure, and for a molecule name (molecular decks are the next increment).
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
    std::vector<std::string> overrides;
    std::string           codeVersion;               //!< git hash (+"-dirty"), supplied by the caller
    std::vector<std::string> ignoredEnvironment;     //!< retired env variables found set (recorded, never honoured)
};
//! Write \a resolved + header + provenance to \a path (a claimed revision).  Header: schema version, \c codeVersion, the input
//! deck's checksum.  The file is the COMPLETE record: loading it reproduces the run with no \c --set.
void WriteRevision(const std::filesystem::path& path, const json& resolved, const Provenance& prov);
//! Load a deck file (hand-written or a revision): strips the header/provenance, WARNS on stderr when its \c codeVersion differs
//! from \a currentCodeVersion, and returns the option payload.  THROWS if the schema version is newer than \c kSchemaVersion.
json LoadDeck(const std::filesystem::path& path, const std::string& currentCodeVersion);

} // namespace qchem::deck
