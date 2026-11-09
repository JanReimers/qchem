//! \file Structure/StructureData.C
//! \brief Pre-defined STRUCTURES -- molecules and crystal cells -- read from
//! src/Structure/Data/{molecules,materials}.json (D-STRUCTDATA, ruled 2026-10-02).
//!
//! WHAT THIS IS AND IS NOT.  It turns the STRUCTURE fields of the data files (lattice, atom basis, site spin
//! flips, Cartesian atom positions) into the CONCRETE \c Molecule or \c UnitCell the caller asked for -- the caller
//! knows what it needs, so there is no polymorphic return and no cast.  It deliberately knows NOTHING about
//! pseudopotentials, valence counts, basis sets or any other recipe for a calculation: that is a property of
//! (element, pseudopotential), not of the structure, and is read at the MoleculeCalculation level.  The point is that
//! every library's unit tests can build the shared geometries (the water, the MnO AFM-II cell, ...).
module;
#include <string>
#include <vector>
export module qchem.StructureData;
export import qchem.UnitCell;    // UnitCell, Bravais
export import qchem.Structure;   // Molecule

export namespace qchem::StructureData
{
//! \brief What a registry name denotes.  The two data files share one namespace of names (a name in both is a data error).
enum class Kind { Molecule, Cell };

//! \brief The kind of \a name; throws (listing every known name) on a miss.
Kind KindOf(const std::string& name);

//! \brief The named entry of \c molecules.json (Cartesian a.u.).  Throws if \a name is unknown or is a CELL.
Molecule GetMolecule(const std::string& name);
//! \brief Every molecule name, in file order.
std::vector<std::string> MoleculeNames();

//! \brief The named crystal/box cell of the SHIPPED \c materials.json, atoms and spin flips included (the primitive, possibly
//! superlattice-re-based cell).  \a a overrides the entry's primary lattice constant (a ladder or equation-of-
//! state scan), 0 = the banked value.  Throws if \a name is unknown or is a MOLECULE.
UnitCell GetCell(const std::string& name, double a=0.0);
//! \brief A cell from ONE entry in the \c materials.json schema (`lattice`, `atoms` with `el`/`frac`/optional `spin`; other
//! keys are the caller's), given as JSON TEXT -- the same parser \c GetCell uses, for an entry stated INLINE (a deck's object
//! `structure`, so a GUI can run a cell it built without a registry entry).  The lattice is the Bravais description (`type`,
//! `a`..`gamma`, `T`) or, exclusively, `vectors_bohr` [[a],[b],[c]] (the three lattice vectors, bohr).  \a who names it in
//! errors.  Throws on malformed JSON, an unknown lattice or atom key, or an unknown element / Bravais type.
UnitCell CellFromEntry(const std::string& entryJson, const std::string& who);
//! \brief Every SHIPPED cell name, in file order.  (A deck names shipped cells only: \c CellEntries lists the other layers.)
std::vector<std::string> CellNames();

//! \brief THE LAYERS of the materials (OpenWork §2, ruled 2026-11-09).  The SHIPPED file stays CURATED: its lattice constants
//! anchor tests, and promoting a material into it is a deliberate commit.  Imports go to the USER's own library (same schema), and
//! a workspace may add its own on top.  The library reads NO environment for them (controller ruling 2026-11-09, the getenv
//! audit): whoever knows where a library lives -- the GUI resolves ~/.qchem/materials.json -- passes its PATH in.  And a DECK
//! never names a non-shipped material: it inlines the entry (the deck is the one source of a run's numbers, D-ENV), so
//! \c KindOf / \c GetCell stay shipped-only.  A name lives in ONE layer: a user or workspace entry may not reuse a shipped (or
//! molecule) name, nor another layer's.
enum class Layer { Shipped, User, Workspace };
//! \brief A library file layered over the shipped one (materials.json schema).  A missing file is an empty layer.
struct LibraryFile
{
    Layer       layer;
    std::string path;
};
//! \brief One entry: its name, its layer and file, and the entry itself (JSON text, the materials.json schema).
struct CellEntry
{
    std::string name;
    Layer       layer;
    std::string file;
    std::string entryJson;
};
//! \brief The shipped entries, then each of \a layers' in order.  Throws, naming the file, on a name two layers hold.
std::vector<CellEntry> CellEntries(const std::vector<LibraryFile>& layers={});
//! \brief The named entry over the shipped file and \a layers; throws (listing the known names) on a miss.
CellEntry GetCellEntry(const std::string& name, const std::vector<LibraryFile>& layers={});
//! \brief Save \a entryJson under \a name in the library \a into (creating it), the file rewritten whole and in order.  Refuses
//! the shipped layer, a shipped or molecule name, a name another of \a others holds, an illegal name (empty, leading '_', a '/'),
//! an entry the parser refuses, and an existing name in \a into unless \a replace.
void SaveMaterial(const LibraryFile& into, const std::string& name, const std::string& entryJson, bool replace=false,
                  const std::vector<LibraryFile>& others={});

//! \brief A parameterised BOX that is not in the data file: one \a element atom at the centre of an \a a-bohr cubic cell
//! (the named boxes in the file are instances of this at the banked sizes).
UnitCell AtomInBox(const std::string& element, double a);
//! \brief ...and a homonuclear dimer at bond length \a d along x, centred; \a afm plants the +m/-m flip.
UnitCell DimerInBox(const std::string& element, double a, double d, bool afm=false);
} // namespace qchem::StructureData
