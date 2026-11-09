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

//! \brief The named crystal/box cell of \c materials.json (any LAYER: shipped, user, workspace), atoms and spin flips included (the primitive, possibly
//! superlattice-re-based cell).  \a a overrides the entry's primary lattice constant (a ladder or equation-of-
//! state scan), 0 = the banked value.  Throws if \a name is unknown or is a MOLECULE.
UnitCell GetCell(const std::string& name, double a=0.0);
//! \brief A cell from ONE entry in the \c materials.json schema (`lattice`, `atoms` with `el`/`frac`/optional `spin`; other
//! keys are the caller's), given as JSON TEXT -- the same parser \c GetCell uses, for an entry stated INLINE (a deck's object
//! `structure`, so a GUI can run a cell it built without a registry entry).  The lattice is the Bravais description (`type`,
//! `a`..`gamma`, `T`) or, exclusively, `vectors_bohr` [[a],[b],[c]] (the three lattice vectors, bohr).  \a who names it in
//! errors.  Throws on malformed JSON, an unknown lattice or atom key, or an unknown element / Bravais type.
UnitCell CellFromEntry(const std::string& entryJson, const std::string& who);
//! \brief Every cell name over every LAYER (shipped, then the user's, then each workspace's), in file order -- the GUI's pick-list.
std::vector<std::string> CellNames();

//! \brief THE LAYERS of the cell registry (OpenWork §2, ruled 2026-11-09).  The SHIPPED file stays CURATED: its lattice constants
//! anchor tests, and promoting a material into it is a deliberate commit.  Imports go to the USER's own library, same schema; a
//! workspace may add its own on top.  A name lives in ONE layer: a user or workspace entry may not reuse a shipped (or molecule)
//! name, and the user's and a workspace's may not collide -- one name, one cell, whichever process reads it.
enum class Layer { Shipped, User, Workspace };
//! \brief One registry entry: its name, its layer and file, and the entry itself (JSON text, the materials.json schema).
struct CellEntry
{
    std::string name;
    Layer       layer;
    std::string file;
    std::string entryJson;
};
//! \brief The user's library: \c $QCHEM_USER_MATERIALS when set (tests, a GUI's sandbox), else \c ~/.qchem/materials.json.  It
//! need not exist (an empty layer).
std::string UserLibraryPath();
//! \brief Layer the workspace library \a path (materials.json schema) on top, for this process.  Throws on a name another layer
//! holds, or on an entry the parser refuses -- at the add, not at the first run that names it.
void AddWorkspaceLibrary(const std::string& path);
//! \brief Every cell entry over every layer, shipped first.  Re-reads the user's file when it changed on disk.
std::vector<CellEntry> CellEntries();
//! \brief The named cell entry, any layer; throws (listing every known name) on a miss.
CellEntry GetCellEntry(const std::string& name);
//! \brief Save \a entryJson under \a name in the USER's library (creating it), the file rewritten whole and in order.  Refuses a
//! shipped, molecule or workspace name, an illegal name (empty, leading '_', a '/'), an entry the parser refuses, and an existing
//! user name unless \a replace.
void SaveUserMaterial(const std::string& name, const std::string& entryJson, bool replace=false);

//! \brief A parameterised BOX that is not in the data file: one \a element atom at the centre of an \a a-bohr cubic cell
//! (the named boxes in the file are instances of this at the banked sizes).
UnitCell AtomInBox(const std::string& element, double a);
//! \brief ...and a homonuclear dimer at bond length \a d along x, centred; \a afm plants the +m/-m flip.
UnitCell DimerInBox(const std::string& element, double a, double d, bool afm=false);
} // namespace qchem::StructureData
