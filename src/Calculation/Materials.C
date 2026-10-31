// File: Calculation/Materials.C  Pre-defined materials and molecules -- DATA, read from
// src/Structure/Data/materials.json (structure: qchem.StructureData; pseudopotential vocabulary: the named set of qchem.Pseudopotential.Set) (doc/OpenWork.md row MD, 2026-09-15).
module;
#include <memory>
#include <string>
#include <vector>
#include <utility>
export module qchem.Materials;
export import qchem.StructureData;   // the STRUCTURE half (GetMolecule/GetCell/...) lives in qcStructure
export import qchem.UnitCell;    // UnitCell, Bravais (the lattice TYPES live in qcStructure; the materials do not)
export import qchem.Structure;
import qchem.Types;

export namespace qchem::Materials
{

//! \brief A crystal: a cell with its atom basis and magnetic decoration, plus the pseudopotential
//! vocabulary that goes with it -- and NOTHING about a run.  The k-mesh, XC grid, kT, basis set and
//! seed are the calculation's axes (\c SolidCalcOptions), not the material's.
//!
//! A material is a USE of a structure, so it lives at the Calculation level as data rather than in
//! \c qcStructure (user ruling 2026-09-15); the same list is what the GUI offers as its pick-list.
struct Material
{
    std::string                              name;
    std::shared_ptr<UnitCell>                cell;      //!< the (primitive, possibly re-based) cell with its atoms + spin flips
    //! The pseudopotential SET this material names (src/Pseudopotential/Data/pseudopotentials.json); empty for a box built from an explicit valence.
    std::string                              pp;
    //! {element, valence electrons} for the elements of the cell in order of first appearance, taken from \c pp --
    //! \c SolidCalcOptions::species verbatim.
    std::vector<std::pair<std::string,int>>  species;
    //! Valence electrons per cell, DERIVED from the atoms and the species list (never stored: a stored
    //! count that disagreed with the species list would be exactly the defect a data file invites).
    int Nelec() const;
};

//! The named entry of \c materials.json; throws with the list of known names on a miss.  \a a overrides
//! the entry's primary lattice constant (a ladder or an equation-of-state scan), 0 = the banked value.
Material Get(const std::string& name, double a=0.0);
//! \brief A material stated INLINE, not in the data file (a deck's object `structure`): one \c materials.json entry -- `lattice`
//! (Bravais or `vectors_bohr`), `atoms`, optional `comment` -- as JSON TEXT, with an optional `name` (a file stem; it names the
//! run's revisions; default "inline") and an optional `pp` (the pseudopotential SET; default \a pp, the run's).  The cell is built
//! by the data file's own parser, so an inline copy of a named entry IS that material.  Throws on an unknown key, a bad name, no
//! set, or anything \c StructureData::CellFromEntry refuses.
Material Inline(const std::string& entryJson, const std::string& pp);
//! Every material name in the data file, in file order -- the GUI's pick-list.
std::vector<std::string> Names();

//! A parameterised BOX that is not in the data file: one \a Z pseudo-atom at the centre of an \a a-bohr
//! cubic cell (the named boxes in the file are instances of this at the banked sizes).
Material AtomInBox(const std::string& element, int valence, double a);
//! ...and a homonuclear dimer at bond length \a d along x, centred; \a afm plants the +m/-m flip.
Material DimerInBox(const std::string& element, int valence, double a, double d, bool afm=false);

} // namespace qchem::Materials
