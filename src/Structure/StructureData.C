// File: Structure/StructureData.C  Pre-defined STRUCTURES -- molecules and crystal cells -- read from
// src/Structure/Data/{molecules,materials}.json (D-STRUCTDATA, ruled 2026-10-02).
//
// WHAT THIS IS AND IS NOT.  It turns the STRUCTURE fields of the data files (lattice, atom basis, site spin
// flips, Cartesian atom positions) into the CONCRETE \c Molecule or \c UnitCell the caller asked for -- the caller
// knows what it needs, so there is no polymorphic return and no cast.  It deliberately knows NOTHING about
// pseudopotentials, valence counts, basis sets or any other recipe for a calculation: that is a property of
// (element, pseudopotential), not of the structure, and is read at the Calculation level.  The point is that
// every library's unit tests can build the shared geometries (the water, the MnO AFM-II cell, ...).
module;
#include <string>
#include <vector>
export module qchem.StructureData;
export import qchem.UnitCell;    // UnitCell, Bravais
export import qchem.Structure;   // Molecule

export namespace qchem::StructureData
{
//! What a registry name denotes.  The two data files share one namespace of names (a name in both is a data error).
enum class Kind { Molecule, Cell };

//! The kind of \a name; throws (listing every known name) on a miss.
Kind KindOf(const std::string& name);

//! The named entry of \c molecules.json (Cartesian a.u.).  Throws if \a name is unknown or is a CELL.
Molecule GetMolecule(const std::string& name);
//! Every molecule name, in file order.
std::vector<std::string> MoleculeNames();

//! The named crystal/box cell of \c materials.json, atoms and spin flips included (the primitive, possibly
//! superlattice-re-based cell).  \a a overrides the entry's primary lattice constant (a ladder or equation-of-
//! state scan), 0 = the banked value.  Throws if \a name is unknown or is a MOLECULE.
UnitCell GetCell(const std::string& name, double a=0.0);
//! Every cell name, in file order -- the GUI's pick-list.
std::vector<std::string> CellNames();

//! A parameterised BOX that is not in the data file: one \a element atom at the centre of an \a a-bohr cubic cell
//! (the named boxes in the file are instances of this at the banked sizes).
UnitCell AtomInBox(const std::string& element, double a);
//! ...and a homonuclear dimer at bond length \a d along x, centred; \a afm plants the +m/-m flip.
UnitCell DimerInBox(const std::string& element, double a, double d, bool afm=false);
} // namespace qchem::StructureData
