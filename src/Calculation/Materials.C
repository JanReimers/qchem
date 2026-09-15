// File: Calculation/Materials.C  Pre-defined materials and molecules -- DATA, read from
// src/Calculation/Data/{materials,molecules}.json (doc/OpenWork.md row MD, 2026-09-15).
module;
#include <memory>
#include <string>
#include <vector>
#include <utility>
export module qchem.Materials;
export import qchem.UnitCell;    // UnitCell, Bravais (the lattice TYPES live in qcStructure; the materials do not)
export import qchem.Structure;   // Molecule
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
    std::vector<std::pair<std::string,int>>  species;   //!< {element, valence electrons} -- \c SolidCalcOptions::species verbatim
    //! Valence electrons per cell, DERIVED from the atoms and the species list (never stored: a stored
    //! count that disagreed with the species list would be exactly the defect a data file invites).
    int Nelec() const;
};

//! The named entry of \c materials.json; throws with the list of known names on a miss.  \a a overrides
//! the entry's primary lattice constant (a ladder or an equation-of-state scan), 0 = the banked value.
Material Get(const std::string& name, double a=0.0);
//! Every material name in the data file, in file order -- the GUI's pick-list.
std::vector<std::string> Names();

//! A parameterised BOX that is not in the data file: one \a Z pseudo-atom at the centre of an \a a-bohr
//! cubic cell (the named boxes in the file are instances of this at the banked sizes).
Material AtomInBox(const std::string& element, int valence, double a);
//! ...and a homonuclear dimer at bond length \a d along x, centred; \a afm plants the +m/-m flip.
Material DimerInBox(const std::string& element, int valence, double a, double d, bool afm=false);

//! The named entry of \c molecules.json (Cartesian a.u.); throws with the known names on a miss.
Molecule GetMolecule(const std::string& name);
std::vector<std::string> MoleculeNames();

} // namespace qchem::Materials
