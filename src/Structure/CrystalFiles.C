//! \file Structure/CrystalFiles.C
//! \brief Crystal structure FILES (CIF, POSCAR, extended xyz) <-> one \c materials.json entry (OpenWork §2, design ruled
//! by the user 2026-11-09).
//!
//! THE ENTRY IS THE INTERNAL FORMAT.  A reader is a TRANSLATOR: it emits one entry in the \c materials.json schema
//! (`comment`, `lattice: {vectors_bohr}`, `atoms: [{el, frac}]`, `provenance`) as JSON text, and the cell is then built
//! by the one parser every entry goes through (\c StructureData::CellFromEntry, \c Materials::Inline) -- so an imported
//! structure IS an entry, never a second path to a cell.  The writers go the other way, for the oracles.
//!
//! WHAT AN ENTRY CANNOT SAY IS REFUSED BY NAME, NEVER DROPPED: a CIF site with partial occupancy, a disorder group, two
//! sites on one point, a POSCAR without species names, an extended xyz that is not periodic in all three directions.
//! Displacement parameters (ADPs) are ignored: they describe the experiment, not the ordered crystal.  Selective-dynamics
//! flags and velocities in a POSCAR are ignored too (relaxation inputs, not structure).
//!
//! A CIF's stated symmetry is kept as provenance and CHECKED against the library's own detection
//! (\c SpaceGroup::Detect on the cell that is stored): a disagreement is a WARNING in \c Imported::warnings and in the
//! entry's `provenance.spaceGroup`, never a refusal (the file may round its coordinates; the user decides).
//!
//! The first version stores `vectors_bohr` (Bravais recognition is later work).  A CENTRED CIF lattice (F, I, C, R, ...)
//! is stored as its PRIMITIVE cell, built from the centring translations the CIF itself lists, and the provenance says
//! so: the conventional cell would run the same crystal at n times the cost.
module;
#include <string>
#include <vector>
export module qchem.CrystalFiles;

export namespace qchem::CrystalFiles
{
//! The three exchange formats.  Extended xyz = the ASE/extxyz convention (`Lattice="..."` and `Properties=species:S:1:pos:R:3`).
enum class Format { CIF, POSCAR, ExtendedXYZ };

//! The format a FILE NAME says: `.cif`; `.xyz` / `.extxyz`; `POSCAR*` / `CONTCAR*` / `.vasp` / `.poscar`.  Throws naming them.
Format FormatOf(const std::string& fileName);
//! The format's name, as the provenance and the errors spell it ("CIF", "POSCAR", "extended xyz").
std::string NameOf(Format);

//! \brief What an import hands back.
struct Imported
{
    std::string              entry;     //!< ONE materials.json entry, JSON text: comment, lattice{vectors_bohr}, atoms, provenance
    std::string              name;      //!< a suggested entry name (the file stem, made a legal name); the caller may rename
    std::vector<std::string> warnings;  //!< what the import checked and found wanting (also in the entry's provenance)
};

//! \brief Translate \a text in format \a f into one entry.  \a source names the file in the provenance and the errors;
//! \a date (YYYY-MM-DD) is the import date the provenance records, empty = today.  Throws \c std::invalid_argument
//! naming what the entry cannot represent, or what the text gets wrong.
Imported Read(const std::string& text, Format f, const std::string& source, const std::string& date="");
//! \brief Read the file at \a path, its format from its name (\c FormatOf).
Imported ReadFile(const std::string& path, const std::string& date="");

//! \brief Write the entry \a entryJson (any materials.json entry: Bravais or vectors_bohr) as format \a f; \a name heads the
//! file (CIF data block, POSCAR comment).  Site spins have no slot in POSCAR or CIF: they are written into the file's
//! comment (POSCAR line 1, a CIF comment) in site order, so they are visible -- the oracle's own magnetic input carries them.
std::string Write(const std::string& entryJson, Format f, const std::string& name);
} // namespace qchem::CrystalFiles
