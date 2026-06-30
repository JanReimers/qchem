// File: BasisSet/Molecule/BasisFiles.C  One owner of the basis-set data directory.
//
// The basis-set data path comes from the build (BASISSET_DATA_PATH, -D on the compile line -- the one
// sanctioned preprocessor use, see CLAUDE.md).  Locating files is this module's job, so callers (the
// Factory, the fit-basis creation in the IBS) ask for a file by name and never handle the path -- and
// they cannot disagree about where the data lives.
module;
#include <string>
export module qchem.BasisSet.Molecule.BasisFiles;

export namespace qchem::BasisSet::Molecule
{
    // Absolute path (as a string, for the readers) to a basis-set data file, e.g. BasisFile("dzvp.bsd").
    std::string BasisFile(const std::string& filename);
}
