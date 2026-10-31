// File: BasisSet/Gaussian/Point/BasisFiles.C  One owner of the basis-set data directory.
//
// The basis-set data path comes from the build (BASISSET_DATA_PATH, -D on the compile line -- the one
// sanctioned preprocessor use, see CLAUDE.md).  Locating files is this module's job, so callers (the
// Factory, the fit-basis creation in the IBS) ask for a file by name and never handle the path -- and
// they cannot disagree about where the data lives.
module;
#include <string>
#include <vector>
export module qchem.BasisSet.Gaussian.Point.BasisFiles;

export namespace qchem::BasisSet::Gaussian
{
    // Absolute path (as a string, for the readers) to a basis-set data file, e.g. BasisFile("dzvp.bsd").
    std::string BasisFile(const std::string& filename);
    //! \brief THE BASIS NAMES (CleanCode D-BASIS-REGISTRY, 2026-10-30): a basis is named by its file STEM in the data directory --
    //! one source of truth, discovered, never registered.  Every \c .bsd stem there, sorted (fit and test files included: a name
    //! is a file, what it is FOR is the caller's business).
    std::vector<std::string> BasisNames();
    //! The file a basis NAME means: \a name (any case) + ".bsd" when that file is in the data directory.  THROWS otherwise,
    //! listing every name there -- a typo is said at once with the alternatives, never a missing-file error from a reader.
    std::string BasisFileOfName(const std::string& name);
    //! An identity checksum (64-bit FNV-1a, hex) of a basis file's CONTENTS (\a filename resolved as \c BasisFile does): what a run
    //! record keeps so "the same name" cannot hide an edited file.  THROWS if the file cannot be read.
    std::string BasisFileChecksum(const std::string& filename);
}
