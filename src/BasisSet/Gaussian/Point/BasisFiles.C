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

    //! \brief WHAT A BASIS FILE IS FOR, as the file itself states it (OpenWork "A basis file states its KIND", lane G 2026-11-21).
    //! Every shipped file says it in its preamble (before the first lone `!`, which the Gaussian94 reader skips) on ONE line:
    //!     KIND all-electron                    an orbital basis for every electron
    //!     KIND pp-valence <pseudopotential set> an orbital basis for the valence electrons of that set (pseudopotentials.json name)
    //!     KIND fit-coulomb | fit-exchange      an auxiliary (fitting) basis, never an orbital menu entry
    //! The property is the FILE's, so a front end (qchem.BasisMenu) and a run judge the file it reads, never a list kept beside it.
    //! A pp-valence file may also declare each block's q (`PP-PROVENANCE`, \c CheckPPProvenance); the two lines say different
    //! things (which set / which variant of it) and are read separately.
    enum class BasisKind { Unstated, AllElectron, PPValence, CoulombFit, ExchangeFit };
    struct BasisFileKind
    {
        BasisKind   kind = BasisKind::Unstated;  //!< Unstated: a file with no KIND line (a user's own .bsd; no shipped file -- a test holds it)
        std::string ppSet;                       //!< the pseudopotential set a PPValence file pairs with; "" otherwise
    };
    //! The KIND line of the basis file \a filename (resolved as \c BasisFile does).  THROWS if the file cannot be read, or on a
    //! malformed KIND line (an unknown word, pp-valence without its set, a set after any other kind, two KIND lines).
    BasisFileKind KindOfBasisFile(const std::string& filename);
    //! The KIND line's spelling of \a k ("all-electron", "pp-valence", "fit-coulomb", "fit-exchange"; "" for Unstated).
    std::string KindName(BasisKind k);
}
