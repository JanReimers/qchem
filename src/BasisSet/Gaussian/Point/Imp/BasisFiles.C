// File: BasisSet/Gaussian/Point/Imp/BasisFiles.C
module;
#include <string>
#include <filesystem>
#ifndef BASISSET_DATA_PATH
#error "BASISSET_DATA_PATH must be defined by CMake"
#endif
module qchem.BasisSet.Gaussian.Point.BasisFiles;

namespace qchem::BasisSet::Gaussian
{
    static const std::filesystem::path theDataDir = BASISSET_DATA_PATH;

    std::string BasisFile(const std::string& filename)
    {
        // A bare name lives in the shipped data directory; anything with a directory part (the deck's \c basis.file
        // override, an absolute or cwd-relative path) is used exactly as given.
        if (filename.find('/')!=std::string::npos) return filename;
        return (theDataDir / filename).string();
    }
}
