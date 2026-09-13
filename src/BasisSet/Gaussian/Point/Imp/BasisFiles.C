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
        return (theDataDir / filename).string();
    }
}
