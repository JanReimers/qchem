// File: BasisSet/Molecule/Imp/BasisFiles.C
module;
#include <string>
#include <filesystem>
#ifndef BASISSET_DATA_PATH
#error "BASISSET_DATA_PATH must be defined by CMake"
#endif
module qchem.BasisSet.Molecule.BasisFiles;

namespace BasisSet::Molecule
{
    static const std::filesystem::path theDataDir = BASISSET_DATA_PATH;

    std::string BasisFile(const std::string& filename)
    {
        return (theDataDir / filename).string();
    }
}
