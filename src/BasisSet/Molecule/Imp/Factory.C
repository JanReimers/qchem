// File::BasisSet/Molecule/Factory.C  Factory function for molecular basis sets.
module;
#include <cassert>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

module qchem.BasisSet.Molecule.Factory;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Readers.Gaussian94;
import qchem.BasisSet.Molecule.PolarizedGaussian;

namespace BasisSetMolecule
{
    BasisSet* Factory(const nlohmann::json& js,const Cluster* cl)
    {
        std::string filepath=js["filepath"].template get<std::string>();
        PolarizedGaussian::Gaussian94Reader reader(filepath);
        return new PolarizedGaussian::BasisSet(&reader,cl);  
    }
}

