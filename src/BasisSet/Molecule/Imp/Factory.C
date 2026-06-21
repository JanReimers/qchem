// File::BasisSet/Molecule/Factory.C  Factory function for molecular basis sets.
module;
#include <cassert>
#include <nlohmann/json.hpp>

module qchem.BasisSet.Molecule.Factory;
import qchem.BasisSet.Molecule.Readers.Gaussian94;
import qchem.BasisSet.Molecule.PolarizedGaussian;

using json = nlohmann::json;

namespace BasisSet::Molecule
{
    BasisSet<double>* Factory(const nlohmann::json& js,const Cluster* cl)
    {
        std::string filepath=js["filepath"].template get<std::string>();
        Gaussian94Reader reader(filepath);
        return new PolarizedGaussian::BasisSet(&reader,cl);
    }
}

