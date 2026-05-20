// File::BasisSet/Molecule/Factory.C  Factory function for molecular basis sets.
module;
#include <cassert>
#include <nlohmann/json.hpp>

module qchem.BasisSet.Molecule.Factory;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Readers.Gaussian94;
import qchem.BasisSet.Molecule.PolarizedGaussian;
import qchem.BasisSet.DB_Cache;

using json = nlohmann::json;

namespace BasisSet::Molecule
{
    BasisSet<double>* Factory(const nlohmann::json& js,const Cluster* cl)
    {
        if (::BasisSet::theGlobalCache==0)
        ::BasisSet::theGlobalCache=new ::BasisSet::IntegralsCache_RAM<double>(true);     

        std::string filepath=js["filepath"].template get<std::string>();
        PolarizedGaussian::Gaussian94Reader reader(filepath);
        return new PolarizedGaussian::BasisSet(&reader,cl);  
    }
}

