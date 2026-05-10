// File::BasisSet/Molecule/Factory.C  Factory function for molecular basis sets.
module;
#include <cassert>
#include <nlohmann/json.hpp>

module qchem.BasisSet1.Molecule.Factory;
import qchem.BasisSet1.Molecule.PolarizedGaussian.Internal.Readers.Gaussian94;
import qchem.BasisSet1.Molecule.PolarizedGaussian;
import qchem.BasisSet1.DB_Cache;

using json = nlohmann::json;

namespace BasisSet1::Molecule
{
    BasisSet<double>* Factory(const nlohmann::json& js,const Cluster* cl)
    {
        if (BasisSet1::theGlobalCache==0)
        BasisSet1::theGlobalCache=new BasisSet1::IntegralsCache_RAM<double>(true);     

        std::string filepath=js["filepath"].template get<std::string>();
        PolarizedGaussian::Gaussian94Reader reader(filepath);
        return new PolarizedGaussian::BasisSet(&reader,cl);  
    }
}

