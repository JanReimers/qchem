// File:: BasisSet/Factory.C  Interfaces for various basis set factories.
module;
#include <nlohmann/json_fwd.hpp>
export module qchem.BasisSet1.Molecule.Factory;
export import qchem.BasisSet1;
export import qchem.Cluster;

export namespace BasisSet1::Molecule
{
    BasisSet1::BasisSet<double>* Factory(const nlohmann::json& js,const Cluster* cl);
}

