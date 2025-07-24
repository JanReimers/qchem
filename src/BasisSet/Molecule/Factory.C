// File:: BasisSet/Factory.C  Interfaces for various basis set factories.
module;
#include <nlohmann/json_fwd.hpp>
export module qchem.BasisSet.Molecule.Factory;
export import qchem.BasisSet;
export import qchem.Cluster;

export namespace BasisSetMolecule
{
    BasisSet* Factory(const nlohmann::json& js,const Cluster* cl);
}

