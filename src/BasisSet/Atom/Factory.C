// File:: BasisSet/Factory.C  Interfaces for various basis set factories.
module;
#include <nlohmann/json_fwd.hpp>
export module qchem.BasisSet.Atom.Factory;
export import qchem.BasisSet;
export import qchem.Cluster;
export import qchem.Symmetry.ElectronConfiguration;

export namespace BasisSetAtom
{
    enum class Type {Slater,Gaussian,BSpline,Slater_RKB,Gaussian_RKB,BSpline_RKB};
    BasisSet* Factory(Type,const nlohmann::json& js,const ElectronConfiguration& ec);
    BasisSet* Factory(Type,const nlohmann::json& js,size_t Z); //Use atomic number Z to deduce the ElectronConfiguration
    BasisSet* Factory(const nlohmann::json& js,size_t Z); //Type is in the json object
}


