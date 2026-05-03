// File:: BasisSet/Factory.C  Interfaces for various basis set factories.
module;
#include <nlohmann/json_fwd.hpp>
export module qchem.BasisSet.Atom.Factory;
export import qchem.BasisSet;
export import qchem.Cluster;
export import qchem.Symmetry.ElectronConfiguration;

export namespace BasisSetAtom
{
    // enum class Type {Slater,Gaussian,BSpline6,BSpline9,BSpliner6,BSpliner9,BSpline16,BSpline19,BSpline1r6,BSpline1r9,Slater_RKB,Gaussian_RKB,BSpline_RKB};
    enum class Type {Slater,Gaussian,BSpline6,BSpline9,BSpliner6,BSpliner9,Slater_RKB,Gaussian_RKB,BSpline_RKB};
    BasisSet* Factory(Type,const nlohmann::json& js,const ElectronConfiguration& ec);
    BasisSet* Factory(Type,const nlohmann::json& js,size_t Z); //Use atomic number Z to deduce the ElectronConfiguration
    BasisSet* Factory(const nlohmann::json& js,size_t Z); //Type is in the json object
}





