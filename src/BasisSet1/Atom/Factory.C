// File:: BasisSet/Factory.C  Interfaces for various basis set factories.
module;
#include <nlohmann/json_fwd.hpp>
export module qchem.BasisSet1.Atom.Factory;
export import qchem.BasisSet1;
export import qchem.Symmetry.ElectronConfiguration;

export namespace BasisSet1::Atom
{
    // enum class Type {Slater,Gaussian,BSpline6,BSpline9,BSpliner6,BSpliner9,BSpline16,BSpline19,BSpline1r6,BSpline1r9,Slater_RKB,Gaussian_RKB,BSpline_RKB};
    enum class Type {Slater,Gaussian,BSpline6,Gaussian_RKB,Slater_RKB,Gaussian2,Slater2,BSpline6_2};
    BasisSet1::Real_BS* Factory(const nlohmann::json& js,const ElectronConfiguration& ec);
    BasisSet1::Real_BS* Factory(const nlohmann::json& js,size_t Z); //Use atomic number Z to deduce the ElectronConfiguration
}





