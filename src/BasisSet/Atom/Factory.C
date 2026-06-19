// File:: BasisSet/Factory.C  Interfaces for various basis set factories.
module;
#include <nlohmann/json_fwd.hpp>
export module qchem.BasisSet.Atom.Factory;
export import qchem.BasisSet;
export import qchem.ElectronConfiguration;


export namespace BasisSet::Atom
{
    enum class Type {Slater,Gaussian,BSpline6,BSpliner6,Gaussian_RKB,Slater_RKB};
    // Use these for custom basis sets.
    Real_BS* Factory(const nlohmann::json&,const ElectronConfiguration&);
    Real_BS* Factory(const nlohmann::json&,size_t Z); //Use atomic number Z to deduce the ElectronConfiguration
    // The object below is to define restricted collection of atom basis sets so that the global DB_Cache to not get too large
    // when running 100s of unit tests.  But these 'canned' versions are also generally usefull for production work.
    //
    // High:   As big and dense as you can get without becoming numercially unstable.
    // Medium: High with one smallest exponent and ~10-20% of the largest exponents removed. 
    //         Giving up core cusp but not exponent density
    // Low:    High with every other expoenent removed.  Giving up density but not core cusp.
    //
    enum class BasisSetAccuracy {N3,N5,Low,Medium,High};
    std::string BasisSetAccuracyStrs[]={"N3     ","N5     ","Low    ", "Medium ","High   "};
    Real_BS* Factory(BasisSetAccuracy, Type,size_t Z);
}





