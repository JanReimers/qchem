// File:: BasisSet/Factory.C  Interfaces for various basis set factories.
module;
#include <nlohmann/json_fwd.hpp>
export module qchem.BasisSet.Atom.Factory;
export import qchem.BasisSet;
export import qchem.ElectronConfiguration;


export namespace qchem::BasisSet::Atom
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
    // (The old N3/N5 tiny fixed even-tempered pools (3/5 exponents, ltrim=0) were REMOVED: they are ill-
    //  conditioned and invalid for SCF, and were repeatedly mistaken for solver/code bugs.  Low is the
    //  smallest SCF-grade pool -- use it for the "small deterministic basis" unit tests that used N3/N5.)
    enum class BasisSetAccuracy {Low,Medium,High};
    std::string BasisSetAccuracyStrs[]={"Low    ", "Medium ","High   "};
    Real_BS* Factory(BasisSetAccuracy, Type,size_t Z);
    //! As above, but the irreps come from an explicit \a ec instead of Atom_EC(Z): \a Z still sets the
    //! exponent pool (the element), while \a ec sets WHICH angular irreps the basis carries.  Needed when
    //! the occupied config differs from the neutral element -- e.g. a pseudo-ion (F- closed p^6 vs neutral
    //! F open p^5); for a neutral atom pass Atom_EC(Z) and the result is identical to the Z overload.
    Real_BS* Factory(BasisSetAccuracy, Type, size_t Z, const ElectronConfiguration& ec);
}





