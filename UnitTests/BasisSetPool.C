// File: UnitTests/BasisSetPool.C  Define a restricted set of basis sets for all unit tests.  
// The object here is to define restrict collection of atom basis sets so that the global DB_Cache to not get too large.
module;
#include <string>
export module qchem.Unittests.BasisSetPool;
export import qchem.BasisSet.Atom.Factory; //To get BS types.

//
// High:   As big and dense as you can get without becoming numercially unstable.
// Medium: High with one smallest exponent and ~10-20% of the largest exponents removed. 
//         Giving up core cusp but not exponent density
// Low:     High with every other expoenent removed.  Giving up density but not core cusp.
//
export enum class BasisSetAccuracy {N3,N5,Low,Medium,High};
export std::string BasisSetAccuracyStrs[]={"N3     ","N5     ","Low    ", "Medium ","High   "};
export BasisSet::Real_BS* PoolFactory(BasisSetAccuracy acc, BasisSet::Atom::Type type,size_t Z);
