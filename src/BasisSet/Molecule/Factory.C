// File:: BasisSet/Factory.C  Interfaces for various basis set factories.
module;
#include <nlohmann/json_fwd.hpp>
export module qchem.BasisSet.Molecule.Factory;
export import qchem.BasisSet;
export import qchem.Cluster;

export namespace BasisSet::Molecule
{
    // Typo-proof orbital-basis selection: the caller names a basis, never a file path (the Factory owns
    // the path, via qchem.BasisSet.Molecule.BasisFiles).  Add a basis = one enum value + one map entry
    // in Imp/Factory.C.
    enum class Basis { DZVP, DZVP2, TZVP, ORB, ORB1 };

    // C++ entry point (compile-time typo-proof).
    Real_BS* Factory(Basis basis, const Cluster* cl);

    // Config-driven entry point: js["basis"] is a basis NAME (e.g. "dzvp"); an unknown name throws with
    // the list of valid names.  (Was js["filepath"] -- a raw path the caller had to build.)
    Real_BS* Factory(const nlohmann::json& js,const Cluster* cl);
}

