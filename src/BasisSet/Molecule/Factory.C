// File: BasisSet/Molecule/Factory.C  Public entry point for building a molecular orbital basis set.
//
// A molecular basis is chosen along THREE ORTHOGONAL axes; every combination is valid.  Mind the word
// "basis": in code a `BasisSet` is an OBJECT (what this Factory returns) -- it is NOT what the first axis
// selects.  The three axes are:
//
//   1. BasisSetData -- WHICH canned data set.  A named, externally generated, optimized table of (mostly
//      radial) exponents + contraction coefficients, shipped as a data file (DZVP, TZVP, ...).  This is the
//      "basis set" a chemist names; it says nothing about how the functions are represented or integrated.
//   2. Engine       -- HOW the integrals over that basis are computed: MnD (McMurchie-Davidson, our native
//      code) or LibCint (the external library).  Same physics, different machinery -- a cross-check /
//      performance axis.  (LibCint is HF-only today: no DFT 3-centre fit.)
//   3. Angular      -- the ANGULAR representation of the in-code basis functions: Cartesian Gaussians
//      ((l+1)(l+2)/2 per shell, carrying l>0 contaminants) or real Spherical harmonics (2l+1 per shell).
//
// Config (json) keys -- "basis" is required, the rest default:
//   { "basis":   "dzvp",                  // BasisSetData, by name
//     "engine":  "mnd" | "libcint",       // Engine,  default "mnd"
//     "angular": "cartesian" | "spherical" }  // Angular, default "cartesian"
module;
#include <nlohmann/json_fwd.hpp>
export module qchem.BasisSet.Molecule.Factory;
export import qchem.BasisSet;
export import qchem.Structure;

export namespace qchem::BasisSet::Molecule
{
    // Axis 1 -- the canned, optimized, externally generated (mostly-radial) basis-set DATA FILE.  The
    // "...Data" suffix disambiguates from the BasisSet OBJECT in code.  Add one = one enum value + one
    // file-map entry in Imp/Factory.C.
    enum class BasisSetData { DZVP, DZVP2, TZVP, ORB, ORB1, SIPP };  // SIPP: valence-only Si Gaussian (PP validation)

    // Axis 2 -- the integral engine (orthogonal to the data set and the angular representation).
    enum class Engine { MnD, LibCint };

    // Axis 3 -- the angular representation of the in-code basis functions (orthogonal to the other two).
    enum class Angular { Cartesian, Spherical };

    // C++ entry point (compile-time typo-proof).  The three axes are independent; all combinations are valid.
    Real_BS* Factory(BasisSetData data, const Structure* cl,
                     Engine engine = Engine::MnD, Angular angular = Angular::Cartesian);

    // Config-driven entry point (see the key table above).  An unknown "basis"/"engine"/"angular" value, or
    // a missing "basis", throws with the list of valid values.
    Real_BS* Factory(const nlohmann::json& js, const Structure* cl);
}
