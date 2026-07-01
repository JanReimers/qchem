// File A_PP.C  Pseudopotential atom regression test (facade-driven).
//
// A pseudo-atom is the all-electron-basis atom with (a) the nuclear attraction replaced by the GTH local +
// KB-separable nonlocal pseudopotential, (b) only the Zion valence electrons (PseudoAtom_EC).  Migrated off
// QchemTester onto qchem::AtomCalculation via {.pseudopotential=true}: valence electrons = Z - charge, the
// element is looked up from Z.  "Did the physics move" anchor (no oracle / Converged() guard -- the virial
// criterion does not apply to a PP and the absolute value is basis/fitting-limited).
#include "gtest/gtest.h"
import qchem.AtomCalculation;        // AtomCalculation, AtomType, BasisSetAccuracy
import qchem.SCFIterator;            // SCFParams
using namespace qchem;
using enum BasisSetAccuracy;

// Silicon, 4 valence electrons (GTH-LDA q4).  The KB nonlocal projectors lift the over-bound local-only s
// state back toward the all-electron valence; total energy + valence charge pinned (Slater/Medium).
TEST(Si_PP_U, Medium)
{
    const int Z=14, val=4;
    AtomCalculation calc(Z, Z-val, {.type = AtomType::Slater, .accuracy = Medium, .pseudopotential = true},
        {.NMaxIter = 120, .MinΔρ = 1e-7, .MinΔFD = 1e-7, .MinVirial = 1e10, .MinFD = 1e-7, .StartingRelaxRo = 0.5, .MergeTol = 1e-7, .Verbose = true});  // virial off (N/A to PP)

    EXPECT_NEAR(calc.Energy(),      -3.336910601, 1e-6);   // pinned regression anchor (Slater/Medium)
    EXPECT_NEAR(calc.TotalCharge(),  4.0,         1e-9);   // valence electron count
}
