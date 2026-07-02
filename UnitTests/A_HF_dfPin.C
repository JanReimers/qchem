// File A_HF_dfPin.C  Tight self-pinned REGRESSION anchors for open-shell d/f atoms, across the three
// non-relativistic atomic basis families: BSpline6, Gaussian (SG) and Slater (SL).
//
// WHY a separate file from A_HF_P.  A_HF_P bounds the signed error vs the Saito near-HF-limit table -- a
// CORRECTNESS check, but only as tight as Saito's own precision plus the basis truncation (~1e-6 rel).  A
// code change that shifts an energy BELOW that floor slips through unnoticed.  These tests add a second,
// much tighter layer: they pin the code's OWN converged energy (a "did E move" self-anchor, exactly like
// the molecular M_* tests), so ANY drift in the ERI / SCF path is caught immediately -- decoupled from how
// physically accurate the number is.  This is the guard against "a change silently breaking the hard-won
// open-shell d/f energies".
//
// Two assertions per (family, element):
//   (1) regression  : E == the pinned self-value to ~1e-9 absolute        [the tight new guard]
//   (2) convergence : IsConverged()
// The looser correctness-vs-Saito bound already lives in A_HF_P and is not duplicated here.
//
// CROSS-FAMILY INVARIANCE.  HF total energy is a physical (basis-independent) quantity, so the three
// families' pinned values must agree to the HF-limit completeness bound.  That single check (below) means
// one hand-won Saito comparison + "the three still agree" replaces three independently hand-won numbers:
// if a family's basis data is ever mis-edited, it desyncs from the other two here.
#include "gtest/gtest.h"
#include <string>
#include <vector>
#include <iostream>
import qchem.AtomCalculation;        // AtomCalculation, AtomType, BasisSetAccuracy, Model, Pol
import qchem.SCFIterator;            // SCFParams
import qchem.Unittests.TestUtils;    // RelativeHFError (Saito oracle) -- used only for the harvest print
using namespace qchem;
using enum BasisSetAccuracy;         // High, Medium, Low

namespace {

//! One open-shell HF atom self-pin: family + accuracy + element + the pinned converged energy + SCF knobs.
//! pin==0.0 means "harvest mode": the case only prints E (no regression assert) so the pin can be read off.
struct PinCase { AtomType type; BasisSetAccuracy acc; int Z; double pin; double atol; };

//! The per-(family,High) SCF knobs, copied verbatim from A_HF_P's High instantiations so these anchors
//! converge on the same footing as the correctness tests.  Z-scaled exactly as MakeParams there.
static SCFParams HFParams(const PinCase& c)
{
    struct K { std::size_t nmax; double droS, dfd, vir, fdS; };
    K k = c.type==AtomType::Gaussian ? K{50,1e-5,1e-7,1e-5,   1e-6}
        : c.type==AtomType::Slater   ? K{32,1e-5,1e-7,1e-6,   1e-6}
        :                              K{50,1e-7,1e-7,2.5e-12,1e-7};   // BSpline6
    return {.NMaxIter = k.nmax, .MinΔρ = c.Z*k.droS, .MinΔFD = k.dfd, .MinVirial = k.vir,
            .MinFD = c.Z*k.fdS, .StartingRelaxRo = c.Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true};
}

static std::string FamilyName(AtomType t)
{ return t==AtomType::Gaussian ? "SG" : t==AtomType::Slater ? "SL" : "BSpline"; }

static std::string CaseName(const testing::TestParamInfo<PinCase>& i)
{ return FamilyName(i.param.type) + "_Z" + std::to_string(i.param.Z); }

} // anonymous namespace

class A_HF_dfPin : public ::testing::TestWithParam<PinCase> {};
TEST_P(A_HF_dfPin, SelfPin)
{
    const PinCase c = GetParam();
    AtomCalculation calc(c.Z, 0, {.type = c.type, .accuracy = c.acc, .model = Model::HF, .pol = Pol::Polarized},
                         HFParams(c));
    const double E = calc.Energy();

    // Harvest print (full precision) -- this is how the pin below was captured; also the ppm-vs-Saito line.
    std::cout.precision(12);
    std::cout << "[df-pin] " << CaseName({c, 0}) << "  E = " << E << "  ";
    RelativeHFError(E, c.Z);          // prints the signed rel. error vs the Saito table (not asserted here)

    EXPECT_TRUE(calc.IsConverged());
    if (c.pin != 0.0) EXPECT_NEAR(E, c.pin, c.atol);   // the tight regression guard (skipped while harvesting)
}

// --- Open-shell d: Sc (Z=21), 3d^1 4s^2, the clean single-unpaired-d ground state --------------------
// Self-pins harvested at High accuracy on this build.  atol 1e-8: ~4 orders below the family-completeness
// spread (BSpline-vs-SL ~8e-5, SG ~1.3e-3) and far below Saito, so it catches code drift while staying
// immune to it.  A legit numeric change re-harvests all three together (and the invariance test below
// confirms they still move as one).
constexpr double E_BSpline_Sc = -759.735911783995;   //  -0.26 ppm vs Saito
constexpr double E_SG_Sc      = -759.734619966359;   //  +1.45 ppm (Gaussian least complete)
constexpr double E_SL_Sc      = -759.735835752635;   //  -0.15 ppm
INSTANTIATE_TEST_SUITE_P(d_Sc, A_HF_dfPin, ::testing::ValuesIn(std::vector<PinCase>{
    {AtomType::BSpline6, High, 21, E_BSpline_Sc, 1e-8},
    {AtomType::Gaussian, High, 21, E_SG_Sc,      1e-8},
    {AtomType::Slater,   High, 21, E_SL_Sc,      1e-8},
}), CaseName);

// Cross-family invariance: HF energy is physical, so the three pins must agree to the completeness bound.
// A compute-free guard on the pins themselves -- if one family's basis data is mis-edited and its pin
// re-harvested, it desyncs from the other two here (whereas a genuine all-family code shift keeps them
// agreeing while the per-family SelfPin asserts flag "update the pins").  SG is looser (Gaussian truncation).
TEST(A_HF_dfPin, ScCrossFamilyInvariance)
{
    EXPECT_NEAR(E_SL_Sc, E_BSpline_Sc, 2e-4);   // the two near-complete families agree tightly
    EXPECT_NEAR(E_SG_Sc, E_BSpline_Sc, 2e-3);   // Gaussian looser (basis incompleteness, not code)
}

// --- Open-shell f: slot for a clean half-filled-f case (Eu Z=63 4f^7, or the U Z=92 already in A_HF_P) --
// Left for the owner to bless the pin: the element (term/config) choice is yours, and these are the
// hand-won numbers.  Wire it exactly like d_Sc once chosen.
