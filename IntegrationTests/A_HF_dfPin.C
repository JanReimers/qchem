// File A_HF_dfPin.C  Tight self-pinned REGRESSION anchors for open-shell d/f atoms, across the three
// non-relativistic atomic basis families: BSpline6, Gaussian (SG) and Slater (SL).
//
// DESIGN.  These pin the code's OWN converged energy (a "did E move" self-anchor, like the molecular M_*
// tests) at a tight 1e-8, so ANY drift in the ERI / SCF path trips them -- DECOUPLED from physical accuracy.
// That decoupling is the whole point: we run at Accuracy=Low (cheap, small basis), yet the tight pin still
// catches code changes.  So the always-on regression tier is fast, and the expensive Accuracy=High runs are
// opt-in physical validation under #ifdef HIGH (the A_HF_U/A_HF_P idiom).
//
//   Low  (default) : E == pinned self-value to 1e-8   -- fast drift guard, NOT an accuracy claim
//   High (#ifdef)  : |rel err vs Saito| < tol          -- expensive correctness validation, opt-in
//
// Elements: Sc(21) 3d^1 (single unpaired d), Eu(63) 4f^7 (half-filled f, max spin), U(92) 5f^3 6d^1.
#include "gtest/gtest.h"
#include <string>
#include <vector>
#include <iostream>
import qchem.AtomCalculation;        // AtomCalculation, AtomType, BasisSetAccuracy, Model, Pol
import qchem.SCFIterator;            // SCFParams
import qchem.Unittests.TestUtils;    // RelativeHFError (Saito oracle)
using namespace qchem;
using enum BasisSetAccuracy;         // High, Medium, Low

namespace {

//! One open-shell HF self-pin case: family + accuracy + element + pinned converged energy + abs tolerance.
//! pin==0.0 == "harvest mode": only prints E (no regression assert), so the pin can be read off a run.
struct PinCase { AtomType type; BasisSetAccuracy acc; int Z; double pin; double atol; };

//! Robust SCF knobs for the (heavy, open-shell f) atoms at Low grade.  Z-scaled thresholds, generous iter
//! budget + slow relaxation for the Z>40 f-shell atoms (near-degenerate, oscillation-prone).
static SCFParams HFParams(const PinCase& c)
{
    const bool heavy = c.Z >= 40;
    return {.NMaxIter = 120, .MinΔρ = c.Z*1e-4, .MinΔFD = 1e-4, .MinVirial = 4e-2,
            .MinFD = c.Z*1e-5, .StartingRelaxRo = heavy ? 0.2 : 0.5, .MergeTol = 1e-7, .Verbose = true};
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

    std::cout.precision(12);
    std::cout << "[df-pin] " << CaseName({c, 0}) << "  E = " << E << "  ";
    RelativeHFError(E, c.Z);                            // prints signed rel. error vs Saito (context only)

    EXPECT_TRUE(calc.IsConverged());
    if (c.pin != 0.0) EXPECT_NEAR(E, c.pin, c.atol);   // tight regression guard (skipped while harvesting)
}

// --- Low tier (always on): tight self-pins on the code's own converged Accuracy=Low energies -----------
// Values harvested on this build; atol 1e-8 (HF is deterministic, so a converged Low run reproduces
// bit-stably -- the crude basis makes these physically inaccurate, but that is irrelevant to a drift guard).
// U(92) on Slater is OMITTED: the Slater-Low basis is too crude to converge Z=92 (IsConverged() fails); the
// BSpline/SG anchors cover U.  Physical validation vs Saito is the opt-in Accuracy=High job (see note below).
INSTANTIATE_TEST_SUITE_P(d_Sc, A_HF_dfPin, ::testing::Values(
    PinCase{AtomType::BSpline6, Low, 21,   -759.732835713717, 1e-8},
    PinCase{AtomType::Gaussian, Low, 21,   -759.488546501283, 1e-8},
    PinCase{AtomType::Slater,   Low, 21,   -758.357843165290, 1e-8}), CaseName);
INSTANTIATE_TEST_SUITE_P(f_Eu, A_HF_dfPin, ::testing::Values(
    PinCase{AtomType::BSpline6, Low, 63, -10423.363067460416, 1e-8},
    PinCase{AtomType::Gaussian, Low, 63, -10422.598808202592, 1e-8},
    PinCase{AtomType::Slater,   Low, 63, -10391.314921936011, 1e-8}), CaseName);
INSTANTIATE_TEST_SUITE_P(f_U,  A_HF_dfPin, ::testing::Values(
    PinCase{AtomType::BSpline6, Low, 92, -25663.352928188484, 1e-8},
    PinCase{AtomType::Gaussian, Low, 92, -25661.494406485683, 1e-8},
    PinCase{AtomType::Slater  , Low, 92, -25465.307086748297, 1e-8}), CaseName);

// Opt-in physical validation (expensive): build with -DHIGH to re-run these elements at Accuracy=High and
// bound the error vs the Saito near-HF-limit table (the A_HF_U/A_HF_P idiom).  Kept OFF by default so the
// always-on tier stays fast; the tight Low self-pins above are what guard against code drift day to day.
#ifdef HIGH
class A_HF_dfPin_HighValidation : public ::testing::TestWithParam<PinCase> {};
TEST_P(A_HF_dfPin_HighValidation, VsSaito)
{
    const PinCase c = GetParam();
    AtomCalculation calc(c.Z, 0, {.type = c.type, .accuracy = High, .model = Model::HF, .pol = Pol::Polarized},
                         HFParams(c));
    EXPECT_TRUE(calc.IsConverged());
    EXPECT_LT(RelativeHFError(calc.Energy(), c.Z), c.atol);   // atol here is the Saito relative-error bound
}
INSTANTIATE_TEST_SUITE_P(dfHigh, A_HF_dfPin_HighValidation, ::testing::Values(
    PinCase{AtomType::BSpline6, High, 21, 0.0, 5e-6}, PinCase{AtomType::BSpline6, High, 63, 0.0, 5e-6},
    PinCase{AtomType::BSpline6, High, 92, 0.0, 5e-6}), CaseName);
#endif
