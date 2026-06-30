// File A_HF_U.C  Atom Hartree-Fock tests for Unpolarized (closed shell) atoms (facade-driven).
//
// Migrated off the QchemTester/TestAtom scaffold onto qchem::AtomCalculation (OpenWork E).
//
// ONE parameterized fixture for the whole file: the case (basis family, accuracy, Z, tolerance) is the
// TEST PARAMETER, so each (basis,accuracy) variant is just another INSTANTIATE_TEST_SUITE_P -- the prefix
// (e.g. Slater_Medium) is the meaningful group name, and there are no empty per-variant fixture classes.
// The per-group tuned SCFParams live in one HFParams() switch; the assertion tolerance stays visible at
// each case.  Anchors/params byte-identical to the pre-refactor per-class tests.
#include "gtest/gtest.h"
#include <string>
import qchem.AtomCalculation;        // AtomCalculation, AtomType, BasisSetAccuracy, Model, Pol
import qchem.SCFIterator;            // SCFParams
import qchem.Unittests.TestUtils;    // RelativeHFError
using namespace qchem;
using enum BasisSetAccuracy;         // High, Medium, Low

//! One closed-shell HF atom test case: basis family + accuracy + element + the NIST relative-error bound.
//! lowerBound adds the BSpline "not too far below" guard (EXPECT_GT(error, -1e-4)).
struct HFCase { AtomType type; BasisSetAccuracy acc; int Z; double tol; bool lowerBound = false; };

// The per-(family,accuracy) tuned SCF params, Z-scaled.  Centralised so the cases stay one-liners; equal to
// the original per-class SCFParams literals.  (BSpline6 and BSpliner6 share params.)
static SCFParams HFParams(const HFCase& c)
{
    auto P = [&](size_t nmax, double dro, double dfd, double vir, double fd) {
        return SCFParams{.NMaxIter = nmax, .MinΔρ = c.Z*dro, .MinΔFD = dfd, .MinVirial = vir,
                         .MinFD = c.Z*fd, .StartingRelaxRo = c.Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true};
    };
    const bool bspline = (c.type==AtomType::BSpline6 || c.type==AtomType::BSpliner6);
    if (bspline) switch (c.acc) {
        case High:   return P(50, 1e-7, 1e-7, 2.5e-13, 1e-7);
        case Medium: return P(30, 1e-7, 1e-7, 2.5e-7,  1e-7);
        default:     return P(30, 1e-7, 1e-7, 5e-5,    1e-7);   // Low
    }
    if (c.type==AtomType::Gaussian) switch (c.acc) {
        case High:   return P(50, 1e-5, 1e-7, 1e-5, 1e-6);
        default:     return P(50, 1e-5, 1e-7, 5e-2, 1e-6);      // Medium
    }
    switch (c.acc) {   // Slater
        case High:   return P(32, 1e-5, 1e-7, 1e-6, 1e-6);
        case Medium: return P(22, 1e-4, 1e-5, 5e-4, 1e-6);
        default:     return P(30, 1e-4, 1e-4, 5e-1, 2e-5);      // Low
    }
}

// Name each case by element (so test ids read .../Z2, .../Z88 instead of .../0, .../1).
static std::string CaseName(const testing::TestParamInfo<HFCase>& i) { return "Z" + std::to_string(i.param.Z); }

class A_HF_U : public ::testing::TestWithParam<HFCase> {};
TEST_P(A_HF_U, Energy)
{
    const HFCase c = GetParam();
    AtomCalculation calc(c.Z, 0, {.type = c.type, .accuracy = c.acc, .model = Model::HF, .pol = Pol::UnPolarized},
                         HFParams(c));
    const double e = RelativeHFError(calc.Energy(), c.Z);
    EXPECT_LT(e, c.tol);
    if (c.lowerBound) EXPECT_GT(e, -1e-4);
    EXPECT_TRUE(calc.IsConverged());
}

#ifdef DEBUG
#define LOW
#else
#define MEDIUM
#define LOW
#endif

#ifdef HIGH
INSTANTIATE_TEST_SUITE_P(BSpline_High,  A_HF_U, ::testing::Values(HFCase{AtomType::BSpline6, High,2,1e-9,true}, HFCase{AtomType::BSpline6, High,88,1e-9,true}), CaseName);
INSTANTIATE_TEST_SUITE_P(BSpliner_High, A_HF_U, ::testing::Values(HFCase{AtomType::BSpliner6,High,2,1e-9,true}, HFCase{AtomType::BSpliner6,High,88,1e-9,true}), CaseName);
INSTANTIATE_TEST_SUITE_P(Gaussian_High, A_HF_U, ::testing::Values(HFCase{AtomType::Gaussian, High,2,2e-6},      HFCase{AtomType::Gaussian, High,36,2e-6}),      CaseName);
INSTANTIATE_TEST_SUITE_P(Slater_High,   A_HF_U, ::testing::Values(HFCase{AtomType::Slater,   High,2,1e-6},      HFCase{AtomType::Slater,   High,88,1e-6}),      CaseName);
#endif

#ifdef MEDIUM
INSTANTIATE_TEST_SUITE_P(BSpline_Medium,  A_HF_U, ::testing::Values(HFCase{AtomType::BSpline6, Medium,2,1e-6,true}, HFCase{AtomType::BSpline6, Medium,4,1e-6,true}), CaseName);
INSTANTIATE_TEST_SUITE_P(BSpliner_Medium, A_HF_U, ::testing::Values(HFCase{AtomType::BSpliner6,Medium,2,1e-6,true}, HFCase{AtomType::BSpliner6,Medium,4,1e-6,true}), CaseName);
INSTANTIATE_TEST_SUITE_P(Gaussian_Medium, A_HF_U, ::testing::Values(HFCase{AtomType::Gaussian, Medium,2,2e-4},      HFCase{AtomType::Gaussian, Medium,4,2e-4}),      CaseName);
INSTANTIATE_TEST_SUITE_P(Slater_Medium,   A_HF_U, ::testing::Values(HFCase{AtomType::Slater,   Medium,2,20e-6},     HFCase{AtomType::Slater,   Medium,88,20e-6}),    CaseName);
#endif

#ifdef LOW
#ifdef DEBUG
INSTANTIATE_TEST_SUITE_P(BSpline_Low,  A_HF_U, ::testing::Values(HFCase{AtomType::BSpline6, Low,2,40e-6,true}, HFCase{AtomType::BSpline6, Low,4,40e-6,true}), CaseName);
INSTANTIATE_TEST_SUITE_P(BSpliner_Low, A_HF_U, ::testing::Values(HFCase{AtomType::BSpliner6,Low,2,40e-6,true}, HFCase{AtomType::BSpliner6,Low,4,40e-6,true}), CaseName);
INSTANTIATE_TEST_SUITE_P(Slater_Low,   A_HF_U, ::testing::Values(HFCase{AtomType::Slater,   Low,2,0.01},       HFCase{AtomType::Slater,   Low,4,0.01}, HFCase{AtomType::Slater,Low,10,0.01}), CaseName);
#else
INSTANTIATE_TEST_SUITE_P(BSpline_Low,  A_HF_U, ::testing::Values(HFCase{AtomType::BSpline6, Low,2,40e-6,true}, HFCase{AtomType::BSpline6, Low,4,40e-6,true}), CaseName);
INSTANTIATE_TEST_SUITE_P(Slater_Low,   A_HF_U, ::testing::Values(HFCase{AtomType::Slater,   Low,2,0.01},       HFCase{AtomType::Slater,   Low,88,0.01}),     CaseName);
#endif
#endif
