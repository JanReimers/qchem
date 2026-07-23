// File A_HF_U.C  Atom Hartree-Fock tests for Unpolarized (closed shell) atoms (facade-driven).
//
// Migrated off the QchemTester/TestAtom scaffold onto qchem::AtomCalculation (OpenWork E).
//
// ONE parameterized fixture for the whole file: the case (basis family, accuracy, Z, tolerance) is the
// TEST PARAMETER, so each (basis,accuracy) variant is just another INSTANTIATE_TEST_SUITE_P -- the prefix
// (e.g. Slater_Medium) is the group name and there are no empty per-variant fixture classes.  A Z-list is a
// plain braced list fed through Cases()/ValuesIn(), so sweeping many elements is one line.  Per-group tuned
// SCFParams live in one HFParams() switch; the assertion tolerance stays visible at each case.
#include "gtest/gtest.h"
#include <string>
#include <vector>
import qchem.AtomCalculation;        // AtomCalculation, AtomType, BasisSetAccuracy, Model, Pol
import qchem.SCFIterator;            // SCFParams
import qchem.Unittests.TestUtils;    // RelativeHFError
using namespace qchem;
using enum BasisSetAccuracy;         // High, Medium, Low

// File-local test helpers in an anonymous namespace: HFCase/the struct vocabulary are intentionally
// per-file (each A_* test tunes its own), so internal linkage avoids ODR clashes with the same names in
// sibling test TUs (e.g. A_HF_P's differently-shaped HFCase).
namespace {

//! One closed-shell HF atom test case: basis family + accuracy + element + the NIST relative-error bound.
//! lowerBound adds the BSpline "not too far below" guard (EXPECT_GT(error, -1e-4)).
struct HFCase { AtomType type; BasisSetAccuracy acc; int Z; double tol; bool lowerBound = false; };

// Closed-(sub)shell "full l-shell" atoms: the well-defined spherical UNpolarized ground states -- He, Be,
// Ne, Mg, Ar, Ca, Zn, Kr, Sr, Pd, Cd, Xe, Ba, Hg, Rn, Ra.  The de-noised version of the old commented
// Z-lists; converges across the periodic table at the standard Slater/Medium grade.  (Yb(70) 4f¹⁴ and the
// actinides/No(102)/Og(118) need a finer basis grade -- add them to a High sweep, not here.)
static const std::vector<int> ClosedShell = {2,4,10,12,18,20,30,36,38,46,48,54,56,80,86,88};

// Expand a Z-list into cases sharing one basis family / accuracy / tolerance / lower-bound flag.
static std::vector<HFCase> Cases(AtomType t, BasisSetAccuracy a, double tol, bool lb, const std::vector<int>& Zs)
{
    std::vector<HFCase> v; v.reserve(Zs.size());
    for (int Z : Zs) v.push_back({t, a, Z, tol, lb});
    return v;
}

// The per-(family,accuracy) tuned SCF params, Z-scaled.  Centralised so the cases stay data; equal to the
// original per-class SCFParams literals.  (BSpline6 and BSpliner6 share params.)
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

} // anonymous namespace

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
INSTANTIATE_TEST_SUITE_P(BSpline_High,  A_HF_U, ::testing::ValuesIn(Cases(AtomType::BSpline6, High,1e-9,true,{2,88})), CaseName);
INSTANTIATE_TEST_SUITE_P(BSpliner_High, A_HF_U, ::testing::ValuesIn(Cases(AtomType::BSpliner6,High,1e-9,true,{2,88})), CaseName);
INSTANTIATE_TEST_SUITE_P(Gaussian_High, A_HF_U, ::testing::ValuesIn(Cases(AtomType::Gaussian, High,2e-6,false,{2,36})), CaseName);
INSTANTIATE_TEST_SUITE_P(Slater_High,   A_HF_U, ::testing::ValuesIn(Cases(AtomType::Slater,   High,1e-6,false,{2,88})), CaseName);
#endif

#ifdef MEDIUM
// Slater/Medium is the cheap, robust grade.  Regular runs spot-check light/mid/heavy; uncomment the
// ClosedShell line for the full periodic sweep (fun, but ~10 s).
INSTANTIATE_TEST_SUITE_P(Slater_Medium,   A_HF_U, ::testing::ValuesIn(Cases(AtomType::Slater,   Medium,20e-6,false,{2,18,88})), CaseName);
// INSTANTIATE_TEST_SUITE_P(Slater_sweep,  A_HF_U, ::testing::ValuesIn(Cases(AtomType::Slater,   Medium,20e-6,false,ClosedShell)), CaseName);
INSTANTIATE_TEST_SUITE_P(Gaussian_Medium, A_HF_U, ::testing::ValuesIn(Cases(AtomType::Gaussian, Medium,2e-4,false,{2,4})), CaseName);
INSTANTIATE_TEST_SUITE_P(BSpline_Medium,  A_HF_U, ::testing::ValuesIn(Cases(AtomType::BSpline6, Medium,1e-6,true,{2,4})), CaseName);
INSTANTIATE_TEST_SUITE_P(BSpliner_Medium, A_HF_U, ::testing::ValuesIn(Cases(AtomType::BSpliner6,Medium,1e-6,true,{2,4})), CaseName);
#endif

#ifdef LOW
INSTANTIATE_TEST_SUITE_P(BSpline_Low, A_HF_U, ::testing::ValuesIn(Cases(AtomType::BSpline6, Low,40e-6,true,{2,4})), CaseName);
#ifdef DEBUG
INSTANTIATE_TEST_SUITE_P(Slater_Low,  A_HF_U, ::testing::ValuesIn(Cases(AtomType::Slater, Low,0.01,false,{2,4,10})), CaseName);
#else
INSTANTIATE_TEST_SUITE_P(Slater_Low,  A_HF_U, ::testing::ValuesIn(Cases(AtomType::Slater, Low,0.01,false,{2,88})), CaseName);
#endif
#endif
