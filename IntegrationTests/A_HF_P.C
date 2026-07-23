// File A_HF_P.C  Atom Hartree-Fock tests for Polarized (open shell) atoms (facade-driven).
//
// Migrated off the QchemTester/TestAtom scaffold onto qchem::AtomCalculation (OpenWork E).  ONE
// parameterized fixture: the case (basis family, accuracy, Z, tolerance, SCF knobs) is the test parameter,
// so each variant is just another INSTANTIATE_TEST_SUITE_P (prefix = group name) -- no empty fixture
// classes.  The per-group SCF knobs ride at the call site (visible, and factored across the Z-list by
// Cases()); they don't key cleanly on (type,accuracy) here, so they are data, not a switch.
#include "gtest/gtest.h"
#include <string>
#include <vector>
import qchem.AtomCalculation;        // AtomCalculation, AtomType, BasisSetAccuracy, Model, Pol
import qchem.SCFIterator;            // SCFParams
import qchem.Unittests.TestUtils;    // RelativeHFError
using namespace qchem;
using enum BasisSetAccuracy;         // High, Medium, Low

// File-local test helpers in an anonymous namespace (internal linkage): the per-file Knobs/HFCase
// vocabulary must not ODR-clash with the same names in sibling test TUs (e.g. A_HF_U's HFCase).
namespace {

//! Per-group SCF knobs (Z-scaled where noted in MakeParams).
struct Knobs { std::size_t nmax; double droScale, dfd, vir, fdScale; };
//! One open-shell HF atom case: basis family + accuracy + element + NIST relative-error bound + SCF knobs.
struct HFCase { AtomType type; BasisSetAccuracy acc; int Z; double tol; Knobs k; };

static SCFParams MakeParams(const HFCase& c)
{
    return {.NMaxIter = c.k.nmax, .MinΔρ = c.Z*c.k.droScale, .MinΔFD = c.k.dfd, .MinVirial = c.k.vir,
            .MinFD = c.Z*c.k.fdScale, .StartingRelaxRo = c.Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true};
}

// Expand a Z-list into cases sharing one basis family / accuracy / tolerance / SCF knobs.
static std::vector<HFCase> Cases(AtomType t, BasisSetAccuracy a, double tol, Knobs k, const std::vector<int>& Zs)
{
    std::vector<HFCase> v; v.reserve(Zs.size());
    for (int Z : Zs) v.push_back({t, a, Z, tol, k});
    return v;
}

static std::string CaseName(const testing::TestParamInfo<HFCase>& i) { return "Z" + std::to_string(i.param.Z); }

} // anonymous namespace

class A_HF_P : public ::testing::TestWithParam<HFCase> {};
TEST_P(A_HF_P, Energy)
{
    const HFCase c = GetParam();
    AtomCalculation calc(c.Z, 0, {.type = c.type, .accuracy = c.acc, .model = Model::HF, .pol = Pol::Polarized},
                         MakeParams(c));
    EXPECT_LT(RelativeHFError(calc.Energy(), c.Z), c.tol);
    EXPECT_TRUE(calc.IsConverged());
}

#ifdef DEBUG
#define LOW
#else
// #define HIGH
#define MEDIUM
#define LOW
#endif

#ifdef HIGH
INSTANTIATE_TEST_SUITE_P(Gaussian_High, A_HF_P, ::testing::ValuesIn(Cases(AtomType::Gaussian, High,2.5e-6,{50,1e-5,1e-7,1e-5,1e-6},   {5,21,92})), CaseName);
INSTANTIATE_TEST_SUITE_P(Slater_High,   A_HF_P, ::testing::ValuesIn(Cases(AtomType::Slater,   High,2.5e-6,{32,1e-5,1e-7,1e-6,1e-6},   {5,21,92})), CaseName);
INSTANTIATE_TEST_SUITE_P(BSpline_High,  A_HF_P, ::testing::ValuesIn(Cases(AtomType::BSpline6, High,2.5e-6,{50,1e-7,1e-7,2.5e-12,1e-7},{5,21})),    CaseName);
INSTANTIATE_TEST_SUITE_P(BSpliner_High, A_HF_P, ::testing::ValuesIn(Cases(AtomType::BSpliner6,High,2.5e-6,{50,1e-7,1e-6,5e-12,1e-7},  {5,21})),    CaseName);
#endif

#ifdef MEDIUM
INSTANTIATE_TEST_SUITE_P(Gaussian_Medium, A_HF_P, ::testing::ValuesIn(Cases(AtomType::Gaussian, Medium,60e-6,{50,1e-5,1e-7,5e-2,1e-6}, {3,5})),   CaseName);
INSTANTIATE_TEST_SUITE_P(Slater_Medium,   A_HF_P, ::testing::ValuesIn(Cases(AtomType::Slater,   Medium,20e-6,{22,1e-4,1e-5,5e-4,1e-6}, {3,5,21})), CaseName);
INSTANTIATE_TEST_SUITE_P(BSpline_Medium,  A_HF_P, ::testing::ValuesIn(Cases(AtomType::BSpline6, Medium,20e-6,{40,1e-7,1e-7,2.5e-7,1e-7},{3,5})),   CaseName);
INSTANTIATE_TEST_SUITE_P(BSpliner_Medium, A_HF_P, ::testing::ValuesIn(Cases(AtomType::BSpliner6,Medium,20e-6,{40,1e-7,1e-7,2.5e-7,1e-7},{3,5})),   CaseName);
#endif

#ifdef LOW
// "Low" grade: a faster (fewer-iter) Slater/Medium-accuracy spot check across a wide Z span.
#ifdef DEBUG
INSTANTIATE_TEST_SUITE_P(Slater_Low, A_HF_P, ::testing::ValuesIn(Cases(AtomType::Slater, Medium,12e-6,{14,1e-4,1e-4,1.1e-4,1e-5}, {3,5,21})),    CaseName);
#else
INSTANTIATE_TEST_SUITE_P(Slater_Low, A_HF_P, ::testing::ValuesIn(Cases(AtomType::Slater, Medium,12e-6,{14,1e-4,1e-4,1.1e-4,1e-5}, {3,5,21,92})), CaseName);
#endif
#endif
