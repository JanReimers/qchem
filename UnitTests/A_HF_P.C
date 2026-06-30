// File A_HF_P.C  Atom Hartree-Fock tests for Polarized (open shell) atoms (facade-driven).
//
// Migrated off the QchemTester/TestAtom scaffold onto qchem::AtomCalculation (retire QchemTester::Init).
// Same atomic bases, per-Z SCFParams, NIST oracle; anchors unchanged.  Polarized (unrestricted) HF.
#include "gtest/gtest.h"
import qchem.AtomCalculation;        // AtomCalculation, AtomType, BasisSetAccuracy, Model, Pol
import qchem.SCFIterator;            // SCFParams
import qchem.Unittests.TestUtils;    // RelativeHFError
using namespace qchem;

using std::cout;
using std::endl;
using enum BasisSetAccuracy;

const bool verbose=true;

// Run an open-shell (polarized) HF atom through the facade; report signed NIST relative error + convergence.
static double HF_P_Eerr(size_t Z, AtomType type, BasisSetAccuracy acc, const SCFParams& p, bool& converged)
{
    AtomCalculation calc(Z, 0, {.type = type, .accuracy = acc, .model = Model::HF, .pol = Pol::Polarized}, p);
    converged = calc.IsConverged();
    return RelativeHFError(calc.Energy(), int(Z));
}

#ifdef DEBUG
#define LOW
#else
// #define HIGH
#define MEDIUM
#define LOW
#endif

#ifdef HIGH
class A_SG_HF_P_High : public ::testing::TestWithParam<size_t> {};
TEST_P(A_SG_HF_P_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_P_Eerr(Z, AtomType::Gaussian, High,
        {.NMaxIter = 50, .MinΔρ = Z*1e-5, .MinΔFD = 1e-7, .MinVirial = 1e-5, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,2.5e-6);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,A_SG_HF_P_High,::testing::Values(5,21,92));

class A_SL_HF_P_High : public ::testing::TestWithParam<size_t> {};
TEST_P(A_SL_HF_P_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_P_Eerr(Z, AtomType::Slater, High,
        {.NMaxIter = 32, .MinΔρ = Z*1e-5, .MinΔFD = 1e-7, .MinVirial = 1e-6, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,2.5e-6);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,A_SL_HF_P_High,::testing::Values(5,21,92));

class A_BS_HF_P_High : public ::testing::TestWithParam<size_t> {};
TEST_P(A_BS_HF_P_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_P_Eerr(Z, AtomType::BSpline6, High,
        {.NMaxIter = 50, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 2.5e-12, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,2.5e-6);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,A_BS_HF_P_High,::testing::Values(5,21));

class A_BSr_HF_P_High : public ::testing::TestWithParam<size_t> {};
TEST_P(A_BSr_HF_P_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_P_Eerr(Z, AtomType::BSpliner6, High,
        {.NMaxIter = 50, .MinΔρ = Z*1e-7, .MinΔFD = 1e-6, .MinVirial = 5e-12, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,2.5e-6);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,A_BSr_HF_P_High,::testing::Values(5,21));
#endif //HIGH

#ifdef MEDIUM
class A_SG_HF_P_Medium : public ::testing::TestWithParam<size_t> {};
TEST_P(A_SG_HF_P_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_P_Eerr(Z, AtomType::Gaussian, Medium,
        {.NMaxIter = 50, .MinΔρ = Z*1e-5, .MinΔFD = 1e-7, .MinVirial = 5e-2, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,60e-6);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,A_SG_HF_P_Medium,::testing::Values(3,5));

class A_SL_HF_P_Medium : public ::testing::TestWithParam<size_t> {};
TEST_P(A_SL_HF_P_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_P_Eerr(Z, AtomType::Slater, Medium,
        {.NMaxIter = 22, .MinΔρ = Z*1e-4, .MinΔFD = 1e-5, .MinVirial = 5e-4, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,20e-6);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,A_SL_HF_P_Medium,::testing::Values(3,5,21));

class A_BS_HF_P_Medium : public ::testing::TestWithParam<size_t> {};
TEST_P(A_BS_HF_P_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_P_Eerr(Z, AtomType::BSpline6, Medium,
        {.NMaxIter = 40, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 2.5e-7, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,20e-6);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,A_BS_HF_P_Medium,::testing::Values(3,5));
class A_BSr_HF_P_Medium : public ::testing::TestWithParam<size_t> {};
TEST_P(A_BSr_HF_P_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_P_Eerr(Z, AtomType::BSpliner6, Medium,
        {.NMaxIter = 40, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 2.5e-7, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,20e-6);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,A_BSr_HF_P_Medium,::testing::Values(3,5));

#endif //MEDIUM

#ifdef LOW
class A_SL_HF_P_Low : public ::testing::TestWithParam<size_t> {};
TEST_P(A_SL_HF_P_Low,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_P_Eerr(Z, AtomType::Slater, Medium,
        {.NMaxIter = 14, .MinΔρ = Z*1e-4, .MinΔFD = 1e-4, .MinVirial = 1.1e-4, .MinFD = Z*1e-5, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,12e-6);
    EXPECT_TRUE(conv);
}
#ifdef DEBUG
INSTANTIATE_TEST_SUITE_P(A_HF,A_SL_HF_P_Low,::testing::Values(3,5,21));
#else
INSTANTIATE_TEST_SUITE_P(A_HF,A_SL_HF_P_Low,::testing::Values(3,5,21,92));
#endif

#endif
