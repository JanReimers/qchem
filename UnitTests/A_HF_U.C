// File A_HF_U.C  Atom Hartree-Fock tests for Unpolarized (closed shell) atoms (facade-driven).
//
// Migrated off the QchemTester/TestAtom scaffold onto qchem::AtomCalculation -- the single-atom sibling
// of the molecular qchem::Calculation front door (OpenWork: retire QchemTester::Init).  Same atomic
// exponent-pool bases, same per-Z SCFParams, same NIST oracle (RelativeHFError now a free Z-keyed helper);
// anchors unchanged.
#include "gtest/gtest.h"
import qchem.AtomCalculation;        // AtomCalculation, AtomCalcOptions, AtomType, BasisSetAccuracy, Model
import qchem.SCFIterator;            // SCFParams
import qchem.Unittests.TestUtils;    // RelativeHFError
using namespace qchem;

using std::cout;
using std::endl;
using enum BasisSetAccuracy;

const bool verbose=true;

// Run a closed-shell HF atom through the facade; report the signed NIST relative error and convergence.
static double HF_Eerr(size_t Z, AtomType type, BasisSetAccuracy acc, const SCFParams& p, bool& converged)
{
    AtomCalculation calc(Z, 0, {.type = type, .accuracy = acc, .model = Model::HF, .pol = Pol::UnPolarized});
    calc.Converge(p);
    converged = calc.IsConverged();
    return RelativeHFError(calc.Energy(), int(Z));
}

#ifdef DEBUG
#define LOW
#else
#define MEDIUM
#define LOW
#endif


#ifdef HIGH
class BS_U_High : public ::testing::TestWithParam<size_t> {};
TEST_P(BS_U_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_Eerr(Z, AtomType::BSpline6, High,
        {.NMaxIter = 50, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 2.5e-13, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,1e-9);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,BS_U_High,::testing::Values(2,88));
// 2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88
class BSr_U_High : public ::testing::TestWithParam<size_t> {};
TEST_P(BSr_U_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_Eerr(Z, AtomType::BSpliner6, High,
        {.NMaxIter = 50, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 2.5e-13, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,1e-9);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,BSr_U_High,::testing::Values(2,88));
class SG_U_High : public ::testing::TestWithParam<size_t> {};
TEST_P(SG_U_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_Eerr(Z, AtomType::Gaussian, High,
        {.NMaxIter = 50, .MinΔρ = Z*1e-5, .MinΔFD = 1e-7, .MinVirial = 1e-5, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,2e-6);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,SG_U_High,::testing::Values(2,36));//));

class SL_U_High : public ::testing::TestWithParam<size_t> {};
TEST_P(SL_U_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_Eerr(Z, AtomType::Slater, High,
        {.NMaxIter = 32, .MinΔρ = Z*1e-5, .MinΔFD = 1e-7, .MinVirial = 1e-6, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,1e-6);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,SL_U_High,::testing::Values(2,88));//));

#endif //HIGH

#ifdef MEDIUM

class BS_U_Medium : public ::testing::TestWithParam<size_t> {};
TEST_P(BS_U_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_Eerr(Z, AtomType::BSpline6, Medium,
        {.NMaxIter = 30, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 2.5e-7, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,1e-6);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,BS_U_Medium,::testing::Values(2,4));//));
class BSr_U_Medium : public ::testing::TestWithParam<size_t> {};
TEST_P(BSr_U_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_Eerr(Z, AtomType::BSpliner6, Medium,
        {.NMaxIter = 30, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 2.5e-7, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,1e-6);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,BSr_U_Medium,::testing::Values(2,4));//));

class SG_U_Medium : public ::testing::TestWithParam<size_t> {};
TEST_P(SG_U_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_Eerr(Z, AtomType::Gaussian, Medium,
        {.NMaxIter = 50, .MinΔρ = Z*1e-5, .MinΔFD = 1e-7, .MinVirial = 5e-2, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,2e-4);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,SG_U_Medium,::testing::Values(2,4));//));

class SL_U_Medium : public ::testing::TestWithParam<size_t> {};
TEST_P(SL_U_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_Eerr(Z, AtomType::Slater, Medium,
        {.NMaxIter = 22, .MinΔρ = Z*1e-4, .MinΔFD = 1e-5, .MinVirial = 5e-4, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,20e-6);
    EXPECT_TRUE(conv);
}

INSTANTIATE_TEST_SUITE_P(A_HF,SL_U_Medium,::testing::Values(2,88));//));

#endif //MEDIUM

#ifdef LOW

class BS_U_Low : public ::testing::TestWithParam<size_t> {};
class BSr_U_Low : public ::testing::TestWithParam<size_t> {};
#ifdef DEBUG
TEST_P(BS_U_Low,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_Eerr(Z, AtomType::BSpline6, Low,
        {.NMaxIter = 30, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 5e-5, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,40e-6);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,BS_U_Low,::testing::Values(2,4));//));
TEST_P(BSr_U_Low,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_Eerr(Z, AtomType::BSpliner6, Low,
        {.NMaxIter = 30, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 5e-5, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,40e-6);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,BSr_U_Low,::testing::Values(2,4));//));
#else
TEST_P(BS_U_Low,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_Eerr(Z, AtomType::BSpline6, Low,
        {.NMaxIter = 30, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 5e-5, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,40e-6);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_HF,BS_U_Low,::testing::Values(2,4));//));
#endif

class SL_U_Low : public ::testing::TestWithParam<size_t> {};
TEST_P(SL_U_Low,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double Eerr=HF_Eerr(Z, AtomType::Slater, Low,
        {.NMaxIter = 30, .MinΔρ = Z*1e-4, .MinΔFD = 1e-4, .MinVirial = 5e-1, .MinFD = Z*2e-5, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(Eerr,0.01); //1%
    EXPECT_TRUE(conv);
}
#ifdef DEBUG
INSTANTIATE_TEST_SUITE_P(A_HF,SL_U_Low,::testing::Values(2,4,10));//));
#else
INSTANTIATE_TEST_SUITE_P(A_HF,SL_U_Low,::testing::Values(2,88));//));
#endif

#endif //LOW
