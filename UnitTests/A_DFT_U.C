// File A_DFT_U.C  Atom DFT tests using a libxc exchange functional (facade-driven).
//
// Migrated off QchemTester onto qchem::AtomCalculation with the new public exchange-functional selector:
// CalcOptions.xc = XCFunctional{.kind=XC::LibXC, .libxcId=7} routes the DFT Hamiltonian through libxc
// (functional id 7) instead of the built-in Slater-Xα.  Same High Gaussian/Slater bases, per-Z SCFParams,
// NIST oracle; anchors unchanged.
#include "gtest/gtest.h"
import qchem.AtomCalculation;        // AtomCalculation, AtomType, Model, Pol
import qchem.Hamiltonian.Factory;    // XCFunctional, XC (the exchange-functional selector)
import qchem.SCFIterator;            // SCFParams
import qchem.Unittests.TestUtils;    // RelativeDFTError
using namespace qchem;
using std::cout;
using std::endl;
using enum BasisSetAccuracy;
using namespace qchem::Hamiltonian;  // XCFunctional, XC

const bool verbose=true;

// libxc LDA (functional id 7) atom DFT; report signed NIST relative error + convergence.
static double Libxc_DFT_err(size_t Z, AtomType type, const SCFParams& p, bool& converged)
{
    AtomCalculation calc(Z, 0, {.type = type, .accuracy = High, .model = Model::Xalpha, .pol = Pol::UnPolarized,
                                .xc = XCFunctional{.kind = XC::LibXC, .libxcId = 7}}, p);
    converged = calc.IsConverged();
    return RelativeDFTError(calc.Energy(), int(Z));
}

class SG_DFT_U_High : public ::testing::TestWithParam<size_t> {};
TEST_P(SG_DFT_U_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double err=Libxc_DFT_err(Z, AtomType::Gaussian,
        {.NMaxIter = 50, .MinΔρ = Z*1e-5, .MinΔFD = 1e-5, .MinVirial = 1e-1, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(err,2e-6);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_DFT,SG_DFT_U_High,::testing::Values(36));

class SL_DFT_U_High : public ::testing::TestWithParam<size_t> {};
TEST_P(SL_DFT_U_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
    bool conv; double err=Libxc_DFT_err(Z, AtomType::Slater,
        {.NMaxIter = 50, .MinΔρ = Z*1e-5, .MinΔFD = 1e-7, .MinVirial = 2e-2, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true}, conv);
    EXPECT_LT(err,2e-6);
    EXPECT_TRUE(conv);
}
INSTANTIATE_TEST_SUITE_P(A_DFT,SL_DFT_U_High,::testing::Values(36));
