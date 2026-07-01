// File A_DFT_U.C  Atom DFT tests using a libxc exchange functional (facade-driven).
//
// One parameterized fixture: CalcOptions.xc = XCFunctional{.kind=XC::LibXC, .libxcId=7} routes the DFT
// Hamiltonian through libxc (functional id 7).  Each (basis) variant is an INSTANTIATE_TEST_SUITE_P; the
// per-group SCF knobs ride at the call site.  Helpers anon-namespaced (per-file vocabulary, no ODR clash).
#include "gtest/gtest.h"
#include <string>
#include <vector>
import qchem.AtomCalculation;        // AtomCalculation, AtomType, Model, Pol
import qchem.Hamiltonian.Factory;    // XCFunctional, XC (the exchange-functional selector)
import qchem.SCFIterator;            // SCFParams
import qchem.Unittests.TestUtils;    // RelativeDFTError
using namespace qchem;
using enum BasisSetAccuracy;         // High
using namespace qchem::Hamiltonian;  // XCFunctional, XC

namespace {

struct Knobs { std::size_t nmax; double droScale, dfd, vir, fdScale; };
struct DFTCase { AtomType type; int Z; double tol; Knobs k; };

static SCFParams MakeParams(const DFTCase& c)
{
    return {.NMaxIter = c.k.nmax, .MinΔρ = c.Z*c.k.droScale, .MinΔFD = c.k.dfd, .MinVirial = c.k.vir,
            .MinFD = c.Z*c.k.fdScale, .StartingRelaxRo = c.Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true};
}
static std::vector<DFTCase> Cases(AtomType t, double tol, Knobs k, const std::vector<int>& Zs)
{
    std::vector<DFTCase> v; v.reserve(Zs.size());
    for (int Z : Zs) v.push_back({t, Z, tol, k});
    return v;
}
static std::string CaseName(const testing::TestParamInfo<DFTCase>& i) { return "Z" + std::to_string(i.param.Z); }

} // anonymous namespace

class A_DFT_U : public ::testing::TestWithParam<DFTCase> {};
TEST_P(A_DFT_U, Libxc)
{
    const DFTCase c = GetParam();
    AtomCalculation calc(c.Z, 0, {.type = c.type, .accuracy = High, .model = Model::Xalpha, .pol = Pol::UnPolarized,
                                  .xc = XCFunctional{.kind = XC::LibXC, .libxcId = 7}}, MakeParams(c));
    EXPECT_LT(RelativeDFTError(calc.Energy(), c.Z), c.tol);
    EXPECT_TRUE(calc.IsConverged());
}

// FittedVcorr cleanliness: the libxc SVWN5 path (Dirac-X + libxc-VWN-C) and the in-house DiracVWN path
// (SlaterExchange + VWN_Correlation) are the SAME functional with the SAME correct correlation-energy
// treatment (both FittedVxc + FittedVcorr now), so on one atom/basis they must agree to the ~1e-9 functional
// implementation gap (LDA_XC_UT), NOT the ~100 ppm the old lumped-energy libxc path differed by.  Kr/Slater.
TEST(A_DFT_U_Consistency, LibxcSVWN5MatchesInHouseDiracVWN)
{
    const DFTCase c{AtomType::Slater, 36, 0.0, {50,1e-5,1e-7,2e-2,1e-6}};
    AtomCalculation libxc(c.Z, 0, {.type = c.type, .accuracy = High, .model = Model::Xalpha,
                                   .xc = XCFunctional{.kind = XC::LibXC, .libxcId = 7}}, MakeParams(c));
    AtomCalculation dirac(c.Z, 0, {.type = c.type, .accuracy = High, .model = Model::LDA}, MakeParams(c));
    EXPECT_NEAR(libxc.Energy(), dirac.Energy(), 1e-4);   // same SVWN5 + same correct E_c => agree (Kr ~ -2750)
}
// Tolerance is the HONEST accuracy of the atomic libxc SVWN5 path vs the NIST LDA reference: ~1.5 ppm
// (Gaussian) / ~0.7 ppm (Slater) at Z=36, with NO fudge.  History: the old `2e-6` match was faked by a
// Z=36-tuned *1.006613 exchange fudge (removed); removing it alone left ~213 ppm -- which turned out to be
// NOT grid/fit but the lumped 3/4 exchange-virial CORRELATION energy (the libxc path used to combine X+C in
// one FittedVxc).  Splitting it into FittedVxc(Dirac-X) + FittedVcorr(libxc-C), so E_c = integral eps_c rho,
// dropped the error ~140x to ~1 ppm -- honestly, and it now generalises off Z=36 (the fudge did not).  The
// functional itself is exact vs libxc (LDA_XC_UT); the residual is the true grid/Vxc-fit floor at High basis.
INSTANTIATE_TEST_SUITE_P(Gaussian_High, A_DFT_U, ::testing::ValuesIn(Cases(AtomType::Gaussian, 5e-6, {50,1e-5,1e-5,1e-1,1e-6}, {36})), CaseName);
INSTANTIATE_TEST_SUITE_P(Slater_High,   A_DFT_U, ::testing::ValuesIn(Cases(AtomType::Slater,   5e-6, {50,1e-5,1e-7,2e-2,1e-6}, {36})), CaseName);
