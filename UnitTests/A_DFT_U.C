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
// Tolerance is the HONEST method accuracy of the atomic LDA-via-libxc path vs the NIST LDA reference: ~213
// ppm at Z=36.  That is genuine numerical/method error, NOT a functional error -- the functional is exact
// vs libxc to machine precision (LDA_XC_UT).  It is the same few-hundred-ppm ballpark as the in-house
// DiracVWN atom path (A_DFT_atom LSDA prints 331 ppm at Z=36), so 213 ppm is real: grid + Vxc fit, plus the
// lumped 3/4 exchange-virial XC energy this libxc path inherits from Ham_DFT_U (the DiracVWN path fixes that
// piece with a separate FittedVcorr; this path does not -- see the Libxc_LDA_Exchange follow-up note).  The
// old 2e-6 "match" was fictitious: a Z=36-tuned *1.006613 fudge on the exchange potential, now removed.  The
// 5e-4 bound reflects the real accuracy, not the fudge.
INSTANTIATE_TEST_SUITE_P(Gaussian_High, A_DFT_U, ::testing::ValuesIn(Cases(AtomType::Gaussian, 5e-4, {50,1e-5,1e-5,1e-1,1e-6}, {36})), CaseName);
INSTANTIATE_TEST_SUITE_P(Slater_High,   A_DFT_U, ::testing::ValuesIn(Cases(AtomType::Slater,   5e-4, {50,1e-5,1e-7,2e-2,1e-6}, {36})), CaseName);
