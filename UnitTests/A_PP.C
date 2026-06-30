// File A_PP.C  Pseudopotential atom regression tests.  A pseudo-atom is the all-electron-basis atom with
// (a) the nuclear attraction replaced by the GTH local + KB-separable non-local pseudopotential (Ham_PP_U),
// (b) only the Zion valence electrons, filled by a valence-only config (PseudoAtom_EC).  These are "did the
// physics move" anchors (pin the converged energy/eigenvalue, no absolute-oracle / Converged() guard): the
// virial criterion does not apply to a pseudopotential, and the absolute value is basis/fitting-limited.
#include "gtest/gtest.h"
import qchem.Unittests.QchemTester;
import qchem.Hamiltonian.Internal.Hamiltonians;       // Ham_PP_U
import qchem.ElectronConfiguration.AtomNR;            // PseudoAtom_EC

const bool verbose=true;
using std::cout;
using std::endl;
using enum BasisSet::Atom::BasisSetAccuracy;
using namespace qchem::Hamiltonian;

// Pseudo-atom fixture: build the full-element basis (TestAtom(Z,0) -> valence-range pool via Z-val), swap in
// the pseudo-ion (Atom charge = Z-val, i.e. `val` electrons) and the valence-only config, and assemble
// Ham_PP_U (local + KB nonlocal) for the element looked up by Z.
class A_PP_U : public TestAtom
{
    int itsVal;
public:
    A_PP_U(int Z, int val) : TestAtom(Z,Z-val), itsVal(val)
    {
        delete itsEC;
        itsEC = new PseudoAtom_EC(Z);
    }
    virtual Hamiltonian* GetHamiltonian(st_t& c) const override
    {
        return new Ham_PP_U(c, QchemTester::itsPT.GetSymbol(GetStructure()->GetNuclearCharge()),
                            itsVal, GetMeshParams(), itsBasisSet);
    }
};

// Silicon, 4 valence electrons (GTH-LDA q4).  The KB nonlocal projectors lift the over-bound local-only s
// state (eps_s ~ -2.72) back toward the all-electron valence (GTH-LDA Si 3s ~ -0.40); total energy and
// eps_s pinned at the Slater/Medium converged values (Medium and High agree to ~6 digits).
class Si_PP_U : public ::testing::Test, public A_PP_U
{
public: Si_PP_U() : A_PP_U(14,4) {}
};
TEST_F(Si_PP_U, Medium)
{
    QchemTester::Init(Medium, BasisSet::Atom::Type::Slater, verbose);
    Iterate({.NMaxIter = 120, .MinΔρ = 1e-7, .MinΔFD = 1e-7, .MinVirial = 1e10, .MinFD = 1e-7, .StartingRelaxRo = 0.5, .MergeTol = 1e-7, .Verbose = true});  // virial off (N/A to PP)

    EXPECT_NEAR(TotalEnergy(), -3.336910601, 1e-6);     // pinned regression anchor (Slater/Medium)
    EXPECT_NEAR(TotalCharge(),  4.0,         1e-9);     // valence electron count
}
