// File: Calculation/tests/MoleculeSession.C  The MOLECULAR potential energy surface (doc/CalcNetworkReview.md §5 step 8, §7):
// H2O, HF/DZVP.  (a) ConvergeAt the reference == the single point, bit for bit; (b) a step from the carried density == a
// fresh run at the new positions, in fewer iterations; and the refusals that are outcomes.
#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>
#include <vector>
import qchem.Calculation;
import qchem.StructureData;
import qchem.Structure;
import qchem.SCFParams;
import qchem.Types;

using namespace qchem;
using StructureData::GetMolecule;

namespace
{
SCFParams Tight()
{
    SCFParams p; p.NMaxIter=100; p.MinΔρPerElectron=1e-9;
    return p;
}
std::vector<rvec3_t> Positions(const Structure& st)
{
    std::vector<rvec3_t> R;
    st.ForEachSite([&](int, const rvec3_t& r, bool){ R.push_back(r); });
    return R;
}
} // namespace

TEST(MoleculeSession, ConvergeAtTheReferenceIsTheSinglePointBitForBit)
{
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"});
    Calculation single(setup, Tight());
    ASSERT_TRUE(single.IsConverged());

    MoleculeSession pes(setup, Tight());
    auto P=pes.ConvergeAt(Positions(pes.Reference()));
    ASSERT_TRUE(P) << P.Error().details;
    EXPECT_DOUBLE_EQ(P->energy, single.Energy());
    EXPECT_EQ(P->iterations, single.IterationCount());
}

TEST(MoleculeSession, ACarriedStepEqualsAFreshRunAtTheNewPositionsInFewerIterations)
{
    const Molecule h2o=GetMolecule("H2O");
    MoleculeSession pes(MoleculeSetup::Resolve(h2o, {.basis="dzvp"}), Tight());
    std::vector<rvec3_t> R=Positions(pes.Reference());
    auto P0=pes.ConvergeAt(R); ASSERT_TRUE(P0) << P0.Error().details;
    R[1].x+=0.05; R[1].y-=0.03;   // bohr: one H off its equilibrium position
    auto P1=pes.ConvergeAt(R); ASSERT_TRUE(P1) << P1.Error().details;

    Molecule moved(h2o);
    for (size_t a=0;a<R.size();a++) moved[a]->itsR=R[a];
    Calculation fresh(moved, {.basis="dzvp"}, Tight());   // seeded from scratch at the new positions
    ASSERT_TRUE(fresh.IsConverged());
    EXPECT_NEAR(P1->energy, fresh.Energy(), 1e-8);
    EXPECT_LT(P1->iterations, fresh.IterationCount()) << "the carried density is the point of the session";
}

TEST(MoleculeSession, RefusalsAreOutcomesAndAWrongCountThrows)
{
    MoleculeSession pes(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .symmetry=true}), Tight());
    std::vector<rvec3_t> R=Positions(pes.Reference());
    auto F=pes.ConvergeAt(R, {.forces=true});
    ASSERT_FALSE(F); EXPECT_EQ(F.Error().why, PESFailure::Why::NoGradient);

    std::vector<rvec3_t> broken=R; broken[1].x+=0.05; broken[1].y-=0.03;   // one H only: C2v -> Cs
    auto B=pes.ConvergeAt(broken);
    ASSERT_FALSE(B); EXPECT_EQ(B.Error().why, PESFailure::Why::SymmetryBroken) << B.Error().details;

    EXPECT_THROW((void)pes.ConvergeAt({R[0]}), std::invalid_argument);
}
