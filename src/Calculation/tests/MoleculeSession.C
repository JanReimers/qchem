// File: Calculation/tests/MoleculeSession.C  The MOLECULAR potential energy surface (doc/CalcNetworkReview.md §5 step 8, §7):
// H2O, HF/DZVP.  (a) ConvergeAt the reference == the single point, bit for bit; (b) a step from the carried density == a
// fresh run at the new positions, in fewer iterations; and the refusals that are outcomes.
#include <gtest/gtest.h>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>
import qchem.MoleculeCalculation;
import qchem.StructureData;
import qchem.Structure;
import qchem.SCFParams;
import qchem.Types;
import qchem.SCFCap;   // BoundedSCF: the SCF-failure arm

using namespace qchem;
using StructureData::GetMolecule;

namespace
{
SCFParams Tight()
{
    SCFParams p; p.NMaxIter=100; p.MinΔρPerElectron=1e-9;
    return p;
}
} // namespace

TEST(MoleculeSession, ConvergeAtTheReferenceIsTheSinglePointBitForBit)
{
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"});
    MoleculeCalculation single(setup, Tight());
    ASSERT_TRUE(single.IsConverged());

    MoleculeSession pes(setup, Tight());
    const PESGeometry ref=pes.Reference();
    EXPECT_FALSE(ref.cell) << "a molecule has no cell";
    auto P=pes.ConvergeAt(std::nullopt, ref.R);
    ASSERT_TRUE(P) << P.Error().details;
    EXPECT_DOUBLE_EQ(P->energy, single.Energy());
    EXPECT_EQ(P->iterations, single.IterationCount());
}

TEST(MoleculeSession, ACarriedStepEqualsAFreshRunAtTheNewPositionsInFewerIterations)
{
    const Molecule h2o=GetMolecule("H2O");
    MoleculeSession pes(MoleculeSetup::Resolve(h2o, {.basis="dzvp"}), Tight());
    std::vector<rvec3_t> R=pes.Reference().R;
    auto P0=pes.ConvergeAt(std::nullopt, R); ASSERT_TRUE(P0) << P0.Error().details;
    R[1].x+=0.05; R[1].y-=0.03;   // bohr, Cartesian: one H off its equilibrium position
    auto P1=pes.ConvergeAt(std::nullopt, R); ASSERT_TRUE(P1) << P1.Error().details;

    Molecule moved(h2o);
    moved.MoveSites(R);
    MoleculeCalculation fresh(MoleculeSetup::Resolve(moved, {.basis="dzvp"}), Tight());   // seeded from scratch at the new positions
    ASSERT_TRUE(fresh.IsConverged());
    EXPECT_NEAR(P1->energy, fresh.Energy(), 1e-8);
    EXPECT_LT(P1->iterations, fresh.IterationCount()) << "the carried density is the point of the session";
}

TEST(MoleculeSession, RefusalsAreOutcomesAndAWrongCountThrows)
{
    MoleculeSession pes(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .symmetry=true}), Tight());
    const std::vector<rvec3_t> R=pes.Reference().R;
    auto F=pes.ConvergeAt(std::nullopt, R, {.forces=true});
    ASSERT_FALSE(F); EXPECT_EQ(F.Error().why, PESFailure::Why::NoGradient);

    std::vector<rvec3_t> broken=R; broken[1].x+=0.05; broken[1].y-=0.03;   // one H only: C2v -> Cs
    auto B=pes.ConvergeAt(std::nullopt, broken);
    ASSERT_FALSE(B); EXPECT_EQ(B.Error().why, PESFailure::Why::SymmetryBroken) << B.Error().details;

    EXPECT_THROW((void)pes.ConvergeAt(std::nullopt, {R[0]}), std::invalid_argument);
    rmat3d_t box; for (size_t i=1;i<=3;i++) for (size_t j=1;j<=3;j++) box(i,j)= i==j ? 20.0 : 0.0;
    EXPECT_THROW((void)pes.ConvergeAt(box, R), std::invalid_argument) << "a molecule has no cell";
}

// A point whose SCF fails is Why::SCF CARRYING the SCF's own failure (the shared SCFFailure, 2026-10-31), so a driver can see
// why and retry with another recipe.
TEST(MoleculeSession, AnSCFFailureIsAnOutcomeCarryingTheSCFsOwnReason)
{
    qchem::BoundedSCF bounded("the failure arm is the subject: one iteration cannot converge water");
    MoleculeSession pes(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"}), {.NMaxIter=1});
    auto F=pes.ConvergeAt(std::nullopt, pes.Reference().R);
    ASSERT_FALSE(F);
    EXPECT_EQ(F.Error().why, PESFailure::Why::SCF);
    ASSERT_TRUE(F.Error().scf.has_value());
    EXPECT_EQ(F.Error().scf->why, SCFFailure::Why::NotConverged);
    EXPECT_EQ(F.Error().scf->iterations, 1u);
    EXPECT_EQ(F.Error().details, F.Error().scf->details);
}
