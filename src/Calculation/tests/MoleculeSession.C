// File: Calculation/tests/MoleculeSession.C  The MOLECULAR potential energy surface (doc/CalcNetworkReview.md §5 step 8, §7):
// H2O, HF/DZVP.  (a) ConvergeAt the reference == the single point, bit for bit; (b) a step from the carried density == a
// fresh run at the new positions, in fewer iterations, each step's first iterate near its answer; and the refusals that are outcomes.
#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>
import qchem.MoleculeCalculation;
import qchem.StructureData;
import qchem.Structure;
import qchem.SCFParams;
import qchem.SCFIterator;
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

// THE STEP'S CARRY starts the step near its answer (lane A, 2026-11-09): every step's FIRST iterate is within 1e-5 Ha of the step's
// converged energy (measured 1.7e-6 for the first step with the AO carry, CarryToMovedAtoms; the N-exact carries measured 1.5e-6 to
// 5.5e-6 here and ten times worse under LDA -- the table is at CarryToMovedAtoms).  DISCRIMINATING: with the carry switched off (each
// point seeded afresh) the first iterate is 4.1 Ha out and the gate fails.  Tr(D S') is a diagnostic: N to first order in the move.
TEST(MoleculeSession, EveryStepsFirstIterateStartsNearItsAnswer)
{
    const Molecule h2o=GetMolecule("H2O");
    std::vector<double> E;
    MoleculeSession pes(MoleculeSetup::Resolve(h2o, {.basis="dzvp"}), Tight(),
                        [&](const SCFIterator::IterationRow& r){ E.push_back(r.eb.GetTotalEnergy()); });
    std::vector<rvec3_t> R=pes.Reference().R;
    auto P0=pes.ConvergeAt(std::nullopt, R); ASSERT_TRUE(P0) << P0.Error().details;
    EXPECT_FALSE(pes.Current()->RestartedFrom()) << "the first point starts from the setup's seed, not a carry";
    const rvec3_t steps[2]={rvec3_t(0.0, 0.05, -0.04), rvec3_t(0.03, -0.02, 0.05)};
    for (const rvec3_t& d : steps)
    {
        R[1]=R[1]+d;   // one hydrogen
        E.clear();
        auto P=pes.ConvergeAt(std::nullopt, R); ASSERT_TRUE(P) << P.Error().details;
        ASSERT_FALSE(E.empty());
        EXPECT_LT(std::abs(E.front()-P->energy), 1e-5) << "the first iterate of a carried step";
        const auto& how=pes.Current()->RestartedFrom();
        ASSERT_TRUE(how) << "a step records its carry";
        EXPECT_FALSE(how->exact);
        EXPECT_NEAR(how->maxDisplacement, norm(d), 1e-12);
        EXPECT_NEAR(how->electrons, 10.0, 1e-2) << "Tr(D S'), a diagnostic: N to first order in the move";
    }
}

TEST(MoleculeSession, RefusalsAreOutcomesAndAWrongCountThrows)
{
    MoleculeSession pes(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .symmetry=true}), Tight());
    const std::vector<rvec3_t> R=pes.Reference().R;
    auto F=pes.ConvergeAt(std::nullopt, R, {.forces=true});
    ASSERT_FALSE(F); EXPECT_EQ(F.Error().why, PESFailure::Why::NoGradient);

    std::vector<rvec3_t> broken=R; broken[1].x+=0.05; broken[1].y-=0.03;   // one H only: C2v -> Cs
    EXPECT_THROW((void)pes.ConvergeAt(std::nullopt, broken), std::invalid_argument) << "an off-group R on an imposed setup is a wiring defect";

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
