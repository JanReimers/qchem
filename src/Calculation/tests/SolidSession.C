// File: Calculation/tests/SolidSession.C  The SOLID potential energy surface (doc/CalcNetworkReview.md §5 step 8, §7): a
// session over one setup, sampled at positions the caller chooses.  Si Γ (two atoms, eight electrons, ~1 s an SCF) carries
// the gates: (a) ConvergeAt the reference == the single point, bit for bit; (b) a step from the carried density == a fresh run
// at the new positions, in fewer iterations; and the refusals that are outcomes, not crashes.
#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>
import qchem.SolidCalculation;
import qchem.Materials;
import qchem.Lattice_3D;
import qchem.Structure;
import qchem.UnitCell;
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.SCFParams;
import qchem.SCFAccelerator.Factory;   // SCFAccelerators::Type
import qchem.Types;

using namespace qchem;

namespace
{
SCFParams Tight()
{
    SCFParams p;
    // Δρ/N 1e-8: above Si Γ's Δρ quadrature floor (ρ_lost/N ~1.8e-9 at densityEcut 20), where a 1e-9 gate wanders.
    p.NMaxIter=100; p.MinΔρPerElectron=1e-8; p.MinΔE=1e30; p.MinΔFD=1e30; p.MinVirial=1e30; p.MinFD=1e30;
    return p;
}
SolidCalcOptions SiGamma(bool imposed=false)
{
    SolidCalcOptions o; o.Nelec=8; o.densityEcut=20.0; o.imposeSymmetry=imposed; o.label="pes";
    return o;
}
const GaussianBasisRecipe kSiBasis{.data=BasisSet::Gaussian::BasisSetData::SIPP_SR};

std::shared_ptr<const SolidSetup> SiSetup(const std::vector<rvec3_t>* at=nullptr, bool imposed=false)
{
    const auto mat=Materials::Get("Si_diamond");
    UnitCell cell(*mat.cell);
    if (at) for (size_t a=0;a<at->size();a++) cell[a]->itsR=(*at)[a];
    return SolidSetup::Resolve(Lattice_3D(cell, ivec3_t(1,1,1)), kSiBasis, SiGamma(imposed));
}
} // namespace

TEST(SolidSession, ConvergeAtTheReferenceIsTheSinglePointBitForBit)
{
    const auto setup=SiSetup();
    SolidCalculation single(setup, Tight());
    auto S=single.Result(); ASSERT_TRUE(S);

    SolidSession pes(setup, {SCFStage{Tight(), SCFAccelerators::Type::Auto}});
    const PESGeometry ref=pes.Reference();
    ASSERT_TRUE(ref.cell) << "a solid's reference has a cell";
    auto P=pes.ConvergeAt(ref.cell, ref.R);
    ASSERT_TRUE(P) << P.Error().details;
    EXPECT_DOUBLE_EQ(P->energy, S->Energy());
    EXPECT_EQ(P->iterations, S->IterationCount());
    EXPECT_FALSE(P->forces);
}

TEST(SolidSession, ACarriedStepEqualsAFreshRunAtTheNewPositionsInFewerIterations)
{
    const auto setup=SiSetup();
    SolidSession pes(setup, {SCFStage{Tight(), SCFAccelerators::Type::Auto}});
    const PESGeometry ref=pes.Reference();
    std::vector<rvec3_t> R=ref.R;
    auto P0=pes.ConvergeAt(ref.cell, R); ASSERT_TRUE(P0) << P0.Error().details;
    R[1].x+=0.05;   // bohr, Cartesian: one atom off its site
    auto P1=pes.ConvergeAt(ref.cell, R); ASSERT_TRUE(P1) << P1.Error().details;

    SolidCalculation fresh(SiSetup(&R), Tight());   // the same recipe resolved AT the new positions, seeded from scratch
    auto F=fresh.Result(); ASSERT_TRUE(F);
    // Each run's stop at Δρ/N 1e-8 is good to ~1e-8 Ha, so two runs agree to their SUM.  Measured 2026-10-30: carried
    // -7.11461659840 (equal to the same run taken to the floor, 29 iterations), fresh 1.1e-8 above it.
    EXPECT_NEAR(P1->energy, F->Energy(), 2e-8);
    EXPECT_GT(P1->energy, P0->energy) << "off the equilibrium site the energy rises";
    EXPECT_LT(P1->iterations, F->IterationCount()) << "the carried density is the point of the session";
}

TEST(SolidSession, RefusalsAreOutcomesAndAWrongCountThrows)
{
    const auto setup=SiSetup(nullptr, /*imposed*/true);
    SolidSession pes(setup, {SCFStage{Tight(), SCFAccelerators::Type::Auto}});
    const PESGeometry ref=pes.Reference();
    const std::vector<rvec3_t>& R=ref.R;

    auto F=pes.ConvergeAt(ref.cell, R, {.forces=true});
    ASSERT_FALSE(F); EXPECT_EQ(F.Error().why, PESFailure::Why::NoGradient);

    std::vector<rvec3_t> broken=R; broken[1].x+=0.05;
    auto B=pes.ConvergeAt(ref.cell, broken);
    ASSERT_FALSE(B); EXPECT_EQ(B.Error().why, PESFailure::Why::SymmetryBroken) << B.Error().details;

    rmat3d_t strained=*ref.cell; for (size_t i=1;i<=3;i++) for (size_t j=1;j<=3;j++) strained(i,j)*=1.01;   // a 1% isotropic strain
    auto C=pes.ConvergeAt(strained, R);
    ASSERT_FALSE(C); EXPECT_EQ(C.Error().why, PESFailure::Why::CellChange) << C.Error().details;

    EXPECT_THROW((void)pes.ConvergeAt(ref.cell, {R[0]}), std::invalid_argument);
    EXPECT_THROW((void)pes.ConvergeAt(std::nullopt, R), std::invalid_argument) << "a solid's geometry has a cell";

    // A setup over a PRE-BUILT basis cannot follow moved atoms: a programming error, said at once.
    const auto mat=Materials::Get("Si_diamond");
    std::shared_ptr<const BasisSet::Real_BS> mol(BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::SIPP_SR, mat.cell.get()));
    SolidSession fixed(SolidSetup::Resolve(Lattice_3D(*mat.cell, ivec3_t(1,1,1)), mol, SiGamma()), {SCFStage{Tight(), SCFAccelerators::Type::Auto}});
    const PESGeometry fr=fixed.Reference();
    std::vector<rvec3_t> moved=fr.R; moved[1].x+=0.05;
    EXPECT_THROW((void)fixed.ConvergeAt(fr.cell, moved), std::invalid_argument);
}
