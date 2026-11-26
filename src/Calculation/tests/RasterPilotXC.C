// File: Calculation/tests/RasterPilotXC.C  The raster pilot's uniform-XC demand rule (qchem.RasterPilotXC), replayed on the rounds MnO Gamma
// uniform cf 4 eps 1e-6 actually measured (`~/Code/qchem6-runs/anchors/lane-a-pin10/mno_cf4_pilot.log`): finest-level x, L0 x, |dE_xc|.
// L0's x is G_0^2 / p_cap(L0) with p_cap 36 (alpha_max) and a = 10.29 bohr: 60^3 -> 9.32, 72^3 -> 13.42, 80^3 -> 16.58.
#include <gtest/gtest.h>
#include <cmath>
import qchem.RasterPilotXC;

using qchem::UniformXCDemand;

namespace
{
constexpr double kTarget=4e-6;   // MnO's eps 1e-6 x 4 atoms
}

// Round 0 has no history: both arms ask, each on its own axis under the rigorous shape.  The rho arm asks what the old rule asked (11.16).
TEST(UniformXCDemand, AFirstRoundOverRaisesBothArms)
{
    UniformXCDemand u;
    const auto v=u.Round(8.289, 9.32, 2.903e-5, kTarget, /*belowRho*/true);
    EXPECT_TRUE(v.over);
    EXPECT_NEAR(v.needRho, 11.16, 0.01) << "the rho arm keeps the old rule's demand";
    EXPECT_GT(v.needL0, 9.32) << "the L0 arm asks for more than L0 has";
    EXPECT_TRUE(v.refuted.empty());
}

// MnO's logged rounds 1 -> 2: the finest level rose (11.94 -> 15.11), L0 stayed 72^3, and the error re-measured to 4 digits.  The rho arm is
// refuted there, so the L0 arm alone asks -- and asks for more than L0 delivers, which the old rule (rho only) never did.
TEST(UniformXCDemand, AnErrorThatDoesNotFallRefutesTheArmThatMovedAlone)
{
    UniformXCDemand u;
    (void)u.Round(11.94, 13.42, 9.981e-6, kTarget, true);
    const auto v=u.Round(15.11, 13.42, 9.981e-6, kTarget, true);
    EXPECT_EQ(v.refuted, "density raster");
    EXPECT_TRUE(u.RhoRefuted());
    EXPECT_FALSE(v.rhoArm);
    EXPECT_EQ(v.needRho, 0.0);
    EXPECT_GT(v.needL0, 13.42);
    // Refuted stays refuted.
    const auto w=u.Round(18.65, 16.58, 5.385e-6, kTarget, true);
    EXPECT_EQ(w.needRho, 0.0);
    EXPECT_TRUE(w.refuted.empty()) << "refuted once, reported once";
}

// Once L0 has moved, its arm fits the power law on L0's own points: MnO's 9.32 -> 13.42 -> 16.58 read d ~ x^-2.9 on both intervals, so from
// (9.32, 2.903e-5) and (13.42, 9.981e-6) the arm predicts 5.39e-6 at 16.58 (measured 5.385e-6) and asks x ~ 18.4 for 4e-6.  The rigorous
// e^-x shape asks far less (the C-growth the OpenWork row names).
TEST(UniformXCDemand, TheL0ArmFitsItsOwnPowerLaw)
{
    UniformXCDemand u;
    (void)u.Round(8.289, 9.32, 2.903e-5, kTarget, true);
    const auto v=u.Round(11.94, 13.42, 9.981e-6, kTarget, true);
    const double k=std::log(2.903e-5/9.981e-6)/std::log(13.42/9.32);
    EXPECT_NEAR(9.981e-6*std::pow(16.58/13.42, -k), 5.385e-6, 0.1e-6) << "the fit predicts MnO's measured round 3";
    EXPECT_NEAR(v.needL0, 13.42*std::pow(9.981e-6/kTarget, 1.0/k), 1e-9);
    EXPECT_GT(v.needL0, 18.0);
    EXPECT_LT(v.needL0, 19.0);
}

// Above the density's rigorous x only L0's sampling is left (Si Gamma): the rho arm never asks.
TEST(UniformXCDemand, AboveTheRigorousDensityXOnlyL0Asks)
{
    UniformXCDemand u;
    const auto v=u.Round(27.0, 21.0, 5e-5, kTarget, /*belowRho*/false);
    EXPECT_FALSE(v.rhoArm);
    EXPECT_EQ(v.needRho, 0.0);
    EXPECT_GT(v.needL0, 21.0);
    EXPECT_LE(v.needL0, UniformXCDemand::kMaxStep*21.0);
}

// Both arms refuted and still over: nothing the pilot sizes answers it, and the verdict says so (the pilot then stops NOT SETTLED).
TEST(UniformXCDemand, BothArmsRefutedIsUnanswered)
{
    UniformXCDemand u;
    (void)u.Round(10.0, 10.0, 1e-5, kTarget, true);
    (void)u.Round(12.0, 10.0, 1e-5, kTarget, true);   // rho alone moved: refuted
    const auto v=u.Round(12.0, 12.0, 1e-5, kTarget, true);   // L0 alone moved: refuted
    EXPECT_TRUE(u.RhoRefuted());
    EXPECT_TRUE(u.L0Refuted());
    EXPECT_TRUE(v.over);
    EXPECT_TRUE(v.unanswered);
}

// A term that meets asks nothing.
TEST(UniformXCDemand, ATermThatMeetsAsksNothing)
{
    UniformXCDemand u;
    const auto v=u.Round(16.0, 16.0, 1e-6, kTarget, true);
    EXPECT_FALSE(v.over);
    EXPECT_FALSE(v.unanswered);
    EXPECT_EQ(v.needRho, 0.0);
    EXPECT_EQ(v.needL0, 0.0);
}
