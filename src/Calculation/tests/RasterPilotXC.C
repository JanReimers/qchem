// File: Calculation/tests/RasterPilotXC.C  The raster pilot's uniform-XC demand rule (qchem.RasterPilotXC), replayed on the rounds MnO Gamma
// uniform cf 4 eps 1e-6 actually measured (`~/Code/qchem6-runs/anchors/lane-a-pin10/mno_cf4_pilot.log`): L0 x, |dE_xc|.
// L0's x is G_0^2 / p_cap(L0) with p_cap 36 (alpha_max) and a = 10.29 bohr: 60^3 -> 9.32, 72^3 -> 13.42, 80^3 -> 16.58.
#include <gtest/gtest.h>
#include <cmath>
import qchem.RasterPilotXC;

using qchem::UniformXCDemand;

namespace
{
constexpr double kTarget=4e-6;   // MnO's eps 1e-6 x 4 atoms
}

// Round 0 has no history: the arm asks under the rigorous shape, for more than L0 has.
TEST(UniformXCDemand, AFirstRoundOverAsksForMoreL0)
{
    UniformXCDemand u;
    const auto v=u.Round(9.32, 2.903e-5, kTarget);
    EXPECT_TRUE(v.over);
    EXPECT_GT(v.needL0, 9.32);
    EXPECT_LE(v.needL0, UniformXCDemand::kMaxStep*9.32);
    EXPECT_FALSE(v.refuted);
}

// MnO's logged rounds 1 -> 2: L0 stayed 72^3 (FFT rounding) and the error re-measured to 4 digits.  An L0 that did not move refutes
// nothing, so the arm still asks and the term is not stalled.
TEST(UniformXCDemand, AnUnmovedL0RefutesNothing)
{
    UniformXCDemand u;
    (void)u.Round(13.42, 9.981e-6, kTarget);
    const auto v=u.Round(13.42, 9.981e-6, kTarget);
    EXPECT_FALSE(v.refuted);
    EXPECT_FALSE(u.L0Refuted());
    EXPECT_GT(v.needL0, 13.42);
    EXPECT_FALSE(u.Stalled()) << "a round whose L0 did not move is not a stall";
}

// Once L0 has moved, the arm fits the power law on L0's own points: MnO's 9.32 -> 13.42 -> 16.58 read d ~ x^-2.9 on both intervals, so from
// (9.32, 2.903e-5) and (13.42, 9.981e-6) the arm predicts 5.39e-6 at 16.58 (measured 5.385e-6) and asks x ~ 18.4 for 4e-6.  The rigorous
// e^-x shape asks far less (the C-growth the OpenWork row names).
TEST(UniformXCDemand, TheL0ArmFitsItsOwnPowerLaw)
{
    UniformXCDemand u;
    (void)u.Round(9.32, 2.903e-5, kTarget);
    const auto v=u.Round(13.42, 9.981e-6, kTarget);
    const double k=std::log(2.903e-5/9.981e-6)/std::log(13.42/9.32);
    EXPECT_NEAR(9.981e-6*std::pow(16.58/13.42, -k), 5.385e-6, 0.1e-6) << "the fit predicts MnO's measured round 3";
    EXPECT_NEAR(v.needL0, 13.42*std::pow(9.981e-6/kTarget, 1.0/k), 1e-9);
    EXPECT_GT(v.needL0, 18.0);
    EXPECT_LT(v.needL0, 19.0);
}

// L0 moved and the error did not fall, still over: raising L0 does not answer this error.  Refuted once, reported once, and STALLED.
TEST(UniformXCDemand, AMovedL0WhoseErrorDidNotFallIsStalled)
{
    UniformXCDemand u;
    (void)u.Round(10.0, 1e-5, kTarget);
    const auto v=u.Round(12.0, 1e-5, kTarget);
    EXPECT_TRUE(v.refuted);
    EXPECT_TRUE(u.L0Refuted());
    EXPECT_EQ(v.needL0, 0.0);
    EXPECT_TRUE(u.Stalled());
    const auto w=u.Round(14.0, 1e-5, kTarget);
    EXPECT_FALSE(w.refuted) << "refuted once, reported once";
    EXPECT_TRUE(u.Stalled());
}

// A term that meets asks nothing and is not stalled.
TEST(UniformXCDemand, ATermThatMeetsAsksNothing)
{
    UniformXCDemand u;
    const auto v=u.Round(16.0, 1e-6, kTarget);
    EXPECT_FALSE(v.over);
    EXPECT_EQ(v.needL0, 0.0);
    EXPECT_FALSE(u.Stalled());
}
