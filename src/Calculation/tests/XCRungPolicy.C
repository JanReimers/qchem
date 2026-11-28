// File: src/Calculation/tests/XCRungPolicy.C  The uniform XC feed's rung route: the structural rule (prod m - prod x < 1), Auto's pick on the
// measured ladders, a stated route winning, and the no-rung case (qchem.XCRungPolicy; user ruling 2026-11-26, "Auto, with B as the fallback").
#include "gtest/gtest.h"
import qchem.XCRungPolicy;
import qchem.Types;

using namespace qchem;
using F=BasisSet::Gaussian::XCRungFeed;
using policy::XCRungRoute;

namespace
{
BasisSet::Gaussian::RungShape Shape(int N0, int Nown)
{
    BasisSet::Gaussian::RungShape s;
    s.N0=ivec3_t(N0,N0,N0); s.Nown=ivec3_t(Nown,Nown,Nown);
    const int m=(Nown+N0-1)/N0;
    s.m=ivec3_t(m,m,m);
    return s;
}
} //anon

// The ladders the cost model was measured on (XCRungPolicy.C's calibration table), and Auto's pick on each.
TEST(XCRungPolicy, AutoOnTheMeasuredLadders)
{
    EXPECT_EQ(policy::DecideXCRung(XCRungRoute::Auto, Shape(90,81)).feed, F::RungMultiple) << "MnO 20/8.3: the rung coarser than L0 (m = 1)";
    EXPECT_EQ(policy::DecideXCRung(XCRungRoute::Auto, Shape(72,50)).feed, F::RungMultiple) << "NaF x32/8: the rung coarser than L0";
    EXPECT_EQ(policy::DecideXCRung(XCRungRoute::Auto, Shape(25,36)).feed, F::PairsOnL0)    << "NaF default: x 1.44 per axis, m 2 would be 2.7x the points";
    EXPECT_EQ(policy::DecideXCRung(XCRungRoute::Auto, Shape(15,24)).feed, F::PairsOnL0)    << "Si Gamma default: x 1.6";
    EXPECT_EQ(policy::DecideXCRung(XCRungRoute::Auto, Shape(24,48)).feed, F::RungMultiple) << "an exact multiple: C costs nothing extra";
    EXPECT_EQ(policy::DecideXCRung(XCRungRoute::Auto, Shape(36,36)).feed, F::RungMultiple) << "equal rasters: C is the plain shared collocation";
}

// The per-axis threshold x = 7^(1/3) = 1.913 on an isotropic raster, from both sides, and the formula's two extras.
TEST(XCRungPolicy, ThresholdAndExtras)
{
    EXPECT_EQ(policy::DecideXCRung(XCRungRoute::Auto, Shape(100,192)).feed, F::RungMultiple);   // 1.92^3 = 7.08: 8 - 7.08 < 1
    EXPECT_EQ(policy::DecideXCRung(XCRungRoute::Auto, Shape(100,190)).feed, F::PairsOnL0);      // 1.90^3 = 6.86: 8 - 6.86 > 1
    EXPECT_EQ(policy::DecideXCRung(XCRungRoute::Auto, Shape(100,297)).feed, F::RungMultiple);   // k = 3: 2.97^3 = 26.20: 27 - 26.20 < 1
    EXPECT_EQ(policy::DecideXCRung(XCRungRoute::Auto, Shape(100,295)).feed, F::PairsOnL0);      // 2.95^3 = 25.67: 27 - 25.67 > 1
    const policy::XCRungCost c=policy::RungCost(Shape(25,36));
    const double x=(36.0/25)*(36.0/25)*(36.0/25);
    EXPECT_NEAR(c.x, x, 1e-12);
    EXPECT_NEAR(c.relB, 1.0/x, 1e-12);
    EXPECT_NEAR(c.relC, (8.0-x)/x, 1e-12);
    for (int N0 : {15,24,25,36,54,72,90})   // the rule IS the comparison of the two extras, whatever R is
        for (int Nr : {20,24,36,48,50,64,81,100,125})
        {
            const policy::XCRungCost k=policy::RungCost(Shape(N0,Nr));
            EXPECT_EQ(policy::XCRungCheaperRoute(k)==F::RungMultiple, k.relC<k.relB) << N0 << " / " << Nr;
        }
}

// A stated route wins; one predicted at more than twice the other's extra runs with a warning naming Auto's pick.
TEST(XCRungPolicy, StatedWins)
{
    const auto b=policy::DecideXCRung(XCRungRoute::PairsOnL0, Shape(90,81));
    EXPECT_EQ(b.feed, F::PairsOnL0);
    EXPECT_TRUE(b.warn) << "B at x 0.73 costs 1.37 R against C's 0.37 R";
    EXPECT_EQ(b.autoPick, F::RungMultiple);
    const auto c=policy::DecideXCRung(XCRungRoute::RungMultiple, Shape(25,36));
    EXPECT_EQ(c.feed, F::RungMultiple);
    EXPECT_TRUE(c.warn) << "C at x 1.44 costs 1.68 R against B's 0.33 R";
    const auto quiet=policy::DecideXCRung(XCRungRoute::RungMultiple, Shape(24,48));
    EXPECT_FALSE(quiet.warn);
}

// No rung: nothing to choose, and the stated name is honoured (the feed is L0 plus the coarse levels either way).
TEST(XCRungPolicy, NoRung)
{
    BasisSet::Gaussian::RungShape s; s.N0=ivec3_t(36,36,36);
    const auto d=policy::DecideXCRung(XCRungRoute::Auto, s);
    EXPECT_FALSE(d.hasRung);
    EXPECT_EQ(d.feed, F::PairsOnL0);
    EXPECT_EQ(policy::DecideXCRung(XCRungRoute::RungMultiple, s).feed, F::RungMultiple);
}
