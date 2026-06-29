// File: UnitTests/AtomicDensityUT.C  Tests for the SAD atomic-density database reader + interpolation.
#include "gtest/gtest.h"
#include <memory>

import qchem.ChargeDensity.AtomicDensity;   // GetAtomicDensity, RadialDensity, RecentredAtomicDensity
import qchem.Types;                          // rvec3_t

using namespace qchem::ChargeDensity;

// The stored grid integrates 4*pi*int r^2 rho dr to the element's electron count (the generator's charge
// check), and the reader's own trapezoid reproduces it.
TEST(AtomicDensity, ChargeIntegratesToNelec)
{
    for (int Z : {1,2,6,8,10})
    {
        RadialDensity rad = GetAtomicDensity(Z, "LDA");
        EXPECT_NEAR(rad.Charge(), double(Z), 3e-3) << "Z=" << Z;
    }
}

// rho(r) is positive, peaks at the core, and decays outward; node interpolation stays in range.
TEST(AtomicDensity, RadialShape)
{
    RadialDensity rad = GetAtomicDensity(8, "LDA");   // oxygen
    EXPECT_GT(rad(0.0), 0.0);
    EXPECT_GT(rad(0.1), rad(1.0));
    EXPECT_GT(rad(1.0), rad(5.0));
    EXPECT_NEAR(rad(1e6), 0.0, 1e-30);                // clamps to 0 far outside the grid
}

// Recentring: rho(|r-R|) -- value at the nucleus equals rad(0); an offset of d equals rad(d).
TEST(AtomicDensity, RecentredMatchesRadial)
{
    auto rad = std::make_shared<const RadialDensity>(GetAtomicDensity(8, "LDA"));
    rvec3_t R(0.3, -1.2, 2.0);
    RecentredAtomicDensity cd(rad, R);
    EXPECT_DOUBLE_EQ(cd(R), (*rad)(0.0));
    const double d = 0.75;
    EXPECT_DOUBLE_EQ(cd(R + rvec3_t(0,0,d)), (*rad)(d));
    EXPECT_DOUBLE_EQ(cd(R + rvec3_t(d,0,0)), (*rad)(d));
}

TEST(AtomicDensity, MissingElementThrows)
{
    EXPECT_THROW(GetAtomicDensity(123, "LDA"), std::runtime_error);
    EXPECT_THROW(GetAtomicDensity(8, "NoSuchFunctional"), std::runtime_error);
}
