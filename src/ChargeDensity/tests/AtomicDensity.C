// File: UnitTests/AtomicDensityUT.C  Tests for the SAD atomic-density database reader + interpolation.
#include "gtest/gtest.h"
#include <memory>

import qchem.ChargeDensity.AtomicDensity;   // GetAtomicDensity, RadialDensity, RecentredAtomicDensity
import qchem.Types;                          // rvec3_t
using namespace qchem;

using namespace qchem::ChargeDensity;

// The stored grid integrates 4*pi*int r^2 rho dr to the element's electron count (the generator's charge
// check), and the reader's own trapezoid reproduces it.
TEST(AtomicDensity, ChargeIntegratesToNelec)
{
    for (int Z : {1,2,6,8,10})
    {
        RadialDensity rad = GetAtomicDensity(Z, "LDA");
        // Relative bound: the log-radial quadrature has a small (~0.04%) error from the sharp core, which
        // grows with Z in absolute terms -- a Z-scaled tolerance is the honest "integrates to ~Nelec" check.
        EXPECT_NEAR(rad.Charge(), double(Z), Z*2e-3) << "Z=" << Z;
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

// ---------------- spin-resolved pairs (doc/SCFSeedingPlan.md sec 10, the spin-SAD tables) ----------------

static constexpr const char* vdb = "atomic_valence_densities.json";

// The Hund entries carry an up-majority pair whose sum is the stored rho and whose difference integrates
// to the Hund moment 2S: Mn q7 neutral (4s^2 3d^5, 2S=5), Mn2+ (3d^5, 2S=5 -- the physical TM-cation
// removal order, 4s first), O neutral (2p^4 triplet, 2S=2).
TEST(AtomicDensity, SpinPairChargesAndMoments)
{
    struct Case { int Z, Nelec; double moment; };
    for (Case c : {Case{25,7,5.0}, Case{25,5,5.0}, Case{8,6,2.0}})
    {
        ASSERT_TRUE(HasAtomicSpinPair(c.Z, "LDA", vdb, c.Nelec)) << "Z="<<c.Z<<" Nelec="<<c.Nelec;
        auto [up,dn] = GetAtomicSpinPair(c.Z, "LDA", vdb, c.Nelec);
        // Nelec-scaled bound, same reason as ChargeIntegratesToNelec: the reader's log-mesh quadrature
        // carries a small sharp-core error on top of the generator's own ~1e-3.
        EXPECT_NEAR(up.Charge()+dn.Charge(), double(c.Nelec), c.Nelec*1e-3) << "Z="<<c.Z; // sum = charge state
        EXPECT_NEAR(up.Charge()-dn.Charge(), c.moment,        c.Nelec*1e-3) << "Z="<<c.Z; // difference = 2S
        EXPECT_GE(up.Charge(), dn.Charge());                                        // up-majority convention

        // The pair sums to the spin-summed rho the unpolarized readers see (same entry, same grid).  The
        // three arrays are emitted independently at 12 significant digits, so consistency is RELATIVE.
        RadialDensity tot = GetAtomicDensity(c.Z, "LDA", vdb, c.Nelec);
        for (double r : {0.05, 0.3, 1.0, 3.0})
            EXPECT_NEAR(up(r)+dn(r), tot(r), 1e-9*tot(r)+1e-15) << "Z="<<c.Z<<" r="<<r;
    }
}

// Closed-shell (or spin-agnostic legacy) entries carry no pair: Has is false and Get throws.
TEST(AtomicDensity, SpinPairAbsentForClosedShells)
{
    EXPECT_FALSE(HasAtomicSpinPair( 8, "LDA", vdb, 8));   // O2- (closed p^6, generated unpolarized)
    EXPECT_FALSE(HasAtomicSpinPair(14, "LDA", vdb, 4));   // Si (legacy spin-agnostic entry)
    EXPECT_THROW(GetAtomicSpinPair(14, "LDA", vdb, 4), std::runtime_error);
}
