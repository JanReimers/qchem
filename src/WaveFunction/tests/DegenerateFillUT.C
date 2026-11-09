// File: WaveFunction/tests/DegenerateFillUT.C  The kT = 0 fill kernel (D-DEGEN-FILL, user ruling 2026-11-08):
// qchem::Orbitals::ShareDegenerateLevels, which both kT = 0 paths (the per-block fill and the ranked reservoir) call.
//
// What these pin: a partly filled DEGENERATE level is shared equally (so the occupations, hence D, do not depend on which
// member the eigensolver listed first); whole levels below it fill whole; the level width is the tolerance, measured from
// the level's LOWEST energy; and a STATED order (a MOM rank, or levelTol <= 0 for a held GDM trial) is the old orbital-by-orbital
// fill, never shared.
#include "gtest/gtest.h"
#include <stdexcept>
import qchem.Orbitals;
import qchem.Types;

using qchem::Orbitals::ShareDegenerateLevels;

namespace
{
rvec_t V(std::initializer_list<double> x) {rvec_t v(x.size()); size_t i=0; for (double d : x) v[i++]=d; return v;}
double Sum(const rvec_t& v) {double s=0.0; for (size_t i=0;i<v.size();++i) s+=v[i]; return s;}
}

// PW_NaF Γ's case: 4 electrons into a 6-fold level (three spatial orbitals, g = 2) -- shared 4/3 each, never 2+2+0.
TEST(DegenerateFill, APartlyFilledLevelIsSharedEqually)
{
    const rvec_t e=V({-0.5, 0.1, 0.1, 0.1, 0.4}), g=V({2,2,2,2,2});
    const rvec_t occ=ShareDegenerateLevels(e, g, rvec_t(), 6.0, 1e-4);
    EXPECT_DOUBLE_EQ(occ[0], 2.0) << "the level below fills whole";
    for (size_t i : {1,2,3}) EXPECT_NEAR(occ[i], 4.0/3.0, 1e-15) << i;
    EXPECT_EQ(occ[4], 0.0);
    EXPECT_NEAR(Sum(occ), 6.0, 1e-14);
}

// The order the eigensolver lists a level's members in must not matter: permute them, get the same occupations back.
TEST(DegenerateFill, TheFillIsInvariantToTheOrderInsideALevel)
{
    const rvec_t e1=V({0.1, 0.1+2e-5, 0.1-3e-5}), e2=V({0.1-3e-5, 0.1, 0.1+2e-5}), g=V({2,2,2});
    const rvec_t a=ShareDegenerateLevels(e1, g, rvec_t(), 1.0, 1e-4), b=ShareDegenerateLevels(e2, g, rvec_t(), 1.0, 1e-4);
    for (size_t i=0;i<3;++i) { EXPECT_NEAR(a[i], 1.0/3.0, 1e-15); EXPECT_NEAR(b[i], 1.0/3.0, 1e-15); }
}

// GPW_Al Γ's case, unequal capacities: 1 electron over a level of g = {2, 4} gets 1/3 and 2/3 (g_i x left / G).
TEST(DegenerateFill, TheShareIsPerUnitCapacity)
{
    const rvec_t occ=ShareDegenerateLevels(V({0.0, 0.0}), V({2, 4}), rvec_t(), 1.0, 1e-4);
    EXPECT_NEAR(occ[0], 1.0/3.0, 1e-15);
    EXPECT_NEAR(occ[1], 2.0/3.0, 1e-15);
}

// The level width is measured from the level's LOWEST member: 0, 0.6e-4, 1.2e-4 at tol 1e-4 is two levels, not one chain.
TEST(DegenerateFill, TheWidthIsMeasuredFromTheLevelsLowestEnergy)
{
    const rvec_t occ=ShareDegenerateLevels(V({0.0, 0.6e-4, 1.2e-4}), V({2,2,2}), rvec_t(), 3.0, 1e-4);
    EXPECT_NEAR(occ[0], 1.5, 1e-15);
    EXPECT_NEAR(occ[1], 1.5, 1e-15);
    EXPECT_EQ(occ[2], 0.0);
}

// A level the budget fills EXACTLY stays integer: an insulator's closed shells are untouched by the sharing.
TEST(DegenerateFill, AWholeLevelStaysInteger)
{
    const rvec_t occ=ShareDegenerateLevels(V({-1.0, 0.2, 0.2, 0.2, 0.9}), V({2,2,2,2,2}), rvec_t(), 8.0, 1e-4);
    EXPECT_EQ(occ[0], 2.0);
    for (size_t i : {1,2,3}) EXPECT_EQ(occ[i], 2.0);
    EXPECT_EQ(occ[4], 0.0);
}

// A STATED order wins exactly: a MOM rank picks orbitals one by one, even inside a degenerate level (no sharing).
TEST(DegenerateFill, AStatedRankIsFilledExactlyNeverShared)
{
    const rvec_t e=V({-0.2, 0.3, 0.3}), g=V({2,2,2});
    const rvec_t occ=ShareDegenerateLevels(e, g, V({0.1, 0.9, 0.5}), 3.0, 1e-4);
    EXPECT_EQ(occ[1], 2.0) << "the highest overlap first";
    EXPECT_EQ(occ[2], 1.0);
    EXPECT_EQ(occ[0], 0.0);
}

// levelTol <= 0: orbital by orbital in stored order (no rank) or rank order -- the held fill, which must reproduce the block
// a geodesic rotated, degenerate or not.
TEST(DegenerateFill, ZeroWidthIsTheOrbitalByOrbitalFill)
{
    const rvec_t e=V({0.1, 0.1, 0.1}), g=V({2,2,2});
    const rvec_t a=ShareDegenerateLevels(e, g, rvec_t(), 4.0, 0.0);
    EXPECT_EQ(a[0], 2.0); EXPECT_EQ(a[1], 2.0); EXPECT_EQ(a[2], 0.0);
    const rvec_t b=ShareDegenerateLevels(e, g, V({0.2, 0.1, 0.9}), 3.0, 0.0);
    EXPECT_EQ(b[2], 2.0); EXPECT_EQ(b[0], 1.0); EXPECT_EQ(b[1], 0.0);
}

//! A rankless orbital-by-orbital fill keeps the STORED order: GDM stores its rotated orbitals with the occupied set leading, and
//! their diag(F) need not ascend (2026-11-08 panel: an energy re-sort moved the MnO dimer's fill off GDM's manifold).  A caller
//! whose orbitals arrive in another order states it as a rank: the molecular reservoir's aufbau ranks by -energy across irreps
//! (2026-11-08 sweep: water C2v HF filled 5A1 over 1B1 when it did not).
TEST(DegenerateFill, ZeroWidthKeepsTheStoredOrderAndARankStatesAnother)
{
    const rvec_t e=V({-0.3, 0.1, -0.7, -0.5}), g=V({2,2,2,2});
    const rvec_t stored=ShareDegenerateLevels(e, g, rvec_t(), 4.0, 0.0);
    EXPECT_EQ(stored[0], 2.0); EXPECT_EQ(stored[1], 2.0); EXPECT_EQ(stored[2], 0.0); EXPECT_EQ(stored[3], 0.0);
    const rvec_t aufbau=ShareDegenerateLevels(e, g, V({0.3, -0.1, 0.7, 0.5}), 4.0, 0.0);   // rank = -e
    EXPECT_EQ(aufbau[2], 2.0); EXPECT_EQ(aufbau[3], 2.0); EXPECT_EQ(aufbau[0], 0.0); EXPECT_EQ(aufbau[1], 0.0);
}

TEST(DegenerateFill, MoreElectronsThanCapacityLeaveTheRestToTheCaller)
{
    const rvec_t occ=ShareDegenerateLevels(V({0.0, 0.5}), V({2,2}), rvec_t(), 5.0, 1e-4);
    EXPECT_EQ(Sum(occ), 4.0);
    EXPECT_THROW((void)ShareDegenerateLevels(V({0.0}), V({2,2}), rvec_t(), 1.0, 1e-4), std::logic_error);
    EXPECT_THROW((void)ShareDegenerateLevels(V({0.0}), V({2}), rvec_t(), -1.0, 1e-4), std::logic_error);
}
