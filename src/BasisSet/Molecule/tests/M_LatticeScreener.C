// File: BasisSet/Molecule/tests/M_LatticeScreener.C
//
// THE COLLOCATION TOLERANCE SEAM, AT THE UNIT LEVEL (doc/ScreeningPlan.md).
//
// The rule that sizes every collocation box used to be one experimental bool consulted twice inside the
// innermost geometry code, and it has the worst defect record in the collocation path (ScreeningPlan §2:
// a 4.1 Ha shifted-MP defect that survived in a SECOND unreachable copy, a 2.3 Ha variational collapse
// from FoldScreenMax's signed average, the union-vs-per-component box, the reduced-vs-full tier split).
// Every one of those was found through a full SCF.  Now that the rule is an OBJECT it can be interrogated
// directly, in microseconds, which is the whole practical point of the seam -- so this file pins the
// contract the box walk relies on rather than inferring it from an integrated energy.
//
// ⚠ THE PROPERTY THAT MATTERS MOST here is GeometryOnlyScreener's WEIGHT-INVARIANCE (see
// GeometryOnlyIgnoresTheWeights).  It is what makes ScreeningPlan §5's prize sound -- if the answer does
// not depend on the density then every BoxGeom, chord and bound is iteration-invariant and hoistable into
// the task list.  The seam ships without an IsGeometryOnly() predicate BECAUSE that claim is a testable
// property of the object, not an assertion it has to make about itself.
module;
#include <gtest/gtest.h>
#include <cmath>
module qchem.BasisSet.Molecule.LatticeScreener;   // the unit under test (tests may import internals)
import qchem.Types;
import qchem.Blaze;

using namespace qchem;
using namespace qchem::BasisSet::Molecule;

namespace
{
constexpr double kEps=1e-10;                       // the floor these fixtures screen at
const double kLnEps=-std::log(kEps);

//! A 2x3 weight block with one term per interesting magnitude class.
rmat_t Block(double w00, double w01, double w02, double w10, double w11, double w12)
{
    rmat_t c(2,3,0.0);
    c(0,0)=w00; c(0,1)=w01; c(0,2)=w02;
    c(1,0)=w10; c(1,1)=w11; c(1,2)=w12;
    return c;
}
}

// ---------------------------------------------------------------------------------------------------
// GeometryOnlyScreener -- CP2K's rule.
// ---------------------------------------------------------------------------------------------------

// THE §5 PROPERTY.  Two blocks whose weights differ by orders of magnitude must produce the IDENTICAL
// plan: same tolerance, same mask.  That -- not a predicate the object asserts about itself -- is what
// licenses hoisting the box geometry into the task list, because the density is the only thing that
// changes between SCF iterations.
TEST(M_LatticeScreener, GeometryOnlyIgnoresTheWeights)
{
    const GeometryOnlyScreener s(kEps);
    ScreenPlan a, b;
    s.Screen({0,0,0}, 0.0, Block( 1e-14, 1.0, 3.7,  -2e-9, 1e6, 0.5), a);
    s.Screen({0,0,0}, 0.0, Block(-1e+14, 0.2, 1e-3,  7.0,  1e-7, 42.0), b);
    EXPECT_EQ(a.Tolerance(), b.Tolerance());
    EXPECT_EQ(a.Tolerance(), kEps);
    for (size_t i=0;i<2;i++)
        for (size_t j=0;j<3;j++) EXPECT_EQ(a.Keeps(i,j), b.Keeps(i,j));
}

// "Geometry only" is NOT "no screening": the prefactor test against the flat eps IS the geometry screen
// (CP2K's fixed radius_list).  ScreeningPlan §3's first cut said this screener "never kills", and that
// would have been a silent behaviour change -- every far-offset task would have been walked.
TEST(M_LatticeScreener, GeometryOnlyStillKillsOnThePrefactor)
{
    const GeometryOnlyScreener s(kEps);
    ScreenPlan p;
    s.Screen({3,0,0}, kLnEps*1.01, Block(1,1,1,1,1,1), p);   // whole box sub-eps
    EXPECT_TRUE(p.Dead());
    s.Screen({3,0,0}, kLnEps*0.99, Block(1,1,1,1,1,1), p);   // just inside
    EXPECT_FALSE(p.Dead());
}

// A zero weight is an ABSENT term (dead under the stream fold, or a zero density element), never a
// present one with no weight.  Both rules must read it the same way or the two directions disagree
// about the active set.
TEST(M_LatticeScreener, ZeroWeightIsAbsentUnderBothRules)
{
    const GeometryOnlyScreener g(kEps);
    const DAwareScreener       d(kEps);
    ScreenPlan pg, pd;
    const rmat_t c=Block(0.0, 1.0, 0.0,  0.0, 0.0, 2.0);
    g.Screen({0,0,0}, 0.0, c, pg);
    d.Screen({0,0,0}, 0.0, c, pd);
    for (const ScreenPlan* p : {&pg,&pd})
    {
        EXPECT_FALSE(p->Dead());
        EXPECT_FALSE(p->Keeps(0,0));  EXPECT_TRUE (p->Keeps(0,1));  EXPECT_FALSE(p->Keeps(0,2));
        EXPECT_FALSE(p->Keeps(1,0));  EXPECT_FALSE(p->Keeps(1,1));  EXPECT_TRUE (p->Keeps(1,2));
    }
}

TEST(M_LatticeScreener, AnAllZeroBlockIsDead)
{
    ScreenPlan p;
    GeometryOnlyScreener(kEps).Screen({0,0,0}, 0.0, Block(0,0,0,0,0,0), p);
    EXPECT_TRUE(p.Dead());
    DAwareScreener(kEps).Screen({0,0,0}, 0.0, Block(0,0,0,0,0,0), p);
    EXPECT_TRUE(p.Dead());
}

// ---------------------------------------------------------------------------------------------------
// DAwareScreener -- this tree's default, and a declared CP2K deviation.
// ---------------------------------------------------------------------------------------------------

// eps_ij = eps/|c_ij|, and the UNION is the tightest survivor's -- the largest box, so no term ever rides
// a box smaller than its own rule asked for.  That is what keeps the no-cut discipline under a shared box.
TEST(M_LatticeScreener, DAwareUnionIsTheTightestSurvivor)
{
    const DAwareScreener s(kEps);
    ScreenPlan p;
    s.Screen({0,0,0}, 0.0, Block(1e-4, 1e-2, 0.0,  0.0, 1e-6, 0.0), p);
    // |c|=1e-2 is the largest weight present -> eps/1e-2 is the SMALLEST tolerance -> the biggest box.
    EXPECT_DOUBLE_EQ(p.Tolerance(), kEps/1e-2);
    EXPECT_TRUE(p.Keeps(0,0)); EXPECT_TRUE(p.Keeps(0,1)); EXPECT_TRUE(p.Keeps(1,1));
}

// CLAMPED AT THE FLOOR.  A |c|>1 must never widen a box BEYOND the geometry screen: the (shell pair,
// offset) task list is enumerated once at the floor, so a tolerance below it would want boxes -- and
// offsets -- that were never built.  The walks assert Floor()>=kDensityEps for the same reason.
TEST(M_LatticeScreener, DAwareNeverAnswersBelowItsFloor)
{
    const DAwareScreener s(kEps);
    ScreenPlan p;
    s.Screen({0,0,0}, 0.0, Block(1e6, 3.0, 0.0, 0.0, 0.0, 0.0), p);
    EXPECT_DOUBLE_EQ(p.Tolerance(), kEps);
    EXPECT_EQ(s.Floor(), kEps);
    EXPECT_EQ(GeometryOnlyScreener(kEps).Floor(), kEps);
}

// ⛔ THE SHIFTED-MP DEFECT'S UNIT GATE.  The rule reads a MAGNITUDE: c and -c are the same term as far as
// the box is concerned.  The 4.1 Ha defect was a screen that read a signed real part and discarded every
// odd-offset term at quarter-integer k; a sign-sensitive tolerance is that bug's shape.
TEST(M_LatticeScreener, DAwareToleranceIsSignBlind)
{
    const DAwareScreener s(kEps);
    ScreenPlan pos, neg;
    s.Screen({1,0,0}, 0.0, Block( 3e-3, -7e-5, 0.0, 0.0, 0.0, 0.0), pos);
    s.Screen({1,0,0}, 0.0, Block(-3e-3,  7e-5, 0.0, 0.0, 0.0, 0.0), neg);
    EXPECT_EQ(pos.Tolerance(), neg.Tolerance());
    EXPECT_EQ(pos.Keeps(0,0), neg.Keeps(0,0));
    EXPECT_EQ(pos.Keeps(0,1), neg.Keeps(0,1));
}

// A small-|c| term is dropped whole while its large-|c| sibling survives on the same offset -- the
// per-term kill that the union box does NOT paper over.
TEST(M_LatticeScreener, DAwareKillsPerTermNotPerBlock)
{
    const DAwareScreener s(kEps);
    ScreenPlan p;
    // pf sits BETWEEN the two terms' -ln(eps/|c|) thresholds.  Note which way the rule runs: a SMALL |c|
    // buys a LOOSER tolerance (eps/|c| is a bigger number) and hence a smaller box, so it is the negligible
    // term that dies first.  |c|=1e-3 -> eps_ij=1e-7, threshold 16.1; |c|=1e-16 -> eps_ij=1e+6, threshold
    // -13.8, i.e. dead at any pf at all.
    const double pf=10.0;
    s.Screen({2,1,0}, pf, Block(1e-3, 1e-16, 0.0, 0.0, 0.0, 0.0), p);
    EXPECT_FALSE(p.Dead());
    EXPECT_TRUE (p.Keeps(0,0));
    EXPECT_FALSE(p.Keeps(0,1));
    EXPECT_DOUBLE_EQ(p.Tolerance(), kEps/1e-3);
}

// THE STATIC-SWEEP IDENTITY.  IntegratePotential with no screenD presents a UNIT weight, and the two
// rules must then agree exactly -- that is what makes the local-PP and V_loc-long sweeps bit-identical
// across the seam whichever screener a run chose.
TEST(M_LatticeScreener, UnitWeightsAgreeAcrossBothRules)
{
    ScreenPlan g, d;
    const rmat_t ones=Block(1,1,1,1,1,1);
    GeometryOnlyScreener(kEps).Screen({0,1,2}, 5.0, ones, g);
    DAwareScreener      (kEps).Screen({0,1,2}, 5.0, ones, d);
    EXPECT_EQ(g.Tolerance(), d.Tolerance());
    EXPECT_EQ(g.Tolerance(), kEps);
    for (size_t i=0;i<2;i++)
        for (size_t j=0;j<3;j++) EXPECT_EQ(g.Keeps(i,j), d.Keeps(i,j));
}

// The plan is REUSED across a shell pair's offsets, so a stale mask bit from the previous offset would be
// a silent extra term in the box.  Reset must clear it -- and must survive a shape change.
TEST(M_LatticeScreener, PlanReuseDoesNotLeakTheLastOffsetsMask)
{
    const DAwareScreener s(kEps);
    ScreenPlan p;
    s.Screen({0,0,0}, 0.0, Block(1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2), p);
    EXPECT_TRUE(p.Keeps(1,2));
    s.Screen({1,0,0}, 0.0, Block(1e-2, 0.0, 0.0, 0.0, 0.0, 0.0), p);   // same shape, one live term
    EXPECT_TRUE (p.Keeps(0,0));
    EXPECT_FALSE(p.Keeps(1,2));
    rmat_t wide(3,4,0.0); wide(2,3)=1.0;                                // a different shape entirely
    s.Screen({1,0,0}, 0.0, wide, p);
    EXPECT_TRUE (p.Keeps(2,3));
    EXPECT_FALSE(p.Keeps(0,0));
}

// ⛔ THE DIAGONAL-SEED PIN (user, 2026-09-04: "the Screener system needs to properly handle diagonal seed
// densities").  The SAD seed's density matrix is DIAGONAL, so a gather that forwarded |D_ij| straight
// through would hand the screener a block that is ZERO everywhere off the diagonal.  A 0 here means
// STRUCTURALLY ABSENT, so every off-diagonal term would be dropped and every off-diagonal h_ij would come
// out 0 -- and h is DIAGONALIZED to make the next density, so the truncation is SELF-FULFILLING: a pair
// with no density can never acquire any.  (h_ij is not weighted by D_ij; only the ENERGY Tr(D h) is blind
// to these terms, and the energy is not the only consumer.)
//
// ⇒ THE CALLER'S CONTRACT is to send the UNIT weight for "no information", never 0.  This test pins both
// halves: what the bug shape would have done, and what the contract does instead.
TEST(M_LatticeScreener, DiagonalSeedKeepsItsOffDiagonals)
{
    const DAwareScreener d(kEps);
    const GeometryOnlyScreener g(kEps);
    ScreenPlan p;

    // (1) THE BUG SHAPE: |D| forwarded raw from a diagonal seed -- off-diagonals arrive as 0.
    d.Screen({0,0,0}, 0.0, Block(0.5, 0.0, 0.0,  0.0, 0.5, 0.0), p);
    EXPECT_TRUE (p.Keeps(0,0)) << "the diagonal survives either way";
    EXPECT_FALSE(p.Keeps(0,1)) << "a raw 0 IS read as structurally absent -- this is why the caller must "
                                  "not send one for a merely-zero density";

    // (2) THE CONTRACT: "no information" travels as the unit weight, so the term is KEPT at the floor.
    const rmat_t seed=Block(0.5, 1.0, 1.0,  1.0, 0.5, 1.0);
    const LatticeScreener* rules[]={&d,&g};
    for (const LatticeScreener* s : rules)
    {
        s->Screen({0,0,0}, 0.0, seed, p);
        EXPECT_FALSE(p.Dead());
        for (size_t a=0;a<2;a++)
            for (size_t b=0;b<3;b++)
                EXPECT_TRUE(p.Keeps(a,b)) << "every present term must ride the box: a=" << a << " b=" << b;
        // The union is the floor: the unit weights ask for eps, and the 0.5 diagonal asks for the LOOSER
        // eps/0.5, so the tightest present tolerance -- the one that sizes the shared box -- is eps itself.
        EXPECT_DOUBLE_EQ(p.Tolerance(), kEps);
    }
}

// The other half of the same contract: on NONZERO weights a geometry-only rule answers IDENTICALLY whether
// the caller sends real magnitudes or unit weights.  That equivalence is what a future cross-k gather memo
// would rest on -- withholding the density screen could then be bit-identical rather than a widening --
// and it holds ONLY away from zero, which is exactly the gap doc/OpenWork.md's k-scaling entry records.
TEST(M_LatticeScreener, GeometryOnlyIsBlindToNonzeroWeights)
{
    const GeometryOnlyScreener g(kEps);
    ScreenPlan withMagnitudes, withUnits;
    g.Screen({1,0,0}, 3.0, Block(0.5, 2e-4, 7.0,  1e-9, 0.5, 3.3), withMagnitudes);
    g.Screen({1,0,0}, 3.0, Block(1.0, 1.0,  1.0,  1.0,  1.0, 1.0), withUnits);
    EXPECT_EQ(withMagnitudes.Tolerance(), withUnits.Tolerance());
    for (size_t a=0;a<2;a++)
        for (size_t b=0;b<3;b++) EXPECT_EQ(withMagnitudes.Keeps(a,b), withUnits.Keeps(a,b));
}
