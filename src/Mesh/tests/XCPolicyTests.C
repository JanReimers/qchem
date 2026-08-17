// File: src/Mesh/tests/XCPolicyTests.C  The XC-quadrature policy: Nyquist arithmetic, cost model, selector.
//
// V1.26.  Two things are pinned here that are easy to get wrong silently:
//   (1) UniformDivisions/UniformCutoff are INVERSES.  They are the bridge between MeshParams' two spellings
//       of one quantity (eCut and nUniform), so every resolution question in the project routes through
//       them; if they drifted apart, "is this mesh fine enough?" would get a different answer from "how fine
//       is this mesh?".
//   (2) The selector's DIRECTION.  The whole point is that the crossover moves with alpha_max, so a soft
//       system must favour uniform and a sharp one must favour Becke -- and Becke's cost must NOT move with
//       alpha_max at all, which is what makes the crossover exist.
#include "gtest/gtest.h"
#include <cmath>
import qchem.Mesh.XCPolicy;
import qchem.Math;                // Pi

using namespace qchem;
using namespace qchem::qcMesh;

//================================================================================================
//  The Nyquist bridge
//================================================================================================

// The formula this replaced was a literal inside UnitCell::CreateIntegrationMesh; pin it so the extraction
// is provably behaviour-preserving rather than merely plausible.
TEST(XCPolicy, UniformDivisionsMatchesTheNyquistFormula)
{
    for (double a : {6.0, 8.7, 10.26})
        for (double E : {4.0, 16.0, 80.0})
            EXPECT_EQ(UniformDivisions(a,E), int(std::ceil(2.0*a*std::sqrt(2.0*E)/Pi)))
                << "a=" << a << " E=" << E;
}

// The inverse must never UNDERSTATE what a mesh resolves -- an optimistic UniformCutoff would silently
// suppress the under-resolution warning it exists to raise.
TEST(XCPolicy, UniformCutoffInvertsUniformDivisions)
{
    for (double a : {6.0, 8.7, 10.26})
        for (double E : {4.0, 16.0, 80.0})
        {
            const int n=UniformDivisions(a,E);
            EXPECT_GE(UniformCutoff(a,n), E) << "a=" << a << " E=" << E << " n=" << n
                << ": the mesh sized for E must report resolving at least E";
            // and it must be TIGHT -- one division coarser must not still suffice (up to the ceiling).
            if (n>1) EXPECT_LT(UniformCutoff(a,n-1), UniformCutoff(a,n));
        }
}

// The measured hazard behind V1.26: nUniform's default of 20 is basis-blind, and at typical cell edges it
// resolves only alpha_max ~ 2-3.  Si/sipp (alpha_max=2) just squeaks under it, which is exactly why the
// uniform route looks fine on Si and on nothing else.
TEST(XCPolicy, DefaultNUniformResolvesOnlyVerySoftBases)
{
    const MeshParams def;                                    // whatever the shipped default is
    for (double a : {8.7, 10.26})
    {
        const double alphaOk=UniformCutoff(a,def.nUniform)/2.0;   // E = 2*alpha_max
        EXPECT_LT(alphaOk, 4.0) << "a=" << a << ": nUniform=" << def.nUniform << " resolves alpha_max="
                                << alphaOk << " -- if this ever exceeds ~4 the warning's premise changed";
        EXPECT_GT(alphaOk, 1.5);
    }
}

//================================================================================================
//  The cost model
//================================================================================================

// The two scalings are the whole argument for a crossover existing: uniform tracks alpha_max^(3/2), Becke
// does not track it at all.  Pin BOTH halves -- a Becke cost that quietly acquired an alpha_max dependence
// would destroy the selector without failing any other test.
TEST(XCPolicy, UniformCostGrowsWithSharpnessAndBeckeDoesNot)
{
    XCMeshSharpness soft {.cellEdge=10.26, .nAtoms=2, .alphaMax=2.0};
    XCMeshSharpness sharp{.cellEdge=10.26, .nAtoms=2, .alphaMax=40.0};
    EXPECT_GT(UniformMeshCost(sharp), 10*UniformMeshCost(soft))
        << "alpha_max 2 -> 40 is a 20x cutoff, so ~20^1.5 ~ 89x in points";
    const MeshParams becke=BeckeXCParams();
    EXPECT_EQ(BeckeMeshCost(becke,soft), BeckeMeshCost(becke,sharp))
        << "the Becke mesh is atom-centred: its size must not depend on the basis exponents";
}

// The Becke cost is nAtoms x nRadial x nDirs, doubled when imposed (the site-adapted invariant mesh prices
// ~2x the free mesh -- measured 1.97x at degree 29).
TEST(XCPolicy, BeckeCostIsAtomsTimesRadialTimesDirections)
{
    XCMeshSharpness s{.cellEdge=10.26, .nAtoms=2, .alphaMax=2.0};
    MeshParams mp=BeckeXCParams(40, 2.0, 29);
    mp.angular=AngularKind::GaussLegendre;
    EXPECT_EQ(AngularDirections(AngularKind::GaussLegendre,29), 450);   // (L+1)^2/2
    EXPECT_EQ(BeckeMeshCost(mp,s), long(2)*40*450);
    s.imposed=true;
    EXPECT_EQ(BeckeMeshCost(mp,s), long(2)*40*450*2);
    mp.angular=AngularKind::Lebedev;                                    // same degree, cheaper rule
    EXPECT_EQ(AngularDirections(AngularKind::Lebedev,29), 302);
}

// AngularDirections must round UP through the menu's principled gaps (degree 25 has no rule -- 230 carries
// negative weights -- so it costs the degree-29 grid), and must stay QUIET while doing it: it is a cost
// probe on a grid the run may never build, so ResolveLebedev's announcement would be noise.
TEST(XCPolicy, LebedevDirectionsRoundUpThroughTheGaps)
{
    EXPECT_EQ(AngularDirections(AngularKind::Lebedev,23), 194);
    EXPECT_EQ(AngularDirections(AngularKind::Lebedev,25), 302) << "degree 25 is a gap -> the degree-29 rule";
    EXPECT_EQ(AngularDirections(AngularKind::Lebedev,3 ),   6) << "6 and 8 share degree 3; cheapest wins";
}

// The PP is an independent sharpness source, which is the point of carrying alpha_pp separately: a soft
// basis with a semicore PP (Na q9, alpha=15.4) needs a finer uniform grid than the basis alone implies.
TEST(XCPolicy, PPSharpnessRaisesTheUniformRequirement)
{
    XCMeshSharpness noPP{.cellEdge=8.7, .nAtoms=2, .alphaMax=2.0};
    XCMeshSharpness withPP=noPP; withPP.alphaPP=15.43;                  // Na q9
    EXPECT_DOUBLE_EQ(RequiredUniformCutoff(noPP  ), 4.0);
    EXPECT_DOUBLE_EQ(RequiredUniformCutoff(withPP), 4.0+15.43);
    EXPECT_GT(UniformMeshCost(withPP), UniformMeshCost(noPP));
}

//================================================================================================
//  The radial-adequacy statistic (V2.7)
//================================================================================================

// Pin the statistic against the V2.6 ladder's own calibration anchors.  These numbers were quoted in the
// measurement record BEFORE the closed form was coded (Si "3.4 at nR=30", Al "3.2 at nR=30"), so agreement
// here means the code computes the quantity the measurement validated, not a lookalike.
TEST(XCPolicy, RadialResolutionRatioMatchesTheMeasuredAnchors)
{
    auto mp=[](int nR){ return BeckeXCParams(nR, 2.0, 29); };           // the production MHL family (m=2)
    EXPECT_NEAR(RadialResolutionRatio(mp(30),  2.0), 10.0/3.0, 1e-9);   // Si alpha_max=2: the quoted 3.4
    EXPECT_NEAR(RadialResolutionRatio(mp(30),  4.0), 3.13, 0.01);       // Al alpha_max=4: the quoted 3.2
    // Linear in nRadial (x is uniform), so the production nR=40 anchors follow:
    EXPECT_NEAR(RadialResolutionRatio(mp(40),  2.0), 40.0/9.0, 1e-9);   // Si production: 4.44
    EXPECT_NEAR(RadialResolutionRatio(mp(40),  4.0), 4.17, 0.01);       // Al production: the metal floor
    EXPECT_NEAR(RadialResolutionRatio(mp(40), 40.0), 3.09, 0.01);       // NaF production: the tightest
}

// The warning's premise, pinned from both sides: every production config clears the floor (no console noise
// on a healthy run), and the V2.6-measured 4-50x-INADEQUATE class (nR=20 on any of the four systems) sits
// below it.  If either side drifts, the floor needs re-fitting -- see kRadialRatioFloor's calibration note.
TEST(XCPolicy, RadialRatioFloorSeparatesProductionFromTheMeasuredInadequateClass)
{
    auto mp=[](int nR){ return BeckeXCParams(nR, 2.0, 29); };
    for (double aMax : {2.0, 4.0, 40.0})                                // Si, Al, NaF sharpness
    {
        EXPECT_GE(RadialResolutionRatio(mp(40), aMax), kRadialRatioFloor)
            << "production nR=40 must stay quiet at alpha_max=" << aMax;
        EXPECT_LT(RadialResolutionRatio(mp(20), aMax), kRadialRatioFloor)
            << "the measured-inadequate nR=20 rung must warn at alpha_max=" << aMax;
    }
    // Off-MHL the statistic makes NO claim (infinity): the closed form is MHL's spacing, nobody else's.
    MeshParams log_=BeckeXCParams(40, 2.0, 29); log_.radial=RadialKind::Log;
    EXPECT_TRUE(std::isinf(RadialResolutionRatio(log_, 40.0)));
}

//================================================================================================
//  The selector
//================================================================================================

// Explicit choices are HONOURED -- the user's ruling is that this is their decision, so no verdict of the
// cost model may override it.  (The diagnostics are console-only; what is pinned here is the return value.)
TEST(XCPolicy, ExplicitChoiceIsAlwaysHonoured)
{
    XCMeshSharpness sharp{.cellEdge=8.7, .nAtoms=2, .alphaMax=40.0};    // Becke territory
    MeshParams u; u.cellKind=UnitCellKind::Uniform; u.nUniform=20;
    EXPECT_EQ(ResolveXCMesh(u,sharp).cellKind, UnitCellKind::Uniform)   // dominated, but honoured
        << "the selector must never override an explicit grid choice";
    MeshParams b; b.cellKind=UnitCellKind::Becke;
    XCMeshSharpness soft{.cellEdge=10.26, .nAtoms=2, .alphaMax=2.0};
    EXPECT_EQ(ResolveXCMesh(b,soft).cellKind, UnitCellKind::Becke);
}

// Auto with NO sharpness information must fall back to Becke, the SAFE grid -- never guess uniform from a
// default-constructed sharpness, which would be an under-resolution dressed as a decision.
TEST(XCPolicy, AutoWithoutSharpnessTakesTheSafeGrid)
{
    MeshParams a; a.cellKind=UnitCellKind::Auto;
    EXPECT_EQ(ResolveXCMesh(a, XCMeshSharpness{}).cellKind, UnitCellKind::Becke);
    EXPECT_EQ(ResolveXCMesh(a).cellKind,                    UnitCellKind::Becke);
}

// ARMED 2026-08-08 (V2.4): Auto now ACTS on the cost verdict.  The margin was validated by converged-run
// A/B on both systems the selector routes to uniform (GPW_SCF.DISABLED_GridRouteAB_*): uniform at its own
// cutoff matched the fine reference at least as well as the PRODUCTION Becke mesh, for 4-15x fewer points.
// This test replaces AutoStillResolvesToBeckeWhileTheSelectorIsDisarmed, deliberately, as that item required.
TEST(XCPolicy, AutoActsOnTheCostVerdict)
{
    MeshParams a; a.cellKind=UnitCellKind::Auto;
    // A genuinely soft PP-free system.  NB this WAS Si/sipp (alphaPP=2.58, uniform 13,824 pts) -- the
    // R2.15 Lebedev flip cheapened the Becke side 450->302 dirs (36,000->24,160 pts), which moves
    // Si/sipp INSIDE the 2x margin: it now routes to Becke, the safe direction, exactly as the cost
    // model should.  What this test pins is the MECHANISM (the verdict is acted on and sized), so the
    // uniform-side case is one that wins on either scheme.
    XCMeshSharpness soft{.cellEdge=10.26, .nAtoms=2, .alphaMax=2.0};
    EXPECT_LE(double(UniformMeshCost(soft))*kUniformMargin, double(BeckeMeshCost(BeckeXCParams(),soft)));
    const MeshParams u=ResolveXCMesh(a,soft);
    EXPECT_EQ(u.cellKind, UnitCellKind::Uniform);
    // ...and it must SIZE what it chose.  Handing back a bare cellKind would leave nUniform's basis-blind
    // default of 20 in charge, which is the under-resolution the selector exists to prevent.
    EXPECT_DOUBLE_EQ(u.eCut, RequiredUniformCutoff(soft));

    XCMeshSharpness sharp{.cellEdge=8.7, .nAtoms=2, .alphaMax=40.0};                  // F-like: Becke wins
    EXPECT_EQ(ResolveXCMesh(a,sharp).cellKind, UnitCellKind::Becke);
}

// When the selector IS armed it must SIZE the uniform mesh it hands back -- returning a bare cellKind would
// leave nUniform's basis-blind default of 20 in charge, i.e. reintroduce the exact bug this prevents.
// (Verified through the arithmetic rather than the env latch, which is a process-lifetime static.)
TEST(XCPolicy, AnArmedUniformVerdictWouldCarryItsOwnCutoff)
{
    XCMeshSharpness soft{.cellEdge=10.26, .nAtoms=2, .alphaMax=2.0, .alphaPP=2.58};
    const double eReq=RequiredUniformCutoff(soft);
    EXPECT_GT(eReq, UniformCutoff(soft.cellEdge, MeshParams{}.nUniform))
        << "if the default nUniform already sufficed, sizing the verdict would be pointless -- and this "
           "system is the SOFTEST case, so it holds a fortiori for sharper ones";
}
