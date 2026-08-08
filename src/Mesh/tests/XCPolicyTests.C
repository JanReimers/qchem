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

// DISARMED: until kUniformMargin is calibrated against a D8 grid-convergence measurement, Auto keeps
// returning the Becke recipe whatever the cost verdict -- arming it would move every pinned anchor on a
// heuristic.  This test is the guard on that promise; it must be UPDATED, deliberately, when the selector
// is armed.
TEST(XCPolicy, AutoStillResolvesToBeckeWhileTheSelectorIsDisarmed)
{
    MeshParams a; a.cellKind=UnitCellKind::Auto;
    XCMeshSharpness soft{.cellEdge=10.26, .nAtoms=2, .alphaMax=2.0, .alphaPP=2.58};   // Si/sipp
    // The verdict itself says uniform -- which is precisely why the disarm matters.
    EXPECT_LE(double(UniformMeshCost(soft))*kUniformMargin, double(BeckeMeshCost(BeckeXCParams(),soft)));
    EXPECT_EQ(ResolveXCMesh(a,soft).cellKind, UnitCellKind::Becke)
        << "Auto must not act on the cost verdict until the margin is measured";
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
