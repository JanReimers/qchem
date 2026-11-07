// File: BeckeClimbTests.C  The one-round Becke climb (qchem.Mesh.BeckeClimb) on a synthetic error model: a pure function, no Hamiltonian.
//
// The model: a rung's V_xc error is an angular part a(L) plus a radial part r(nR), each falling geometrically along its axis, so the rung
// difference on one axis is exactly that axis's error minus the next rung's -- and the geometric tail estimate is EXACT.  The climb must
// then stop at the cheapest rung that meets the target on both axes, never earlier.
#include "gtest/gtest.h"
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>
import qchem.Mesh.BeckeClimb;

using namespace qchem;
using qcMesh::MeshParams;

namespace
{
const int kDeg[]={29,35,41,47,53,65,77,89};
const int kRad[]={40,60,80,100,120,140};
int AngIndex(int L) {for (int i=0;i<8;++i) if (kDeg[i]==L) return i; ADD_FAILURE()<<"degree "<<L<<" is not a ladder rung"; return 0;}
int RadIndex(int n) {for (int i=0;i<6;++i) if (kRad[i]==n) return i; ADD_FAILURE()<<"count "<<n<<" is not a ladder rung"; return 0;}

//! The synthetic surface: error(L, nR) = A f^iL + B g^iR (both parts positive, so a difference is the parts' difference).
struct Model
{
    double A, f, B, g;
    double Ang(int L) const {return A*std::pow(f, AngIndex(L));}
    double Rad(int n) const {return B*std::pow(g, RadIndex(n));}
    double Error(const MeshParams& m) const {return Ang(m.angularDegree)+Rad(m.nRadial);}
    //! max|ΔH| between two rungs: |error(a) - error(b)|, and the climb asks only for one-axis neighbours
    double Distance(const MeshParams& a, const MeshParams& b) const {return std::fabs(Error(a)-Error(b));}
};

MeshParams Rung(int nR, int L)
{
    MeshParams m=qcMesh::BeckeXCParams();   // the production recipe (a Becke cellKind is recipe-stated)
    m.nRadial=nR;
    m.angularDegree=L;
    m.angular=qcMesh::AngularKind::Lebedev;
    return m;
}
}

TEST(BeckeClimb, TheTailEstimateIsTheGeometricTailAndFallsBackToTheDifference)
{
    EXPECT_DOUBLE_EQ(qcMesh::BeckeTailEstimate(1e-5, std::nullopt), 1e-5) << "first measurement on an axis: no ratio yet";
    EXPECT_NEAR(qcMesh::BeckeTailEstimate(1e-5, 4e-5), 1e-5/(1.0-0.25), 1e-20) << "r = 1/4: the tail is d/(1 - r)";
    EXPECT_DOUBLE_EQ(qcMesh::BeckeTailEstimate(4e-5, 1e-5), 4e-5) << "r >= 1 (not converging): the difference itself";
    EXPECT_GT(qcMesh::BeckeTailEstimate(2.9e-5, 3.0e-5), 1e-4) << "r near 1 forces a raise, never an early stop";
}

TEST(BeckeClimb, ItStopsWhereBothAxesEstimatesMeetTheTarget)
{
    const qcMesh::BeckeLadderCalibration cal;
    const Model m{.A=3e-4, .f=0.4, .B=2e-4, .g=0.2};
    const double target=1e-5;
    const auto R=qcMesh::ClimbBeckeLadder(Rung(40,29), target, cal, false, [&](const MeshParams& a, const MeshParams& b){return m.Distance(a,b);}, 0);
    ASSERT_EQ(R.end, qcMesh::BeckeClimbResult::End::Settled) << R.why;
    ASSERT_GE(R.steps.size(), 2u);
    for (size_t k=R.steps.size()-2;k<R.steps.size();++k)
    {
        EXPECT_FALSE(R.steps[k].over) << R.steps[k].axis;
        EXPECT_TRUE(R.steps[k].rung.nRadial==R.mesh.nRadial && R.steps[k].rung.angularDegree==R.mesh.angularDegree) << "the last round audits the stop";
    }
    // What the rule GUARANTEES on a geometric surface: each axis's own error within target/(1 - r) -- the estimate is the exact tail when it has a
    // ratio, and the bare difference (an under-read by 1 - r) when it has none (see the next test).
    EXPECT_LE(m.Ang(R.mesh.angularDegree), target/(1.0-m.f));
    EXPECT_LE(m.Rad(R.mesh.nRadial), target/(1.0-m.g));
    EXPECT_EQ(R.raises, AngIndex(R.mesh.angularDegree)+RadIndex(R.mesh.nRadial)) << "every raise is one rung on one axis";
    EXPECT_EQ(R.steps.size(), size_t(2*(R.raises+1))) << "two audits per rung visited, each axis against its rung above";
    for (const auto& s : R.steps) EXPECT_EQ(s.over, s.tailEstimate>target);
}

TEST(BeckeClimb, WithoutARatioTheEstimateUnderReadsTheTailByOneMinusR)
{
    // The pilot's rule, kept as it is (and replayed on 6 banked cases, all meeting eps): an axis's FIRST measurement, or one after the
    // OTHER axis moved (no same-axis neighbour below), has no ratio, so the estimate is the bare difference d = error x (1 - r).  Here the
    // angular error at the start is 1.5x the target and d reads 0.9x: the climb stops.  Bounded by 1/(1 - r); a named fallback ratio from
    // BeckeLadderCalibration would close it (OpenWork, the seed row's follow-up).
    const qcMesh::BeckeLadderCalibration cal;
    const double target=1e-5;
    const Model m{.A=1.5e-5, .f=0.4, .B=1e-12, .g=0.1};
    const auto R=qcMesh::ClimbBeckeLadder(Rung(80,29), target, cal, false, [&](const MeshParams& a, const MeshParams& b){return m.Distance(a,b);}, 0);
    EXPECT_EQ(R.end, qcMesh::BeckeClimbResult::End::Settled);
    EXPECT_EQ(R.raises, 0);
    EXPECT_GT(m.Error(R.mesh), target) << "the known under-read";
    EXPECT_LE(m.Ang(R.mesh.angularDegree), target/(1.0-m.f));
}

TEST(BeckeClimb, OneAxisAloneWouldStopEarlyBothAxesDoNot)
{
    // The per-axis caveat (doc/Records/SeedTestReport.md addendum, LiF nR 30): the angular part is converged, the radial part is not.  The
    // angular difference reads tiny; only the radial audit sees the error.
    const qcMesh::BeckeLadderCalibration cal;
    const Model m{.A=1e-9, .f=0.5, .B=2e-3, .g=0.1};
    const auto R=qcMesh::ClimbBeckeLadder(Rung(40,53), 1e-5, cal, false, [&](const MeshParams& a, const MeshParams& b){return m.Distance(a,b);}, 0);
    ASSERT_EQ(R.end, qcMesh::BeckeClimbResult::End::Settled);
    ASSERT_GE(R.steps.size(), 2u);
    EXPECT_EQ(R.steps[0].axis, "angular");
    EXPECT_FALSE(R.steps[0].over) << "the angular audit alone would accept the starting rung";
    EXPECT_TRUE(R.steps[1].over) << "the radial audit refuses it";
    EXPECT_GT(R.mesh.nRadial, 40);
    EXPECT_LE(m.Error(R.mesh), 1e-5);
}

TEST(BeckeClimb, AfterTheOtherAxisMovesARatioNearOneRaisesAConvergedAxis)
{
    // The pilot's rule, kept as it is: an axis's ratio is taken against its last measurement at ANY other rung.  When the OTHER axis moved,
    // this axis's difference is unchanged up to noise, so r ~ 1 and d/(1 - r) explodes: a converged axis is raised.  A SAFE-side cost (extra
    // rungs), never an early stop.  The obvious repair -- a ratio only between same-axis neighbours, the calibrated factor otherwise --
    // FAILS on banked data (NaF eps 1e-4: stops at (40, 29), 1.16 eps; its low-L angular tail falls at r ~ 0.74, not the calibrated 0.4),
    // so it is not made here (doc/Records/SeedTestReport.md addendum; OpenWork, the seed row).
    const qcMesh::BeckeLadderCalibration cal;
    const double target=1e-5;
    const Model m{.A=1e-9, .f=0.5, .B=2e-3, .g=0.1};
    double noise=1.0;   // the angular differences at successive radial rungs, each a hair smaller (rounding, in a real run)
    const auto R=qcMesh::ClimbBeckeLadder(Rung(40,53), target, cal, false, [&](const MeshParams& a, const MeshParams& b)
    {
        const double d=m.Distance(a,b);
        if (a.nRadial==b.nRadial) { noise*=1.0-1e-9; return d*noise; }
        return d;
    }, 0);
    // With a difference that keeps shrinking by a hair, every re-measure reads r ~ 1, so the artefact walks the converged axis to its TOP:
    // the climb ends unsettled at degree 89 (honestly reported, never an early stop) -- the cost hazard the OpenWork row names.
    EXPECT_EQ(R.end, qcMesh::BeckeClimbResult::End::AxisEnd);
    EXPECT_EQ(R.mesh.angularDegree, 89) << "the converged angular axis climbed to its end";
    EXPECT_LE(m.Error(R.mesh), target) << "the stop is safe";
}

TEST(BeckeClimb, AConfirmingRoundRaisesAtMostOnceAndSaysSo)
{
    const qcMesh::BeckeLadderCalibration cal;
    const Model m{.A=3e-4, .f=0.4, .B=2e-4, .g=0.2};
    auto dist=[&](const MeshParams& a, const MeshParams& b){return m.Distance(a,b);};
    const auto settled=qcMesh::ClimbBeckeLadder(Rung(40,29), 1e-5, cal, false, dist, 0);
    const auto confirm=qcMesh::ClimbBeckeLadder(settled.mesh, 1e-5, cal, false, dist, 1);
    EXPECT_EQ(confirm.end, qcMesh::BeckeClimbResult::End::Settled) << "the settled rung confirms at the same D";
    EXPECT_EQ(confirm.raises, 0);
    EXPECT_EQ(confirm.steps.size(), 2u) << "one audit per axis, no more";
    // At a D where the settled rung is over (the confirming audit DISAGREES), the round raises once and stops, unaudited -- never loops.
    const auto disagree=qcMesh::ClimbBeckeLadder(Rung(40,29), 1e-5, cal, false, dist, 1);
    EXPECT_EQ(disagree.end, qcMesh::BeckeClimbResult::End::RaiseCap);
    EXPECT_EQ(disagree.raises, 1);
    EXPECT_EQ(disagree.steps.size(), 2u);
    EXPECT_TRUE(disagree.mesh.nRadial>40 || disagree.mesh.angularDegree>29) << "the raise is one rung on one axis";
}

TEST(BeckeClimb, AnAxisAtItsEndIsUnauditedTheOtherStillClimbs)
{
    const qcMesh::BeckeLadderCalibration cal;
    const Model m{.A=1e-9, .f=0.5, .B=2e-4, .g=0.2};
    std::map<std::pair<int,int>, int> calls;
    const auto R=qcMesh::ClimbBeckeLadder(Rung(40,89), 1e-5, cal, false, [&](const MeshParams& a, const MeshParams& b)
    {
        ++calls[{a.nRadial*1000+a.angularDegree, b.nRadial*1000+b.angularDegree}];
        return m.Distance(a,b);
    }, 0);
    EXPECT_EQ(R.end, qcMesh::BeckeClimbResult::End::AxisEnd) << "degree 89 is the angular axis's top";
    EXPECT_NE(R.why.find("angular"), std::string::npos) << R.why;
    EXPECT_GT(R.mesh.nRadial, 40) << "the radial axis climbed regardless";
    EXPECT_LE(m.Rad(R.mesh.nRadial), 1e-5);
    for (const auto& s : R.steps) EXPECT_EQ(s.axis, "radial");
    for (const auto& [k, n] : calls) EXPECT_EQ(n, 1) << "each pair is measured once";
}

TEST(BeckeClimb, AStateOfNoRungsLeftOnEitherAxisIsReportedNotLooped)
{
    const qcMesh::BeckeLadderCalibration cal;
    const Model m{.A=1e-3, .f=0.99, .B=1e-3, .g=0.99};   // never meets 1e-9 on this ladder
    const auto R=qcMesh::ClimbBeckeLadder(Rung(40,29), 1e-9, cal, false, [&](const MeshParams& a, const MeshParams& b){return m.Distance(a,b);}, 0);
    EXPECT_EQ(R.end, qcMesh::BeckeClimbResult::End::AxisEnd);
    EXPECT_EQ(R.mesh.angularDegree, 89);
    EXPECT_LE(R.raises, 7+5);
}

TEST(BeckeClimb, ABadArgumentThrows)
{
    const qcMesh::BeckeLadderCalibration cal;
    auto zero=[](const MeshParams&, const MeshParams&){return 0.0;};
    EXPECT_THROW(qcMesh::ClimbBeckeLadder(Rung(40,29), 1e-5, cal, false, zero, -1), std::invalid_argument);
    EXPECT_THROW(qcMesh::ClimbBeckeLadder(Rung(40,29), 0.0,  cal, false, zero, 0),  std::invalid_argument);
}
