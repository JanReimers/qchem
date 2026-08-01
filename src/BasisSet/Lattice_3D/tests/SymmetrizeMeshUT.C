// File: SymmetrizeMeshUT.C  Mesh symmetry folding (doc/SymmetryUpgradePlan.md §6 T2 groundwork).
//
// The {r}-quadrature layer of the T2 plan: torus point folding, the invariance checker, the
// group-averaged invariant-mesh constructor, the compact weight fold, and the pointwise
// star-average -- everything the site-adapted Becke grid needs BEFORE any XC wiring.
#include "gtest/gtest.h"
#include <vector>
#include <cmath>

import qchem.BasisSet.Lattice_3D.Internal.SymmetrizeMesh;   // FoldMesh/Symmetrize/UnmatchedCounts/MakeInvariant
import qchem.Types;

using namespace qchem;
using qchem::qcMesh::Mesh;
using namespace qchem::BasisSet::Lattice_3D::Internal;
namespace SL = qchem::Symmetry::Lattice_3D;

// Simple-cubic monatomic crystal: O_h (48 symmorphic ops), identity cell -> fractional == Cartesian.
static SL::SpaceGroup Cubic()
{
    std::vector<SL::AtomSite> basis = {{14, rvec3_t(0,0,0)}};
    return SL::SpaceGroup::Detect(Matrix3D<double>(), basis);
}

// A totally symmetric (O_h-invariant, torus-periodic) test integrand.
static double SymF(const rvec3_t& r)
{
    auto c=[](double x){return std::cos(2*M_PI*x);};
    return c(r.x)+c(r.y)+c(r.z) + c(r.x)*c(r.y)*c(r.z);
}

//---------------------------------------------------------------------------------------
TEST(SymmetrizeMesh, GenericPointIsNotInvariantUntilMadeInvariant)
{
    SL::SpaceGroup sg = Cubic();

    Mesh m;
    m.Append(rvec3_t(0.30, 0.10, 0.00), 2.0);   // a generic-ish point: most O_h ops move it off-mesh

    auto bad = UnmatchedCounts(m, sg);
    ASSERT_EQ(bad.size(), 48u);
    int nBad = 0;
    for (int b : bad) if (b > 0) ++nBad;
    EXPECT_GT(nBad, 30);                        // today's situation: the mesh is NOT op-invariant

    Mesh inv = MakeInvariant(m, sg);
    for (int b : UnmatchedCounts(inv, sg)) EXPECT_EQ(b, 0);   // ...and this is the cure

    // The orbit of (0.30,0.10,0) under O_h on the torus: |O_h|/|stabilizer| = 48/2 = 24 points
    // (the z=0 mirror is its 2-element site group).
    EXPECT_EQ(inv.size(), 24u);

    // Group-averaging conserves the total weight exactly.
    double w0 = 0, w1 = 0;
    for (size_t i = 0; i < m.size();   ++i) w0 += m.Weights()[i];     // indexed: Blaze iterator ops
    for (size_t i = 0; i < inv.size(); ++i) w1 += inv.Weights()[i];   // are not visible cross-module
    EXPECT_NEAR(w1, w0, 1e-14);
}

TEST(SymmetrizeMesh, FoldOfInvariantMeshIsExactForSymmetricIntegrand)
{
    SL::SpaceGroup sg = Cubic();

    // A mildly generic mesh: three seed points with distinct weights, made invariant.
    Mesh seed;
    seed.Append(rvec3_t(0.30, 0.10, 0.00), 2.0);
    seed.Append(rvec3_t(0.25, 0.25, 0.25), 1.0);
    seed.Append(rvec3_t(0.40, 0.00, 0.00), 0.5);
    Mesh inv = MakeInvariant(seed, sg);

    Mesh folded = Symmetrize(inv, sg);
    EXPECT_LT(folded.size(), inv.size()/4);     // large stars fold away (the T2 lever)

    // Sum w conserved; Sum w*f identical for a totally symmetric integrand (~1e-13 reordering tier).
    double w0=0, w1=0, i0=0, i1=0;
    for (size_t i = 0; i < inv.size();    ++i) { w0 += inv.Weights()[i];    i0 += inv.Weights()[i]   *SymF(inv.Points()[i]);    }
    for (size_t i = 0; i < folded.size(); ++i) { w1 += folded.Weights()[i]; i1 += folded.Weights()[i]*SymF(folded.Points()[i]); }
    EXPECT_NEAR(w1, w0, 1e-13);
    EXPECT_NEAR(i1, i0, 1e-13);
}

TEST(SymmetrizeMesh, SummedWeightsFoldExactlyEvenWhenUnequalInStar)
{
    // The fold uses Sum_members w (NOT w_rep*n_star): exact for a symmetric integrand even when the
    // weights are NOT equal within a star (a builder wobble must not become an integration error).
    SL::SpaceGroup sg = Cubic();
    Mesh inv = MakeInvariant([&]{Mesh s; s.Append(rvec3_t(0.30,0.10,0.00),2.0); return s;}(), sg);

    Mesh wobbled;                                        // same points, one weight perturbed
    for (size_t i = 0; i < inv.size(); ++i)
        wobbled.Append(inv.Points()[i], inv.Weights()[i]*(i==3 ? 1.17 : 1.0));

    Mesh folded = Symmetrize(wobbled, sg);
    double i0=0, i1=0;
    for (size_t i = 0; i < wobbled.size(); ++i) i0 += wobbled.Weights()[i]*SymF(wobbled.Points()[i]);
    for (size_t i = 0; i < folded.size();  ++i) i1 += folded.Weights()[i] *SymF(folded.Points()[i]);
    EXPECT_NEAR(i1, i0, 1e-13);
}

TEST(SymmetrizeMesh, PointwiseStarAverageIsTheOrbitMean)
{
    SL::SpaceGroup sg = Cubic();
    Mesh inv = MakeInvariant([&]{Mesh s; s.Append(rvec3_t(0.40,0.00,0.00),1.0); return s;}(), sg);
    ASSERT_EQ(inv.size(), 6u);                           // the axis star {±0.4,0,0} etc.

    SL::Fold f = FoldMesh(inv, sg);
    ASSERT_EQ(f.Reps(), 1u);

    std::vector<double> rho(inv.size(), 1.0);
    rho[2] = 1.6;                                        // a symmetry-broken sample
    SL::SymmetrizeValues(f, rho);
    for (double v : rho) EXPECT_NEAR(v, 1.1, 1e-14);     // the orbit mean (5*1.0 + 1.6)/6

    // Idempotent (a projector on an invariant set under the full group).
    SL::SymmetrizeValues(f, rho);
    for (double v : rho) EXPECT_NEAR(v, 1.1, 1e-14);
}

TEST(SymmetrizeMesh, TorusFoldMatchesAcrossTheCellBoundary)
{
    // {E|(3/4,3/4,3/4)} maps (1/2,1/2,1/2) -> (5/4,..) == (1/4,..) mod 1: only the TORUS fold merges.
    std::vector<rvec3_t> pts = {rvec3_t(0.50,0.50,0.50), rvec3_t(0.25,0.25,0.25)};
    Matrix3D<double> E;
    std::vector<SL::SymOp> ops = {{E, rvec3_t(0,0,0)}, {E, rvec3_t(0.75,0.75,0.75)}};

    SL::Fold torus = SL::FoldPointsPeriodic(pts, ops, 1e-9);
    EXPECT_EQ(torus.Reps(), 1u);
    EXPECT_EQ(torus.starSize[0], 2);

    SL::Fold plain = SL::FoldPoints(pts, ops, 1e-9);     // no wrap: no merge
    EXPECT_EQ(plain.Reps(), 2u);
}
