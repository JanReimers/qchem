// File: SymmetrizeMeshUT.C  Mesh symmetry folding (doc/SymmetryUpgradePlan.md §6 T2 groundwork).
//
// The {r}-quadrature layer of the T2 plan: torus point folding, the invariance checker, the
// group-averaged invariant-mesh constructor, the compact weight fold, and the pointwise
// star-average -- everything the site-adapted Becke grid needs BEFORE any XC wiring.
#include "gtest/gtest.h"
#include <vector>
#include <cmath>

import qchem.SymmetrizeMesh;   // FoldMesh/Symmetrize/UnmatchedCounts/MakeInvariant/MakeInvariantAngularMesh
import qchem.Lattice_3D;       // Lattice_3D::GetSpaceGroup + FCCUnitCell/UnitCell (the site-adapted Becke gate)
import qchem.Mesh.Angular;     // MakeAngular (the GL direction-count baseline for the growth measure)
import qchem.Types;

using namespace qchem;
using qchem::qcMesh::Mesh;
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

//---------------------------------------------------------------------------------------
//  MakeInvariantAngularMesh (§6a W2): the site-group-invariant angular quadrature.
//
// The full O_h as CARTESIAN orthogonal ops: the 48 signed permutation matrices.
static std::vector<Matrix3D<double>> OhOps()
{
    std::vector<Matrix3D<double>> ops;
    int perm[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    for (auto& pr : perm)
        for (int s = 0; s < 8; ++s)
        {
            Matrix3D<double> R(0,0,0, 0,0,0, 0,0,0);
            for (int r = 0; r < 3; ++r) R(r+1, pr[r]+1) = (s>>r & 1) ? -1.0 : 1.0;
            ops.push_back(R);
        }
    return ops;
}

//! Closed-form sphere integral of x^a y^b z^c: 0 if any exponent odd, else
//! 4π (a-1)!!(b-1)!!(c-1)!!/(a+b+c+1)!!.
static double MonomialIntegral(int a, int b, int c)
{
    if (a%2 || b%2 || c%2) return 0.0;
    auto dfac=[](int n){ double f=1; for (int k=n; k>1; k-=2) f*=k; return f; };
    return 4.0*M_PI*dfac(a-1)*dfac(b-1)*dfac(c-1)/dfac(a+b+c+1);
}

static void CheckAngularQuadrature(const Mesh& q, const std::vector<Matrix3D<double>>& ops, int L)
{
    // Positive weights summing to 4π; every point on the unit sphere.
    double wsum = 0;
    for (size_t i = 0; i < q.size(); ++i)
    {
        EXPECT_GT(q.Weights()[i], 0.0);
        wsum += q.Weights()[i];
        const rvec3_t& p = q.Points()[i];
        EXPECT_NEAR(p.x*p.x+p.y*p.y+p.z*p.z, 1.0, 1e-12);
    }
    EXPECT_NEAR(wsum, 4.0*M_PI, 1e-8);

    // INVARIANCE by construction: every op maps every point onto a point (the T2 precondition).
    for (const auto& R : ops)
        for (size_t i = 0; i < q.size(); ++i)
        {
            rvec3_t im = R*q.Points()[i];
            bool found = false;
            for (size_t j = 0; j < q.size(); ++j)
            {
                const rvec3_t& p = q.Points()[j];
                double dx=p.x-im.x, dy=p.y-im.y, dz=p.z-im.z;
                if (dx*dx+dy*dy+dz*dz < 1e-16) {found = true; break;}
            }
            EXPECT_TRUE(found);
        }

    // EXACTNESS to degree L: every monomial x^a y^b z^c, a+b+c <= L, against the closed form.
    for (int a = 0; a <= L; ++a)
        for (int b = 0; a+b <= L; ++b)
            for (int c = 0; a+b+c <= L; ++c)
            {
                double Q = 0;
                for (size_t i = 0; i < q.size(); ++i)
                {
                    const rvec3_t& p = q.Points()[i];
                    Q += q.Weights()[i]*std::pow(p.x,a)*std::pow(p.y,b)*std::pow(p.z,c);
                }
                EXPECT_NEAR(Q, MonomialIntegral(a,b,c), 1e-8)
                    << "monomial ("<<a<<","<<b<<","<<c<<")";
            }
}

TEST(InvariantAngularMesh, OhInvariantAndExactToDegree9)
{
    auto ops = OhOps();
    ASSERT_EQ(ops.size(), 48u);
    Mesh q = MakeInvariantAngularMesh(ops, 9);
    ASSERT_GT(q.size(), 0u);
    CheckAngularQuadrature(q, ops, 9);
    // Deterministic: a second build is identical.
    Mesh q2 = MakeInvariantAngularMesh(ops, 9);
    ASSERT_EQ(q2.size(), q.size());
    EXPECT_NEAR(q2.Points()[0].x, q.Points()[0].x, 0.0);
}

TEST(InvariantAngularMesh, TrivialGroupStillExact)
{
    std::vector<Matrix3D<double>> onlyE = {Matrix3D<double>()};
    Mesh q = MakeInvariantAngularMesh(onlyE, 5);
    ASSERT_GT(q.size(), 0u);
    CheckAngularQuadrature(q, onlyE, 5);
}

//---------------------------------------------------------------------------------------
//  the ops-taking CreateIntegrationMesh (§6a W2b): the whole periodic Becke mesh invariant BY CONSTRUCTION.
//
TEST(SiteAdaptedBeckeMesh, DiamondInvariantByConstruction)
{
    FCCUnitCell cell(10.26);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});          // non-symmorphic: partner grids need the glide's rotation
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    const SL::SpaceGroup& sg = lat.GetSpaceGroup();
    ASSERT_EQ(sg.Order(), 48u);

    std::vector<SL::SymOp> ops;
    for (const auto& op : sg.Ops()) ops.push_back({op.W, op.tau});

    qcMesh::MeshParams mp;
    mp.cellKind=qcMesh::UnitCellKind::Becke;
    mp.radial=qcMesh::RadialKind::MHL; mp.nRadial=6; mp.mhl_m=2; mp.mhl_alpha=2.0;
    mp.angular=qcMesh::AngularKind::GaussLegendre; mp.angularDegree=5;

    Mesh adapted = cell.CreateIntegrationMesh(mp, ops);
    ASSERT_GT(adapted.size(), 0u);

    // THE gate: every op maps every mesh point onto a mesh point -- zero unmatched (the T2
    // precondition today's shared-orientation grids fail).
    for (int bad : UnmatchedCounts(adapted, sg))
        EXPECT_EQ(bad, 0);
    int nBadPlain = 0;                            // ...and the plain mesh really does fail it
    Mesh plain = cell.CreateIntegrationMesh(mp);
    for (int bad : UnmatchedCounts(plain, sg)) if (bad > 0) ++nBadPlain;
    EXPECT_GT(nBadPlain, 30);

    // Same quadrature class: both integrate 1 to the cell volume within radial-tail accuracy.
    double wa=0, wp=0;
    for (size_t i = 0; i < adapted.size(); ++i) wa += adapted.Weights()[i];
    for (size_t i = 0; i < plain.size();   ++i) wp += plain.Weights()[i];
    EXPECT_NEAR(wa/wp, 1.0, 0.05);
}

// §7 step 5 REMAINING: the PRODUCTION-L growth re-measure (plan §6a W2b close-out).  The W2b
// measurement at coarse L=9/T_d was 240 dirs vs GL-9's 50 (the low-L generic-orbit premium);
// the claim to test is the ~2x asymptote at the production recipe L=29 (GL-29 = 450 dirs).
// This is a MEASUREMENT gate: it prints the numbers and pins only the asymptote class.
TEST(InvariantAngularMesh, ProductionL29Growth)
{
    FCCUnitCell cell(10.26);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    const SL::SpaceGroup& sg = lat.GetSpaceGroup();

    // Diamond's T_d site group as Cartesian ops (A W A^-1), exactly as the ops-taking CreateIntegrationMesh does.
    const Matrix3D<double>& A = cell.GetCellMatrix();
    const Matrix3D<double> Ainv = Invert(A);
    std::vector<Matrix3D<double>> siteOps;
    for (const auto& op : sg.SiteStabilizer(rvec3_t(0,0,0))) siteOps.push_back(A*op.W*Ainv);
    ASSERT_EQ(siteOps.size(), 24u);

    const int L = 29;
    Mesh q = MakeInvariantAngularMesh(siteOps, L);
    qcMesh::MeshParams mp;
    mp.angular=qcMesh::AngularKind::GaussLegendre; mp.angularDegree=L;
    const size_t nGL = qcMesh::MakeAngular(mp).size();

    const double growth = double(q.size())/double(nGL);
    std::cout << "[MEASURE] L=" << L << " T_d invariant dirs=" << q.size()
              << " GL dirs=" << nGL << " growth=" << growth << "x" << std::endl;

    // Invariant by construction, positive weights, Sum w = 4pi (the cheap subset of
    // CheckAngularQuadrature -- the full monomial sweep at L=29 is the expensive part).
    double wsum = 0;
    for (size_t i = 0; i < q.size(); ++i) {EXPECT_GT(q.Weights()[i], 0.0); wsum += q.Weights()[i];}
    EXPECT_NEAR(wsum, 4.0*M_PI, 1e-8);

    // The asymptote claim: the production-L premium must be well under the coarse-L 4.8x
    // (240/50), in the ~2x class.
    EXPECT_LT(growth, 2.5);
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

// The invariance contract AT SCALE (production angular recipe L=29): with tens of thousands of
// points, eps-BORDERLINE tail-drop decisions in the Becke builder flip inconsistently across orbit
// partners (bit-different rotated distances) -- measured 212 orphaned points at nR=40/L=29 before
// the orbit-consistency filter in the ops-taking CreateIntegrationMesh.  This gate holds the cured contract:
// UnmatchedCounts identically zero at a recipe big enough to exercise the borderline channel.
// (The full production growth measure, one-off 2026-08-02: adapted 49172 vs plain 25007 = 1.97x,
// the dirs-ratio asymptote of the ProductionL29Growth angular gate above.)
TEST(SiteAdaptedBeckeMesh, ProductionAngularRecipeInvariantAfterTailDrop)
{
    FCCUnitCell cell(10.26);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    const SL::SpaceGroup& sg = lat.GetSpaceGroup();
    std::vector<SL::SymOp> ops;
    for (const auto& op : sg.Ops()) ops.push_back({op.W, op.tau});

    qcMesh::MeshParams mp;
    mp.cellKind=qcMesh::UnitCellKind::Becke;
    mp.radial=qcMesh::RadialKind::MHL; mp.nRadial=15; mp.mhl_m=2; mp.mhl_alpha=2.0;
    mp.angular=qcMesh::AngularKind::GaussLegendre; mp.angularDegree=29;

    Mesh adapted = cell.CreateIntegrationMesh(mp, ops);
    ASSERT_GT(adapted.size(), 10000u);
    for (int bad : UnmatchedCounts(adapted, sg)) EXPECT_EQ(bad, 0);
}
