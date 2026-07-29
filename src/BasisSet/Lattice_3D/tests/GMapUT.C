// File: GMapUT.C  Unit tests for the G-space density symmetrizer (doc/GPWPlan1.md item 3, IBZ).
//
// SymmetrizeGMap replaces rho-tilde(G) by its star average over a reciprocal point group -- the step that
// makes an IBZ-reduced density EXACT (the star weights alone are exact only for the band sum).  These pin the
// core operation: the trivial group is a no-op, charge (the G=0 / total) is conserved, and a cubic point
// group spreads a single G onto its full star with equal weight.

#include "gtest/gtest.h"
#include <vector>
#include <cstdlib>   // std::abs(int)
#include <cmath>     // std::lround

import qchem.BasisSet.Internal.GMap;              // ΔG_Map, SymmetrizeGMap, Matrix3D
import qchem.Symmetry.Lattice_3D.SpaceGroup;      // SpaceGroup::Detect / ReciprocalPointOps (the cubic ops)
import qchem.Types;

using namespace qchem;
using qchem::Symmetry::Lattice_3D::SpaceGroup;
using qchem::Symmetry::Lattice_3D::AtomSite;

// The reciprocal point ops of a simple-cubic monatomic crystal (O_h, 48 ops; centrosymmetric so time
// reversal adds nothing).
static std::vector<Matrix3D<double>> CubicOps()
{
    std::vector<AtomSite> basis = {{14, rvec3_t(0,0,0)}};
    SpaceGroup sg = SpaceGroup::Detect(Matrix3D<double>(), basis);   // identity cell = simple cubic a=1
    return sg.ReciprocalPointOps(/*timeReversal*/true, /*symmorphicOnly*/true);
}

// The reciprocal ops of a PRIMITIVE FCC monatomic crystal -- the actual case the GPW metal tests use.  The
// primitive reciprocal lattice is BCC (non-orthogonal), so these ops are integer matrices that are NOT signed
// axis permutations; the G-index symmetrization must still permute integer G-indices exactly.
static std::vector<Matrix3D<double>> FCCPrimitiveOps()
{
    // FCC primitive cell: columns are ½a(0,1,1),(1,0,1),(1,1,0) with a=1.
    Matrix3D<double> A(0.0,0.5,0.5, 0.5,0.0,0.5, 0.5,0.5,0.0);
    std::vector<AtomSite> basis = {{13, rvec3_t(0,0,0)}};
    SpaceGroup sg = SpaceGroup::Detect(A, basis);
    return sg.ReciprocalPointOps(/*timeReversal*/true, /*symmorphicOnly*/true);
}

// Apply every op to G0, collect the DISTINCT integer images (its star) -- the reconstruction target.
static std::vector<ivec3_t> StarOf(const ivec3_t& G0, const std::vector<Matrix3D<double>>& ops)
{
    std::vector<ivec3_t> star;
    for (const auto& U : ops)
    {
        rvec3_t u = U * rvec3_t(double(G0.x),double(G0.y),double(G0.z));
        ivec3_t m((int)std::lround(u.x),(int)std::lround(u.y),(int)std::lround(u.z));
        bool seen=false; for (auto& s:star) if (s.x==m.x&&s.y==m.y&&s.z==m.z) {seen=true;break;}
        if (!seen) star.push_back(m);
    }
    return star;
}

// THE RECONSTRUCTION IDENTITY (the IBZ correctness core): a single irreducible rep carrying its STAR-SIZE
// weight symmetrizes to its full star at UNIT weight -- i.e. P({G0: |star|}) == {each star member: 1}.  This
// is exactly how P(ρ_reduced) must reconstruct ρ_full: the star-weighted rep spreads back to unit-weighted
// partners.  Tested on the real FCC primitive ops (the metal case), where U·m is a non-trivial integer map.
TEST(SymmetrizeGMap, FCCReconstructsStarFromWeightedRep)
{
    auto ops = FCCPrimitiveOps();
    ASSERT_EQ(ops.size(), 48u);                       // FCC point group is still O_h
    const ivec3_t G0(1,0,0);
    auto star = StarOf(G0, ops);
    ΔG_Map rep = {{G0, dcmplx((double)star.size(), 0.0)}};   // the rep carries its star-size weight
    ΔG_Map out = SymmetrizeGMap(rep, ops);

    EXPECT_EQ(out.size(), star.size()) << "reconstruction must populate exactly the star";
    for (const auto& m : star)
    {
        ASSERT_TRUE(out.count(m)) << "star member missing from reconstruction: "<<m.x<<","<<m.y<<","<<m.z;
        EXPECT_NEAR(out.at(m).real(), 1.0, 1e-12) << "star member should reconstruct at UNIT weight";
        EXPECT_NEAR(out.at(m).imag(), 0.0, 1e-12);
    }
}

// The trivial group (empty ops = {E}) leaves the map untouched -- the "no symmetry is the trivial instance"
// invariant (molecules / Γ / unreduced crystals pass through unchanged).
TEST(SymmetrizeGMap, TrivialGroupIsNoOp)
{
    ΔG_Map rg = {{ivec3_t(1,0,0), dcmplx(2.0,-1.0)}, {ivec3_t(0,0,0), dcmplx(3.0,0.0)}};
    ΔG_Map out = SymmetrizeGMap(rg, {});
    ASSERT_EQ(out.size(), rg.size());
    for (const auto& [m,v] : rg)
    {
        EXPECT_NEAR(out.at(m).real(), v.real(), 1e-15);
        EXPECT_NEAR(out.at(m).imag(), v.imag(), 1e-15);
    }
}

// A single off-centre G is spread over its full O_h star (the six {±1,0,0} points) with equal weight 1/6,
// and the total (the integral / G=0-preserving sum) is conserved.
TEST(SymmetrizeGMap, CubicSpreadsGOverItsStar)
{
    auto ops = CubicOps();
    ASSERT_EQ(ops.size(), 48u);                       // O_h
    ΔG_Map rg = {{ivec3_t(1,0,0), dcmplx(1.0,0.0)}};
    ΔG_Map out = SymmetrizeGMap(rg, ops);

    EXPECT_EQ(out.size(), 6u);                         // the star of (1,0,0): {±1,0,0},{0,±1,0},{0,0,±1}
    for (const auto& [m,v] : out)
    {
        EXPECT_EQ(std::abs(m.x)+std::abs(m.y)+std::abs(m.z), 1);   // a star member
        EXPECT_NEAR(v.real(), 1.0/6.0, 1e-12);
        EXPECT_NEAR(v.imag(), 0.0,     1e-12);
    }
    dcmplx total(0.0); for (const auto& [m,v] : out) total += v;
    EXPECT_NEAR(total.real(), 1.0, 1e-12);            // conserved
}

// The G=0 component is invariant (its own star), and symmetrization conserves the total for a generic field.
TEST(SymmetrizeGMap, ConservesTotalAndFixesGamma)
{
    auto ops = CubicOps();
    ΔG_Map rg = {{ivec3_t(0,0,0), dcmplx(5.0,0.0)}, {ivec3_t(2,1,0), dcmplx(1.0,0.5)},
                 {ivec3_t(1,1,1), dcmplx(-2.0,0.0)}};
    dcmplx before(0.0); for (const auto& [m,v] : rg) before += v;
    ΔG_Map out = SymmetrizeGMap(rg, ops);
    dcmplx after(0.0);  for (const auto& [m,v] : out) after += v;
    EXPECT_NEAR(after.real(), before.real(), 1e-12);
    EXPECT_NEAR(after.imag(), before.imag(), 1e-12);
    EXPECT_NEAR(out.at(ivec3_t(0,0,0)).real(), 5.0, 1e-12);   // Γ is its own star -- unchanged
}
