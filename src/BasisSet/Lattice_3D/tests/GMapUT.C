// File: GMapUT.C  Unit tests for the G-space density symmetrizer (doc/GPWPlan1.md items 3 + 5, IBZ).
//
// SymmetrizeGMap replaces rho-tilde(G) by its star average over the crystal point group -- the step that
// makes an IBZ-reduced density EXACT (the star weights alone are exact only for the band sum).  These pin the
// core operation: the trivial group is a no-op, charge (the G=0 / total) is conserved, a cubic point group
// spreads a single G onto its full star, and -- the NON-SYMMORPHIC case -- a diamond structure factor stays
// invariant under the glide ops only when the e^{+2pi i m.tau} phase is carried (the whole point of item 5).

#include "gtest/gtest.h"
#include <vector>
#include <cstdlib>   // std::abs(int)
#include <cmath>     // std::lround, std::acos
#include <complex>   // std::polar

import qchem.BasisSet.Internal.GMap;              // ΔG_Map, SymmetrizeGMap, ReciprocalOp
import qchem.Symmetry.Lattice_3D.SpaceGroup;      // SpaceGroup::Detect / ReciprocalOps (the crystal ops)
import qchem.Types;

using namespace qchem;
using qchem::Symmetry::Lattice_3D::SpaceGroup;
using qchem::Symmetry::Lattice_3D::AtomSite;
using qchem::Symmetry::Lattice_3D::ReciprocalOp;

static const double PI = std::acos(-1.0);

// The reciprocal {U|τ} ops of a simple-cubic monatomic crystal (O_h, 48 ops; τ=0 -- symmorphic).
static std::vector<ReciprocalOp> CubicOps()
{
    std::vector<AtomSite> basis = {{14, rvec3_t(0,0,0)}};
    SpaceGroup sg = SpaceGroup::Detect(Matrix3D<double>(), basis);   // identity cell = simple cubic a=1
    return sg.ReciprocalOps();
}

// The reciprocal ops of a PRIMITIVE FCC monatomic crystal -- the actual case the GPW metal tests use.  The
// primitive reciprocal lattice is BCC (non-orthogonal), so these ops are integer matrices that are NOT signed
// axis permutations; the G-index symmetrization must still permute integer G-indices exactly.
static std::vector<ReciprocalOp> FCCPrimitiveOps()
{
    Matrix3D<double> A(0.0,0.5,0.5, 0.5,0.0,0.5, 0.5,0.5,0.0);   // ½a(0,1,1),(1,0,1),(1,1,0), a=1
    std::vector<AtomSite> basis = {{13, rvec3_t(0,0,0)}};
    return SpaceGroup::Detect(A, basis).ReciprocalOps();
}

// The reciprocal {U|τ} ops of DIAMOND (Fd-3m, NON-symmorphic): FCC lattice + 2-atom basis at (0,0,0),(¼,¼,¼).
static std::vector<ReciprocalOp> DiamondOps()
{
    Matrix3D<double> A(0.0,0.5,0.5, 0.5,0.0,0.5, 0.5,0.5,0.0);
    std::vector<AtomSite> basis = {{14, rvec3_t(0,0,0)}, {14, rvec3_t(0.25,0.25,0.25)}};
    return SpaceGroup::Detect(A, basis).ReciprocalOps();
}

// Apply the scatter matrix U (=Wᵀ) of every op to G0, collect the DISTINCT integer images (its star).
static std::vector<ivec3_t> StarOf(const ivec3_t& G0, const std::vector<ReciprocalOp>& ops)
{
    std::vector<ivec3_t> star;
    for (const auto& op : ops)
    {
        rvec3_t u = op.U * rvec3_t(double(G0.x),double(G0.y),double(G0.z));
        ivec3_t m((int)std::lround(u.x),(int)std::lround(u.y),(int)std::lround(u.z));
        bool seen=false; for (auto& s:star) if (s.x==m.x&&s.y==m.y&&s.z==m.z) {seen=true;break;}
        if (!seen) star.push_back(m);
    }
    return star;
}

// The closure of a seed set under the scatter matrices (a G-set closed under the ops -- exact for a projector).
static std::vector<ivec3_t> Orbit(const std::vector<ivec3_t>& seeds, const std::vector<ReciprocalOp>& ops)
{
    std::vector<ivec3_t> orb = seeds;
    auto has=[&](const ivec3_t& m){ for (auto& s:orb) if (s.x==m.x&&s.y==m.y&&s.z==m.z) return true; return false; };
    for (size_t q=0; q<orb.size(); ++q)
        for (const auto& op : ops)
        {
            rvec3_t u = op.U * rvec3_t(double(orb[q].x),double(orb[q].y),double(orb[q].z));
            ivec3_t m((int)std::lround(u.x),(int)std::lround(u.y),(int)std::lround(u.z));
            if (!has(m)) orb.push_back(m);
        }
    return orb;
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

// NON-SYMMORPHIC (item 5): the diamond structure factor ρ̃(m)=Σ_atoms e^{-2πi m·R} = 1 + e^{-2πi m·τ0}
// (τ0=(¼,¼,¼)) is EXACTLY Fd-3m invariant.  On a G-set closed under the ops, a CORRECT non-symmorphic
// symmetrizer (carrying the e^{+2πi m·τ} glide phase) leaves it UNCHANGED; dropping the phase (the old W-only
// guard) or getting its sign wrong corrupts it.  This pins the glide phase decisively without an SCF.
TEST(SymmetrizeGMap, NonSymmorphicDiamondStructureFactorInvariant)
{
    auto ops = DiamondOps();
    ASSERT_EQ(ops.size(), 48u);                        // full O_h point part of Fd-3m

    const rvec3_t tau0(0.25,0.25,0.25);
    auto sf=[&](const ivec3_t& m){ double p=-2.0*PI*(m.x*tau0.x+m.y*tau0.y+m.z*tau0.z);
                                   return 1.0 + std::polar(1.0,p); };   // 1 + e^{-2πi m·τ0}

    // A G-set closed under the ops (union of a couple of orbits) so the projector is exact on it.
    auto S = Orbit({ivec3_t(0,0,0), ivec3_t(1,0,0), ivec3_t(1,1,1), ivec3_t(2,0,0)}, ops);
    ASSERT_GT(S.size(), 10u);
    ΔG_Map rg; for (const auto& m : S) rg[m] = sf(m);

    ΔG_Map out = SymmetrizeGMap(rg, ops);
    for (const auto& m : S)
    {
        ASSERT_TRUE(out.count(m));
        EXPECT_NEAR(out.at(m).real(), rg.at(m).real(), 1e-10)
            << "glide phase must leave the diamond structure factor invariant at "<<m.x<<","<<m.y<<","<<m.z;
        EXPECT_NEAR(out.at(m).imag(), rg.at(m).imag(), 1e-10);
    }
}

// CONTROL: without the τ phase (τ forced to 0, the old W-only guard) the SAME diamond structure factor is NOT
// invariant -- proving the phase is doing real work (the test above is not vacuously passing).
TEST(SymmetrizeGMap, NonSymmorphicWithoutPhaseCorruptsStructureFactor)
{
    auto ops = DiamondOps();
    std::vector<ReciprocalOp> wOnly; for (auto op : ops) { op.tau = rvec3_t(0,0,0); wOnly.push_back(op); }

    const rvec3_t tau0(0.25,0.25,0.25);
    auto sf=[&](const ivec3_t& m){ double p=-2.0*PI*(m.x*tau0.x+m.y*tau0.y+m.z*tau0.z); return 1.0+std::polar(1.0,p); };
    auto S = Orbit({ivec3_t(1,0,0), ivec3_t(1,1,1)}, ops);
    ΔG_Map rg; for (const auto& m : S) rg[m] = sf(m);

    ΔG_Map out = SymmetrizeGMap(rg, wOnly);
    double maxdev=0.0; for (const auto& m : S) maxdev=std::max(maxdev, std::abs(out.at(m)-rg.at(m)));
    EXPECT_GT(maxdev, 1e-3) << "W-only (τ=0) symmetrization MUST corrupt the non-symmorphic structure factor";
}
