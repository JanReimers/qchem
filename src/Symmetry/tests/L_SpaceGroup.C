// File: Symmetry/tests/L_SpaceGroup.C  Crystal space-group detection + IBZ reduction (Tier A).
#include <gtest/gtest.h>
#include <vector>
import qchem.Symmetry.Lattice_3D.SpaceGroup;
import qchem.Symmetry.Lattice_3D.BZReduction;
import qchem.Math;   // fabs
using namespace qchem;
using namespace qchem::Symmetry::Lattice_3D;

// FCC primitive cell (columns = half-face diagonals); scale is irrelevant to symmetry.
static Matrix3D<double> FCC(double a=1.0)
{
    return Matrix3D<double>(0.0,   a/2, a/2,
                            a/2,   0.0, a/2,
                            a/2,   a/2, 0.0);
}

//---------------------------------------------------------------------------------------
//  Detection.
//
TEST(SpaceGroup, FCC_Si_diamond_is_Oh_nonsymmorphic)
{
    // Diamond Si: FCC lattice, two-atom basis at (0,0,0) and (1/4,1/4,1/4) (primitive coords).
    std::vector<AtomSite> basis = {
        {14, rvec3_t(0.0,  0.0,  0.0 )},
        {14, rvec3_t(0.25, 0.25, 0.25)},
    };
    SpaceGroup sg = SpaceGroup::Detect(FCC(), basis);

    // Crystal point group m-3m = O_h has 48 operations.
    EXPECT_EQ(sg.Order(), 48u);
    EXPECT_EQ(sg.PointGroupOps().size(), 48u);

    // Fd-3m is non-symmorphic: inversion etc. carry a (1/4,1/4,1/4) glide/screw translation.
    EXPECT_FALSE(sg.isSymmorphic());

    // Centrosymmetric -> time reversal adds nothing to the reciprocal-space ops.
    EXPECT_EQ(sg.ReciprocalPointOps(false).size(), 48u);
    EXPECT_EQ(sg.ReciprocalPointOps(true ).size(), 48u);

    // The full {U|τ} / {W|τ} density-symmetrization ops: all 48, with 24 glide ops carrying τ=(¼,¼,¼)
    // (the non-symmorphic half) and 24 symmorphic ops with τ=0.  U is the G-index scatter matrix Wᵀ.
    auto rops = sg.ReciprocalOps();
    auto dops = sg.DirectOps();
    ASSERT_EQ(rops.size(), 48u);
    ASSERT_EQ(dops.size(), 48u);
    int nGlide=0, nSymm=0;
    for (const auto& op : dops)
    {
        bool glide = fabs(op.tau.x)>1e-9 || fabs(op.tau.y)>1e-9 || fabs(op.tau.z)>1e-9;
        if (glide) { ++nGlide;
            EXPECT_NEAR(op.tau.x, 0.25, 1e-9); EXPECT_NEAR(op.tau.y, 0.25, 1e-9); EXPECT_NEAR(op.tau.z, 0.25, 1e-9); }
        else ++nSymm;
    }
    EXPECT_EQ(nGlide, 24);
    EXPECT_EQ(nSymm,  24);
    // ReciprocalOp.U is the transpose of DirectOp.W (same op ordering), the exact G-index scatter map.
    for (size_t i=0;i<rops.size();++i)
        for (int a=1;a<=3;++a) for (int b=1;b<=3;++b)
            EXPECT_NEAR(rops[i].U(a,b), dops[i].W(b,a), 1e-12);
}

TEST(SpaceGroup, SimpleCubic_monatomic_is_Oh_symmorphic)
{
    std::vector<AtomSite> basis = {{14, rvec3_t(0,0,0)}};
    SpaceGroup sg = SpaceGroup::Detect(Matrix3D<double>(), basis);  // identity cell = simple cubic a=1
    EXPECT_EQ(sg.Order(), 48u);
    EXPECT_TRUE(sg.isSymmorphic());   // all tau = 0
}

TEST(SpaceGroup, LowSymmetry_generic_cell_has_only_inversion)
{
    // A generic (triclinic) cell with a single atom: holohedry is only { E, i }.
    Matrix3D<double> A(1.0, 0.13, 0.09,
                       0.0, 1.31, 0.21,
                       0.0, 0.0,  1.07);
    std::vector<AtomSite> basis = {{14, rvec3_t(0,0,0)}};
    SpaceGroup sg = SpaceGroup::Detect(A, basis);
    EXPECT_EQ(sg.Order(), 2u);   // E and inversion
}

//---------------------------------------------------------------------------------------
//  Site stabilizers (the T2 site-group <-> angular-grid connection, SymmetryUpgradePlan §6).
//
TEST(SpaceGroup, SiteStabilizersPartitionTheGroup)
{
    // Diamond: both atom sites have the T_d site group (24 ops; orbit 2 = 48/24); a generic
    // position has only E; |orbit| * |stabilizer| == |group| throughout.
    std::vector<AtomSite> basis = {
        {14, rvec3_t(0.0,  0.0,  0.0 )},
        {14, rvec3_t(0.25, 0.25, 0.25)},
    };
    SpaceGroup sg = SpaceGroup::Detect(FCC(), basis);
    EXPECT_EQ(sg.SiteStabilizer(rvec3_t(0.0 ,0.0 ,0.0 )).size(), 24u);
    EXPECT_EQ(sg.SiteStabilizer(rvec3_t(0.25,0.25,0.25)).size(), 24u);
    EXPECT_EQ(sg.SiteStabilizer(rvec3_t(0.13,0.29,0.41)).size(), 1u);   // generic: E only

    // Monatomic simple cubic: the atom site carries the FULL O_h (orbit 1 = 48/48).
    std::vector<AtomSite> sc = {{14, rvec3_t(0,0,0)}};
    SpaceGroup cubic = SpaceGroup::Detect(Matrix3D<double>(), sc);
    EXPECT_EQ(cubic.SiteStabilizer(rvec3_t(0,0,0)).size(), 48u);
}

//---------------------------------------------------------------------------------------
//  IBZ reduction.
//
TEST(BZReduction, SimpleCubic_2x2x2_gamma_folds_to_four)
{
    std::vector<AtomSite> basis = {{14, rvec3_t(0,0,0)}};
    SpaceGroup sg = SpaceGroup::Detect(Matrix3D<double>(), basis);   // O_h

    IBZMesh ibz = ReduceToIBZ(ivec3_t(2,2,2), rvec3_t(0,0,0), sg.ReciprocalPointOps());

    // Orbits of {0,1/2}^3 under O_h: {000}, {00½}x3, {0½½}x3, {½½½} -> 4 reps.
    ASSERT_EQ(ibz.points.size(), 4u);

    std::vector<int> stars;
    for (const auto& p : ibz.points) stars.push_back(p.starSize);
    // scan order gives representatives (000),(00½),(0½½),(½½½) with sizes 1,3,3,1.
    EXPECT_EQ(stars[0], 1);
    EXPECT_EQ(stars[1], 3);
    EXPECT_EQ(stars[2], 3);
    EXPECT_EQ(stars[3], 1);

    EXPECT_NEAR(ibz.WeightSum(), 1.0, 1e-12);
    EXPECT_NEAR(ibz.points[0].weight, 1.0/8, 1e-12);
    EXPECT_NEAR(ibz.points[1].weight, 3.0/8, 1e-12);
}

TEST(BZReduction, FCC_4x4x4_gamma_partition_and_count)
{
    std::vector<AtomSite> basis = {
        {14, rvec3_t(0.0,  0.0,  0.0 )},
        {14, rvec3_t(0.25, 0.25, 0.25)},
    };
    SpaceGroup sg = SpaceGroup::Detect(FCC(), basis);

    IBZMesh ibz = ReduceToIBZ(ivec3_t(4,4,4), rvec3_t(0,0,0), sg.ReciprocalPointOps());

    // Partition invariants: every grid point owned, star sizes sum to the full grid.
    ASSERT_EQ(ibz.FullSize(), 64u);
    for (int owner : ibz.ownerOfGrid)
    {
        EXPECT_GE(owner, 0);
        EXPECT_LT(owner, int(ibz.points.size()));
    }
    int total = 0;
    for (const auto& p : ibz.points) total += p.starSize;
    EXPECT_EQ(total, 64);
    EXPECT_NEAR(ibz.WeightSum(), 1.0, 1e-12);

    // FCC 4x4x4 Gamma-centred reduces to the textbook 8 irreducible k-points.
    EXPECT_EQ(ibz.points.size(), 8u);
}

TEST(BZReduction, IdentityOps_do_not_fold)
{
    // With only E, no folding: every grid point is its own star.
    std::vector<Matrix3D<double>> onlyE = { Matrix3D<double>() };
    IBZMesh ibz = ReduceToIBZ(ivec3_t(2,2,2), rvec3_t(0,0,0), onlyE);
    EXPECT_EQ(ibz.points.size(), 8u);
    EXPECT_NEAR(ibz.WeightSum(), 1.0, 1e-12);
    for (const auto& p : ibz.points) EXPECT_EQ(p.starSize, 1);
}

//---------------------------------------------------------------------------------------
//  §7-step-7 S1: the collinear magnetic (Shubnikov) group (doc/SymmetryUpgradePlan.md).
//
//  Fixture: the MnO AFM-II rhombohedral cell -- the chemical rocksalt cell DOUBLED along [111];
//  Mn sublattices at (0,0,0) [+m] and (1/2,1/2,1/2) [-m], O (non-magnetic, spin=0) at 1/4 and 3/4.
//  The doubling is the point: Detect keeps ONE tau coset per W, so the pure ANTI-TRANSLATION
//  {E|1/2,1/2,1/2} -- which paired with Flip IS the sublattice mirror m1=-m2 -- exists only through
//  ShubnikovOps' own coset re-enumeration.
namespace
{
Matrix3D<double> MnOCell(double a=8.40)
{
    return Matrix3D<double>(a,   a/2, a/2,
                            a/2, a,   a/2,
                            a/2, a/2, a);
}
std::vector<AtomSite> MnOBasis(int s1, int s2)   // the two Mn sublattice spins; O stays spin=0
{
    return { {25, rvec3_t(0.0 ,0.0 ,0.0 ), s1},
             {25, rvec3_t(0.5 ,0.5 ,0.5 ), s2},
             { 8, rvec3_t(0.25,0.25,0.25), 0 },
             { 8, rvec3_t(0.75,0.75,0.75), 0 } };
}
bool IsOp(const SymOp& op, const Matrix3D<double>& W, const rvec3_t& tau, SpinAction s)
{
    if (op.sigma!=s) return false;
    for (int i=1;i<=3;i++) for (int j=1;j<=3;j++) if (fabs(op.W(i,j)-W(i,j))>1e-9) return false;
    for (int c=0;c<3;c++)
    {
        double d=(&op.tau.x)[c]-(&tau.x)[c]; d-=floor(d+0.5);
        if (fabs(d)>1e-6) return false;
    }
    return true;
}
} //anon

TEST(Shubnikov, MnO_AFM2_SplitsIntoSixNoneSixFlipPerCoset)
{
    std::vector<AtomSite> afm = MnOBasis(+1,-1);
    SpaceGroup sg = SpaceGroup::Detect(MnOCell(), afm);   // detection ignores the spin decoration
    const size_t nGrey = sg.Order();                       // one tau coset per W (12 for this cell)
    EXPECT_EQ(nGrey, 12u);

    auto M = sg.ShubnikovOps(afm);
    // Every detected W admits BOTH chemical-lattice cosets (the cell is a doubled chemical cell),
    // and for the AFM decoration each W contributes exactly one None and one Flip op.
    ASSERT_EQ(M.size(), 2*nGrey);
    size_t nNone=0, nFlip=0;
    for (const auto& op : M) (op.sigma==SpinAction::None ? nNone : nFlip)++;
    EXPECT_EQ(nNone, nGrey);
    EXPECT_EQ(nFlip, nGrey);

    // The two named ops the machinery exists for: the identity, and the ANTI-TRANSLATION
    // {E|1/2,1/2,1/2}*Flip = the sublattice mirror m1=-m2 (invisible to Detect's one-coset rule).
    Matrix3D<double> E; // identity
    bool haveId=false, haveAnti=false;
    for (const auto& op : M)
    {
        haveId   = haveId   || IsOp(op, E, rvec3_t(0,0,0),       SpinAction::None);
        haveAnti = haveAnti || IsOp(op, E, rvec3_t(0.5,0.5,0.5), SpinAction::Flip);
    }
    EXPECT_TRUE(haveId);
    EXPECT_TRUE(haveAnti) << "the anti-translation (the m1=-m2 mirror) must be in the Shubnikov group";
}

TEST(Shubnikov, GroupClosure_AFM2)
{
    // {W1|t1,s1}.{W2|t2,s2} = {W1 W2 | W1 t2 + t1, s1 xor s2} must land in the set (M is a group).
    std::vector<AtomSite> afm = MnOBasis(+1,-1);
    SpaceGroup sg = SpaceGroup::Detect(MnOCell(), afm);
    auto M = sg.ShubnikovOps(afm);
    for (const auto& p : M)
        for (const auto& q : M)
        {
            Matrix3D<double> W = p.W*q.W;
            rvec3_t          t = p.W*q.tau + p.tau;
            SpinAction       s = (p.sigma==q.sigma) ? SpinAction::None : SpinAction::Flip;
            bool member=false;
            for (const auto& r : M) if (IsOp(r, W, t, s)) { member=true; break; }
            ASSERT_TRUE(member) << "Shubnikov set is not closed under composition";
        }
}

TEST(Shubnikov, UndecoratedAndFMBasesGiveAllNone)
{
    // No decoration (all spin=0): the grey group, every coset sigma=None, nothing dropped.
    std::vector<AtomSite> grey = MnOBasis(0,0);
    SpaceGroup sg = SpaceGroup::Detect(MnOCell(), grey);
    auto G = sg.ShubnikovOps(grey);
    ASSERT_EQ(G.size(), 2*sg.Order());
    for (const auto& op : G) EXPECT_EQ(op.sigma, SpinAction::None);

    // FM decoration (+,+): every op still maps + onto +, so again all None -- including the plain
    // translation {E|1/2,1/2,1/2} (None here, Flip under AFM: the SAME spatial op classifies
    // differently per ordering, which is exactly what CommonOps below must respect).
    std::vector<AtomSite> fm = MnOBasis(+1,+1);
    auto F = SpaceGroup::Detect(MnOCell(), fm).ShubnikovOps(fm);
    ASSERT_EQ(F.size(), G.size());
    for (const auto& op : F) EXPECT_EQ(op.sigma, SpinAction::None);
}

TEST(Shubnikov, CommonOpsOfAFMAndFMKeepOnlySublatticePreservingOps)
{
    // The §3 ordering-search subgroup: impose only what BOTH candidate orderings share, so the SCF
    // still chooses between them.  AFM-II and FM share the 12 sublattice-PRESERVING ops (None in
    // both); the 12 sublattice-swapping ops differ in sigma (Flip vs None) and must be excluded.
    std::vector<AtomSite> afm = MnOBasis(+1,-1), fm = MnOBasis(+1,+1);
    SpaceGroup sg = SpaceGroup::Detect(MnOCell(), afm);
    auto Mafm = sg.ShubnikovOps(afm);
    auto Mfm  = sg.ShubnikovOps(fm);
    auto C    = CommonOps({Mafm, Mfm});
    ASSERT_EQ(C.size(), 12u);
    for (const auto& op : C) EXPECT_EQ(op.sigma, SpinAction::None);
    // ...and the excluded half really is the swap half: no common op moves Mn1 off its sublattice.
    for (const auto& op : C)
    {
        rvec3_t g = op.W*rvec3_t(0,0,0) + op.tau;   // image of Mn1
        double d[3]={g.x-floor(g.x+0.5), g.y-floor(g.y+0.5), g.z-floor(g.z+0.5)};
        EXPECT_LT(fabs(d[0])+fabs(d[1])+fabs(d[2]), 1e-6) << "a common op must fix the Mn1 sublattice";
    }
}

TEST(Shubnikov, AnIncompatibleDecorationDropsOps)
{
    // A FERRI-style decoration on a cell whose geometry demands equivalence: tag the two Mn with
    // UNEQUAL magnitudes' signs is not expressible here (spin is a sign), so instead break the
    // pattern with a THIRD magnetic site: simple-cubic 3-atom chain cell, spins (+,+,-) -- the
    // translation by 1/3 maps + onto + and + onto -, satisfying NEITHER rule: it must be dropped,
    // while the identity survives.  (Geometry: a=1 cubic cell, one species at z=0,1/3,2/3.)
    Matrix3D<double> A(1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 3.0);
    std::vector<AtomSite> ferri = { {25, rvec3_t(0,0,0.0      ), +1},
                                    {25, rvec3_t(0,0,1.0/3.0  ), +1},
                                    {25, rvec3_t(0,0,2.0/3.0  ), -1} };
    SpaceGroup sg = SpaceGroup::Detect(A, ferri);
    auto M = sg.ShubnikovOps(ferri);
    // The 1/3 and 2/3 translations must be gone; the identity must remain.
    Matrix3D<double> E;
    bool haveId=false;
    for (const auto& op : M)
    {
        haveId = haveId || IsOp(op, E, rvec3_t(0,0,0), SpinAction::None);
        EXPECT_FALSE(IsOp(op, E, rvec3_t(0,0,1.0/3.0), SpinAction::None));
        EXPECT_FALSE(IsOp(op, E, rvec3_t(0,0,1.0/3.0), SpinAction::Flip));
        EXPECT_FALSE(IsOp(op, E, rvec3_t(0,0,2.0/3.0), SpinAction::None));
        EXPECT_FALSE(IsOp(op, E, rvec3_t(0,0,2.0/3.0), SpinAction::Flip));
    }
    EXPECT_TRUE(haveId);
}
