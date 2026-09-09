// File: Symmetry/tests/L_Fold.C  Shared orbit-fold primitive (doc/SymmetryUpgradePlan.md §7 step 2).
#include <gtest/gtest.h>
#include <vector>
import qchem.Symmetry.Lattice_3D.Fold;
import qchem.Symmetry.Lattice_3D.SpaceGroup;
import qchem.Symmetry.Lattice_3D.BZReduction;
import qchem.Math;   // fabs, lround
using namespace qchem;
using namespace qchem::Symmetry::Lattice_3D;

// FCC primitive cell (columns = half-face diagonals); scale is irrelevant to symmetry.
static Matrix3D<double> FCC(double a=1.0)
{
    return Matrix3D<double>(0.0,   a/2, a/2,
                            a/2,   0.0, a/2,
                            a/2,   a/2, 0.0);
}

static SpaceGroup DiamondSi()
{
    std::vector<AtomSite> basis = {
        {14, rvec3_t(0.0,  0.0,  0.0 )},
        {14, rvec3_t(0.25, 0.25, 0.25)},
    };
    return SpaceGroup::Detect(FCC(), basis);
}

static int LinearIndex(const ivec3_t& i, const ivec3_t& N)
{
    return (i.x * N.y + i.y) * N.z + i.z;
}

//! Test-side replica of the grid action, for edge-op verification.  Note the CONJUGATED matrix
//! M = D U D^{-1} (D=diag(N)): a grid point is k=(i+s)/N componentwise, so i' = M(i+s)-s mod N.
//! On an isotropic mesh M==U; on an anisotropic one a non-integral M means U is not a symmetry of
//! the MESH at all and the op does not act (see FoldGrid's contract).
static bool GridImage(const Matrix3D<double>& U, const ivec3_t& i,
                      const ivec3_t& N, const rvec3_t& shift, ivec3_t& out)
{
    const int    n[3] = {N.x, N.y, N.z};
    const double s[3] = {shift.x, shift.y, shift.z};
    const double ii[3]= {double(i.x), double(i.y), double(i.z)};
    int r[3];
    for (int a = 0; a < 3; ++a)
    {
        double t = -s[a];
        for (int c = 0; c < 3; ++c)
        {
            const double m = double(n[a])*U(a+1,c+1)/double(n[c]);
            if (fabs(m - double(lround(m))) > 1e-6) return false;   // not a mesh symmetry
            t += m*(ii[c] + s[c]);
        }
        long v = lround(t);
        if (fabs(t - double(v)) > 1e-6) return false;
        r[a] = int(((v % n[a]) + n[a]) % n[a]);
    }
    out = ivec3_t(r[0], r[1], r[2]);
    return true;
}

//! The partition invariants every fold must satisfy.
static void CheckPartition(const Fold& f, size_t n)
{
    ASSERT_EQ(f.owner.size(), n);
    ASSERT_EQ(f.members.size(), f.Reps());
    ASSERT_EQ(f.starSize.size(), f.Reps());
    int total = 0;
    for (size_t r = 0; r < f.Reps(); ++r)
    {
        EXPECT_EQ(f.owner[f.repRaw[r]], int(r));                    // rep owns itself
        EXPECT_EQ(int(f.members[r].size()), f.starSize[r]);         // one edge per member
        total += f.starSize[r];
        for (auto [m, o] : f.members[r]) EXPECT_EQ(f.owner[m], int(r));
    }
    EXPECT_EQ(total, int(n));                                       // stars partition the set
    for (size_t i = 0; i < n; ++i) { EXPECT_GE(f.owner[i], 0); EXPECT_LT(f.owner[i], int(f.Reps())); }
}

//---------------------------------------------------------------------------------------
//  FoldGrid: the k-star, bit-identical to ReduceToIBZ.
//
TEST(Fold, KMeshFoldBitIdenticalToReduceToIBZ)
{
    SpaceGroup sg = DiamondSi();
    const ivec3_t N(4,4,4);
    const rvec3_t s(0,0,0);

    IBZMesh ibz = ReduceToIBZ(N, s, sg.ReciprocalPointOps());
    Fold    f   = sg.FoldKMesh(N, s);
    CheckPartition(f, ibz.FullSize());

    ASSERT_EQ(f.Reps(), ibz.points.size());
    for (size_t r = 0; r < f.Reps(); ++r)
    {
        EXPECT_EQ(f.repRaw[r],   LinearIndex(ibz.points[r].index, N));
        EXPECT_EQ(f.starSize[r], ibz.points[r].starSize);
    }
    for (size_t i = 0; i < ibz.FullSize(); ++i)
        EXPECT_EQ(f.owner[i], ibz.ownerOfGrid[i]);
}

TEST(Fold, GridEdgeOpsMapRepToMember)
{
    SpaceGroup sg = DiamondSi();
    const ivec3_t N(4,4,4);
    const rvec3_t s(0,0,0);
    auto U = sg.ReciprocalPointOps(true);
    std::vector<SymOp> ops;
    for (const auto& u : U) ops.push_back({u, rvec3_t(0,0,0)});

    Fold f = FoldGrid(N, s, ops);
    for (size_t r = 0; r < f.Reps(); ++r)
    {
        int lin = f.repRaw[r];
        ivec3_t rep(lin / (N.y*N.z), (lin / N.z) % N.y, lin % N.z);
        for (auto [m, o] : f.members[r])
        {
            ASSERT_GE(o, 0);                       // a full group: every edge has its op
            ivec3_t img;
            ASSERT_TRUE(GridImage(ops[o].W, rep, N, s, img));
            EXPECT_EQ(LinearIndex(img, N), m);     // ops[o] maps rep -> member, exactly
        }
    }
}

//---------------------------------------------------------------------------------------
//  ⛔ KP-0 (2026-09-09): THE ANISOTROPIC MESH.  A cubic op that PERMUTES AXES is not a symmetry of a
//  2x1x1 grid, but the old action applied U to the index vector and let the mod-N wrap land the image
//  on a grid point anyway -- a non-bijection, so the BFS orbits OVERLAPPED and the star sizes came out
//  1 and 2 on a TWO-point grid (Sum(w)=1.5).  Every BZ-summed quantity carries those weights, so the
//  Si multi-k SCF ran on 12 electrons instead of 8.  The action is M=DUD^{-1}, and an op whose M is
//  non-integral simply is not a mesh symmetry.
//
//  These meshes are the ones the ISOTROPIC 4x4x4 gates above cannot see: there M==U identically.
//
TEST(Fold, AnisotropicKMeshStarsStillPartition)
{
    SpaceGroup sg = DiamondSi();
    const rvec3_t s(0,0,0);
    for (const ivec3_t N : {ivec3_t(2,1,1), ivec3_t(1,2,1), ivec3_t(1,1,2),
                            ivec3_t(2,2,1), ivec3_t(4,2,1), ivec3_t(3,2,1)})
    {
        Fold f = sg.FoldKMesh(N, s);
        SCOPED_TRACE(testing::Message() << "N=" << N.x << "x" << N.y << "x" << N.z);
        CheckPartition(f, size_t(N.x)*N.y*N.z);          // stars PARTITION -- Sum starSize == Ntot
        ReduceToIBZ(N, s, sg.ReciprocalPointOps());      // and the weight-sum invariant holds (it throws)
    }
}

// THE ORIGINAL DEFECT, PINNED EXACTLY.  Si 2x1x1 has Gamma and the zone boundary k=(1/2,0,0); both are
// TRIM (k = -k modulo a reciprocal vector), so each is its OWN star and the two weights are 1/2 each.
// The reported symptom was starSize {1,2} and Sum(w)=1.5 = 12/8 electrons.
TEST(Fold, Si2x1x1WeightsSumToOne)
{
    SpaceGroup sg = DiamondSi();
    IBZMesh ibz = ReduceToIBZ(ivec3_t(2,1,1), rvec3_t(0,0,0), sg.ReciprocalPointOps());

    ASSERT_EQ(ibz.points.size(), 2u);                    // no reduction is possible: 2 TRIM points
    EXPECT_EQ(ibz.points[0].starSize, 1);
    EXPECT_EQ(ibz.points[1].starSize, 1);
    EXPECT_NEAR(ibz.WeightSum(), 1.0, 1e-14);
}

// THE SHIFTED mesh adds the second per-op condition, (M-I)s integral: an op may map the GRID onto
// itself and still miss the OFFSET point set.  A half-shifted 2x2x2 mesh is the standard MP case, and
// its weights must partition just the same.
TEST(Fold, ShiftedKMeshStarsStillPartition)
{
    SpaceGroup sg = DiamondSi();
    for (const rvec3_t s : {rvec3_t(0.5,0.5,0.5), rvec3_t(0.5,0,0), rvec3_t(0,0.5,0.5)})
    for (const ivec3_t N : {ivec3_t(2,2,2), ivec3_t(2,1,1), ivec3_t(4,2,2)})
    {
        SCOPED_TRACE(testing::Message() << "N=" << N.x << "x" << N.y << "x" << N.z
                                        << " shift=" << s.x << "," << s.y << "," << s.z);
        CheckPartition(sg.FoldKMesh(N, s), size_t(N.x)*N.y*N.z);
        ReduceToIBZ(N, s, sg.ReciprocalPointOps());
    }
}

//---------------------------------------------------------------------------------------
//  FoldGVectors: the T1 {G}-sum bookkeeping -- representatives x n_star reproduce the
//  full sum exactly (integer arithmetic, so equality is exact, not 1e-13).
//
TEST(Fold, GBallFoldsToShellStarsAndSumsExactly)
{
    std::vector<AtomSite> basis = {{14, rvec3_t(0,0,0)}};
    SpaceGroup sg = SpaceGroup::Detect(Matrix3D<double>(), basis);   // simple cubic, O_h

    std::vector<ivec3_t> ball;                                       // |m|^2 <= 4
    for (int x = -2; x <= 2; ++x)
    for (int y = -2; y <= 2; ++y)
    for (int z = -2; z <= 2; ++z)
        if (x*x + y*y + z*z <= 4) ball.push_back(ivec3_t(x,y,z));

    Fold f = sg.FoldGVectors(ball);
    CheckPartition(f, ball.size());

    // O_h orbits of the ball: (000), (100)x6, (110)x12, (111)x8, (200)x6 -> 5 stars.
    EXPECT_EQ(f.Reps(), 5u);

    // A totally symmetric G-sum folds exactly: sum |m|^2 over the ball == sum over reps of n_star*|m_rep|^2.
    long full = 0, reduced = 0;
    for (const auto& m : ball) full += m.x*m.x + m.y*m.y + m.z*m.z;
    for (size_t r = 0; r < f.Reps(); ++r)
    {
        const ivec3_t& m = ball[f.repRaw[r]];
        reduced += long(f.starSize[r]) * (m.x*m.x + m.y*m.y + m.z*m.z);
    }
    EXPECT_EQ(reduced, full);

    // Edge ops reproduce every member exactly: member = U_op * rep.
    for (size_t r = 0; r < f.Reps(); ++r)
    {
        const ivec3_t& rep = ball[f.repRaw[r]];
        auto rops = sg.ReciprocalOps();
        for (auto [mi, o] : f.members[r])
        {
            ASSERT_GE(o, 0);
            rvec3_t im = rops[o].U * rvec3_t(rep.x, rep.y, rep.z);
            EXPECT_EQ(lround(im.x), ball[mi].x);
            EXPECT_EQ(lround(im.y), ball[mi].y);
            EXPECT_EQ(lround(im.z), ball[mi].z);
        }
    }
}

//---------------------------------------------------------------------------------------
//  FoldPoints: the tolerance path, and the non-symmorphic review pin -- tau ACTS on {r},
//  so a bare-W fold cannot merge glide-related points.
//
TEST(Fold, TauActsOnPointsBareWCannotFold)
{
    std::vector<rvec3_t> pts = {rvec3_t(0.10, 0.20, 0.30),
                                rvec3_t(0.35, 0.45, 0.55)};          // = pts[0] + (1/4,1/4,1/4)
    Matrix3D<double> E;                                              // identity

    std::vector<SymOp> withTau = {{E, rvec3_t(0,0,0)}, {E, rvec3_t(0.25,0.25,0.25)}};
    Fold f = FoldPoints(pts, withTau, 1e-9);
    CheckPartition(f, pts.size());
    ASSERT_EQ(f.Reps(), 1u);                                         // {E|tau} merges the pair
    EXPECT_EQ(f.starSize[0], 2);

    std::vector<SymOp> bareW = {{E, rvec3_t(0,0,0)}};                // tau stripped: no fold
    Fold g = FoldPoints(pts, bareW, 1e-9);
    EXPECT_EQ(g.Reps(), 2u);
}

TEST(Fold, DiamondBasisSitesFoldToOneOrbit)
{
    // The two diamond Si sites are symmetry-equivalent (swapped by the inversion glide
    // {-I|(1/4,1/4,1/4)}): the atom-permutation reality of plan §2b item 4.
    SpaceGroup sg = DiamondSi();
    std::vector<rvec3_t> sites = {rvec3_t(0,0,0), rvec3_t(0.25,0.25,0.25)};
    Fold f = sg.FoldPoints(sites);
    CheckPartition(f, sites.size());
    ASSERT_EQ(f.Reps(), 1u);
    EXPECT_EQ(f.starSize[0], 2);
}

TEST(Fold, IdentityOnlyGivesSingletonStars)
{
    std::vector<rvec3_t> pts = {rvec3_t(0.1,0,0), rvec3_t(0,0.2,0), rvec3_t(0,0,0.3)};
    std::vector<SymOp> onlyE = {{Matrix3D<double>(), rvec3_t(0,0,0)}};
    Fold f = FoldPoints(pts, onlyE, 1e-9);
    CheckPartition(f, pts.size());
    EXPECT_EQ(f.Reps(), 3u);
    for (size_t r = 0; r < f.Reps(); ++r) EXPECT_EQ(f.starSize[r], 1);
}

//---------------------------------------------------------------------------------------
//  Tier 4a (plan §4): the op type is spin-native from the start -- a pure spatial op has
//  sigma = None; the fold geometry ignores sigma (it acts on field channels, not points).
//
TEST(Fold, SymOpSpinActionDefaultsToNone)
{
    SymOp op{Matrix3D<double>(), rvec3_t(0,0,0)};
    EXPECT_EQ(op.sigma, SpinAction::None);
    SymOp tr{Matrix3D<double>(), rvec3_t(0,0,0), SpinAction::Flip};
    EXPECT_EQ(tr.sigma, SpinAction::Flip);
}
