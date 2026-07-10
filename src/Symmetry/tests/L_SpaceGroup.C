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
