// File: UnitTests/M_PointGroup.C  Point-group symmetry primitives (stage 1a).
#include <gtest/gtest.h>
#include <vector>
#include <cmath>
import qchem.Symmetry.PointGroup;

using namespace Symmetry;

// Water, C2v, lying in the yz-plane with its C2 axis along z.  (Geometry in arbitrary units;
// only the symmetry matters.)  The centroid is NOT at an atom, which exercises the
// origin-centering in IsSymmetryOf.
static std::vector<SymPoint> Water()
{
    return {
        {8, rvec3_t(0.0,  0.000,  0.117)},   // O on the C2 (z) axis
        {1, rvec3_t(0.0,  0.757, -0.467)},   // H
        {1, rvec3_t(0.0, -0.757, -0.467)},   // H
    };
}

TEST(PointGroup, Water_is_C2v)
{
    auto w = Water();
    rvec3_t o = Centroid(w);
    const double tol = 1e-6;
    rvec3_t x(1,0,0), y(0,1,0), z(0,0,1);

    // The four C2v operations { E, C2(z), sigma_v(xz), sigma_v'(yz) } are symmetries.
    EXPECT_TRUE(IsSymmetryOf(SymOp::E(),      w, o, tol));
    EXPECT_TRUE(IsSymmetryOf(SymOp::Cn(z,2),  w, o, tol)); // C2 about z
    EXPECT_TRUE(IsSymmetryOf(SymOp::Sigma(y), w, o, tol)); // sigma(xz), normal along y
    EXPECT_TRUE(IsSymmetryOf(SymOp::Sigma(x), w, o, tol)); // sigma(yz), normal along x

    // Operations water does NOT have.
    EXPECT_FALSE(IsSymmetryOf(SymOp::Inversion(), w, o, tol)); // no inversion centre
    EXPECT_FALSE(IsSymmetryOf(SymOp::Cn(z,3),     w, o, tol)); // no C3
    EXPECT_FALSE(IsSymmetryOf(SymOp::Cn(x,2),     w, o, tol)); // no C2 perpendicular to z
    EXPECT_FALSE(IsSymmetryOf(SymOp::Sigma(z),    w, o, tol)); // no sigma_h (xy plane)
}

// --- Molecular-top classification (stage 1b-1) ---------------------------------------

static double AxisAlongZ(const PrincipalAxes& pa)            // |unique axis . z|, ~1 if along z
{
    rvec3_t z(0,0,1);
    return std::fabs(pa.axis[pa.uniqueAxis] * z);
}

TEST(PointGroup, Top_water_asymmetric)
{
    auto w = Water();
    PrincipalAxes pa = ClassifyTop(w, Centroid(w));
    EXPECT_EQ(pa.top, TopType::Asymmetric);
}

TEST(PointGroup, Top_ammonia_symmetric_axis_z)
{
    // NH3, C3 axis along z: N near the apex, 3 H in a lower plane 120 deg apart.
    const double c=std::cos(2.0*M_PI/3.0), s=std::sin(2.0*M_PI/3.0), rho=0.94;
    std::vector<SymPoint> nh3 = {
        {7, rvec3_t(0.0, 0.0, 0.10)},
        {1, rvec3_t(rho,        0.0,        -0.30)},
        {1, rvec3_t(rho*c,      rho*s,      -0.30)},
        {1, rvec3_t(rho*c,     -rho*s,      -0.30)},
    };
    PrincipalAxes pa = ClassifyTop(nh3, Centroid(nh3));
    EXPECT_EQ(pa.top, TopType::Symmetric);
    EXPECT_NEAR(AxisAlongZ(pa), 1.0, 1e-3);   // unique (C3) axis lies along z
}

TEST(PointGroup, Top_methane_spherical)
{
    // CH4, tetrahedral H at alternate cube corners.
    std::vector<SymPoint> ch4 = {
        {6, rvec3_t( 0, 0, 0)},
        {1, rvec3_t( 1, 1, 1)},
        {1, rvec3_t( 1,-1,-1)},
        {1, rvec3_t(-1, 1,-1)},
        {1, rvec3_t(-1,-1, 1)},
    };
    PrincipalAxes pa = ClassifyTop(ch4, Centroid(ch4));
    EXPECT_EQ(pa.top, TopType::Spherical);
}

TEST(PointGroup, Top_co2_linear_axis_z)
{
    // CO2 on the z-axis.
    std::vector<SymPoint> co2 = {
        {6, rvec3_t(0,0, 0.00)},
        {8, rvec3_t(0,0, 1.16)},
        {8, rvec3_t(0,0,-1.16)},
    };
    PrincipalAxes pa = ClassifyTop(co2, Centroid(co2));
    EXPECT_EQ(pa.top, TopType::Linear);
    EXPECT_NEAR(AxisAlongZ(pa), 1.0, 1e-6);   // C_inf axis = z
}

TEST(PointGroup, SymOp_basics)
{
    rvec3_t z(0,0,1);
    EXPECT_TRUE (SymOp::Cn(z,2).IsProper());
    EXPECT_FALSE(SymOp::Sigma(z).IsProper());
    EXPECT_FALSE(SymOp::Inversion().IsProper());
    EXPECT_EQ   (SymOp::Cn(z,3).GetOrder(), 3);
    EXPECT_EQ   (SymOp::Cn(z,2).Label(), "C2");
    EXPECT_EQ   (SymOp::Inversion().Label(), "i");

    // C2 about z maps (x,y,z) -> (-x,-y,z).
    rvec3_t r(1.0, 2.0, 3.0);
    rvec3_t Rr = SymOp::Cn(z,2).Apply(r);
    EXPECT_NEAR(Rr.x, -1.0, 1e-12);
    EXPECT_NEAR(Rr.y, -2.0, 1e-12);
    EXPECT_NEAR(Rr.z,  3.0, 1e-12);
}
