// File: UnitTests/M_PointGroup.C  Point-group symmetry primitives (stage 1a).
#include <gtest/gtest.h>
#include <vector>
import qchem.Symmetry.PointGroup;
import qchem.Math;            // Pi, sin, cos, fabs (project-wide, for test geometry)

using namespace Symmetry;

// Count axes of a given order.
static int CountOrder(const std::vector<RotationAxis>& a, int order)
{
    int n=0; for (const auto& ax : a) if (ax.order==order) ++n; return n;
}

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
    return fabs(pa.axis[pa.uniqueAxis] * z);
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
    const double c=cos(2.0*Pi/3.0), s=sin(2.0*Pi/3.0), rho=0.94;
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

// --- Proper rotation-axis finder (stage 1b-2a) ---------------------------------------

static std::vector<SymPoint> Benzene()  // D6h, ring in the xy-plane, C6 axis = z
{
    std::vector<SymPoint> b;
    const double Rc=1.39, Rh=2.47;
    for (int k=0;k<6;k++)
    {
        double a = k*Pi/3.0;                          // 60 degrees apart
        b.push_back({6, rvec3_t(Rc*cos(a), Rc*sin(a), 0.0)});
        b.push_back({1, rvec3_t(Rh*cos(a), Rh*sin(a), 0.0)});
    }
    return b;
}

TEST(PointGroup, Axes_water)
{
    auto w = Water();
    auto ax = FindRotationAxes(w, Centroid(w), 1e-6);
    EXPECT_EQ(ax.size(), 1u);
    EXPECT_EQ(ax[0].order, 2);
}

TEST(PointGroup, Axes_ammonia_single_C3)
{
    const double c=cos(2.0*Pi/3.0), s=sin(2.0*Pi/3.0), rho=0.94;
    std::vector<SymPoint> nh3 = {
        {7, rvec3_t(0.0, 0.0, 0.10)},
        {1, rvec3_t(rho,   0.0,   -0.30)},
        {1, rvec3_t(rho*c, rho*s, -0.30)},
        {1, rvec3_t(rho*c,-rho*s, -0.30)},
    };
    auto ax = FindRotationAxes(nh3, Centroid(nh3), 1e-6);
    EXPECT_EQ(ax.size(), 1u);     // C3v has only the one C3 axis (no perpendicular C2)
    EXPECT_EQ(ax[0].order, 3);
}

TEST(PointGroup, Axes_benzene_D6h)
{
    auto b = Benzene();
    auto ax = FindRotationAxes(b, Centroid(b), 1e-6);
    EXPECT_EQ(ax[0].order, 6);          // principal axis C6
    EXPECT_EQ(CountOrder(ax, 6), 1);
    EXPECT_EQ(CountOrder(ax, 2), 6);    // six C2 axes perpendicular to C6
    EXPECT_EQ(ax.size(), 7u);
}

TEST(PointGroup, Axes_methane_Td)
{
    std::vector<SymPoint> ch4 = {
        {6, rvec3_t( 0, 0, 0)},
        {1, rvec3_t( 1, 1, 1)},
        {1, rvec3_t( 1,-1,-1)},
        {1, rvec3_t(-1, 1,-1)},
        {1, rvec3_t(-1,-1, 1)},
    };
    auto ax = FindRotationAxes(ch4, Centroid(ch4), 1e-6);
    EXPECT_EQ(CountOrder(ax, 3), 4);    // four C3 axes through the C-H bonds
    EXPECT_EQ(CountOrder(ax, 2), 3);    // three C2 axes through opposite H-H edge midpoints
    EXPECT_EQ(ax.size(), 7u);
}

// --- Mirror planes / inversion / improper axes inventory (stage 1b-2b) ---------------

static std::vector<SymPoint> Ammonia()
{
    const double c=cos(2.0*Pi/3.0), s=sin(2.0*Pi/3.0), rho=0.94;
    return {
        {7, rvec3_t(0.0, 0.0, 0.10)},
        {1, rvec3_t(rho,   0.0,   -0.30)},
        {1, rvec3_t(rho*c, rho*s, -0.30)},
        {1, rvec3_t(rho*c,-rho*s, -0.30)},
    };
}
static std::vector<SymPoint> Methane()
{
    return {
        {6, rvec3_t( 0, 0, 0)},
        {1, rvec3_t( 1, 1, 1)}, {1, rvec3_t( 1,-1,-1)},
        {1, rvec3_t(-1, 1,-1)}, {1, rvec3_t(-1,-1, 1)},
    };
}

TEST(PointGroup, Inventory_water_C2v)
{
    auto w = Water(); rvec3_t o = Centroid(w);
    EXPECT_EQ(FindMirrorPlanes(w,o,1e-6).size(), 2u);   // 2 sigma_v
    EXPECT_FALSE(HasInversion(w,o,1e-6));
    EXPECT_EQ(FindImproperAxes(w,o,1e-6).size(), 0u);
}

TEST(PointGroup, Inventory_ammonia_C3v)
{
    auto a = Ammonia(); rvec3_t o = Centroid(a);
    EXPECT_EQ(FindMirrorPlanes(a,o,1e-6).size(), 3u);   // 3 sigma_v
    EXPECT_FALSE(HasInversion(a,o,1e-6));
    EXPECT_EQ(FindImproperAxes(a,o,1e-6).size(), 0u);
}

TEST(PointGroup, Inventory_benzene_D6h)
{
    auto b = Benzene(); rvec3_t o = Centroid(b);
    EXPECT_EQ(FindMirrorPlanes(b,o,1e-6).size(), 7u);   // sigma_h + 3 sigma_v + 3 sigma_d
    EXPECT_TRUE(HasInversion(b,o,1e-6));
    auto s = FindImproperAxes(b,o,1e-6);
    EXPECT_EQ(s.size(), 1u);                            // S6 along the C6 axis
    EXPECT_EQ(s[0].order, 6);
}

TEST(PointGroup, Inventory_methane_Td)
{
    auto m = Methane(); rvec3_t o = Centroid(m);
    EXPECT_EQ(FindMirrorPlanes(m,o,1e-6).size(), 6u);   // 6 sigma_d
    EXPECT_FALSE(HasInversion(m,o,1e-6));
    auto s = FindImproperAxes(m,o,1e-6);
    EXPECT_EQ(s.size(), 3u);                            // 3 S4 axes
    EXPECT_EQ(s[0].order, 4);
}

// --- Full detection + Schoenflies symbol + abelian descent (stage 1b-2c) -------------

TEST(PointGroup, Detect_water_C2v)
{
    auto g = DetectPointGroup(Water(), 1e-6);
    EXPECT_EQ(g.symbol, "C2v");
    EXPECT_EQ(g.abelian, "C2v");
    EXPECT_EQ(g.order, 4);
    EXPECT_EQ(g.principalOrder, 2);
}

TEST(PointGroup, Detect_ammonia_C3v)
{
    auto g = DetectPointGroup(Ammonia(), 1e-6);
    EXPECT_EQ(g.symbol, "C3v");
    EXPECT_EQ(g.abelian, "Cs");       // C3 has no real abelian subgroup beyond a mirror
    EXPECT_EQ(g.order, 6);
}

TEST(PointGroup, Detect_benzene_D6h)
{
    auto g = DetectPointGroup(Benzene(), 1e-6);
    EXPECT_EQ(g.symbol, "D6h");
    EXPECT_EQ(g.abelian, "D2h");
    EXPECT_EQ(g.order, 24);
    EXPECT_TRUE(g.hasInversion);
    EXPECT_EQ(g.nC2perp, 6);
}

TEST(PointGroup, Detect_methane_Td)
{
    auto g = DetectPointGroup(Methane(), 1e-6);
    EXPECT_EQ(g.symbol, "Td");
    EXPECT_EQ(g.abelian, "C2v");
    EXPECT_EQ(g.order, 24);
    EXPECT_FALSE(g.hasInversion);
}

TEST(PointGroup, Detect_co2_Dinfh)
{
    std::vector<SymPoint> co2 = {
        {6, rvec3_t(0,0, 0.00)}, {8, rvec3_t(0,0, 1.16)}, {8, rvec3_t(0,0,-1.16)},
    };
    auto g = DetectPointGroup(co2, 1e-6);
    EXPECT_EQ(g.symbol, "D∞h");
    EXPECT_EQ(g.abelian, "D2h");
    EXPECT_EQ(g.top, TopType::Linear);
}

TEST(PointGroup, Detect_rectangle_D2h)
{
    // Four identical atoms at the corners of a (non-square) rectangle in the xy-plane.
    std::vector<SymPoint> rect = {
        {5, rvec3_t( 1.5,  1.0, 0)}, {5, rvec3_t(-1.5,  1.0, 0)},
        {5, rvec3_t(-1.5, -1.0, 0)}, {5, rvec3_t( 1.5, -1.0, 0)},
    };
    auto g = DetectPointGroup(rect, 1e-6);
    EXPECT_EQ(g.symbol, "D2h");
    EXPECT_EQ(g.abelian, "D2h");
    EXPECT_EQ(g.order, 8);
    EXPECT_TRUE(g.hasInversion);
}

TEST(PointGroup, Detect_scalene_planar_Cs)
{
    // Three distinct atoms all in the y=0 plane and no other symmetry -> Cs.
    std::vector<SymPoint> cs = {
        {8, rvec3_t(0,0,0)}, {1, rvec3_t(1,0,0)}, {6, rvec3_t(0,0,1)},
    };
    auto g = DetectPointGroup(cs, 1e-6);
    EXPECT_EQ(g.symbol, "Cs");
    EXPECT_EQ(g.abelian, "Cs");
    EXPECT_EQ(g.order, 2);
    EXPECT_EQ(g.nSigma, 1);
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
