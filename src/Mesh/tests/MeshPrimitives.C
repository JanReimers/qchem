// File: src/Mesh/tests/MeshPrimitives.C  Smoke + acceptance tests for the qcMesh primitives.
#include "gtest/gtest.h"
#include <cmath>
import qchem.Mesh.Product;       // Mesh, ProductMesh, MakeRadial/MakeAngular, MeshParams, enums
import qchem.Mesh.Quadrature;    // free functions (over qcMath Scalar/VectorFunction)
import qchem.Mesh.GaussLegendre;
import qchem.Blaze;
import qchem.Math;                // Pi, Pi32, FourPi
using namespace qchem;

using namespace qchem::qcMesh;          // Mesh, MakeRadial/MakeAngular, quadrature, ...

//================================================================================================
//  Test fields: a single normalised-ish s-function phi(r) = exp(-2|r|), and some scalar potentials.
//================================================================================================
class ExpBasis : public VectorFunction<double>
{
public:
    size_t GetVectorSize() const override {return 1;}
    rvec_t operator()(const rvec3_t& r) const override {return {std::exp(-2*norm(r))};}
    rvec3vec_t Gradient(const rvec3_t& r) const override
    {
        double mr=norm(r);
        if (mr==0.0) return {rvec3_t(0,0,0)};
        return {(r/mr)*(-2*std::exp(-2*mr))};
    }
};

class ExpScalar : public ScalarFunction<double>     // exp(-2r)
{
public:
    double   operator()(const rvec3_t& r) const override {return std::exp(-2*norm(r));}
    rvec3_t  Gradient  (const rvec3_t&  ) const override {return rvec3_t(0,0,0);}
};

class GaussScalar : public ScalarFunction<double>   // exp(-r^2)
{
public:
    double   operator()(const rvec3_t& r) const override {double m=norm(r); return std::exp(-m*m);}
    rvec3_t  Gradient  (const rvec3_t&  ) const override {return rvec3_t(0,0,0);}
};

class OneOverR : public ScalarFunction<double>      // 1/r, finite (0) at the origin
{
public:
    double   operator()(const rvec3_t& r) const override {double m=norm(r); return m==0.0 ? 0.0 : 1.0/m;}
    rvec3_t  Gradient  (const rvec3_t&  ) const override {return rvec3_t(0,0,0);}
};

static double Sum(const rvec_t& v) {double s=0; for (size_t i=0;i<v.size();i++) s+=v[i]; return s;}

//================================================================================================
//  1. Shared 1D Gauss-Legendre: known table + polynomial exactness to degree 2n-1.
//================================================================================================
TEST(Mesh_GaussLegendre, TwoPointTable)
{
    GaussLegendre gl(2,-1.0,1.0);
    double r=1.0/std::sqrt(3.0);
    EXPECT_NEAR(gl.x[0],-r,1e-14);
    EXPECT_NEAR(gl.x[1], r,1e-14);
    EXPECT_NEAR(gl.w[0],1.0,1e-14);
    EXPECT_NEAR(gl.w[1],1.0,1e-14);
}

TEST(Mesh_GaussLegendre, PolynomialExactness)
{
    const int n=5;                              // exact up to degree 2n-1 = 9
    GaussLegendre gl(n,-1.0,1.0);
    for (int k=0; k<=2*n-1; k++)
    {
        double q=0;
        for (int i=0;i<n;i++) q+=gl.w[i]*std::pow(gl.x[i],k);
        double exact=(k%2==1) ? 0.0 : 2.0/(k+1);   // integral_{-1}^{1} x^k dx
        EXPECT_NEAR(q,exact,1e-12) << "degree " << k;
    }
}

//================================================================================================
//  2. Radial meshes: integral_0^inf r^2 exp(-2r) dr = 0.25  (the r^2 jacobian is folded into W).
//================================================================================================
static double RadialIntegralExp(const MeshParams& p)
{
    RadialMesh m=MakeRadial(p);
    double I=0;
    for (size_t i=0;i<m.size();i++) I+=m.W()[i]*std::exp(-2*m.R()[i]);
    return I;
}

TEST(Mesh_Radial, MHL)
{
    MeshParams p; p.radial=RadialKind::MHL; p.nRadial=200; p.mhl_m=2; p.mhl_alpha=3.0;
    EXPECT_NEAR(RadialIntegralExp(p),0.25,4e-15);
}
TEST(Mesh_Radial, Log)
{
    MeshParams p; p.radial=RadialKind::Log; p.nRadial=400; p.logStart=1e-4; p.logStop=40;
    EXPECT_NEAR(RadialIntegralExp(p),0.25,1e-3);
}
TEST(Mesh_Radial, Linear)
{
    MeshParams p; p.radial=RadialKind::Linear; p.nRadial=4000; p.logStop=40;
    EXPECT_NEAR(RadialIntegralExp(p),0.25,1e-3);
}

//================================================================================================
//  3. Angular meshes: sum of weights = 4*pi, and integral z^2 dOmega = 4*pi/3.
//================================================================================================
// sum W = 4*pi, and integral z^2 dOmega = 4*pi/3 to the scheme's tolerance ztol.
static void CheckAngular(const MeshParams& p, double sumtol, double ztol)
{
    AngularMesh a=MakeAngular(p);
    EXPECT_NEAR(Sum(a.W()),FourPi,sumtol);
    double zz=0;
    for (size_t j=0;j<a.size();j++) zz+=a.W()[j]*a.Dirs()[j].z*a.Dirs()[j].z;
    EXPECT_NEAR(zz,FourPi/3.0,ztol);
}
TEST(Mesh_Angular, Gauss)
{
    // numDir=32 was removed (inherited weight-sum bug, see LebedevAngularMesh.C).
    // 24 and 30 have ~7-figure direction constants (not unit vectors) -> exact only to ~1e-7.
    CheckAngular({.angular=AngularKind::Lebedev, .nAngular=12}, 1e-10, 1e-10);
    CheckAngular({.angular=AngularKind::Lebedev, .nAngular=24}, 1e-10, 1e-7);
    CheckAngular({.angular=AngularKind::Lebedev, .nAngular=30}, 1e-10, 1e-7);
    CheckAngular({.angular=AngularKind::Lebedev, .nAngular=50}, 1e-10, 1e-10);
}
TEST(Mesh_Angular, GaussLegendre)   // Gauss-exact in cos(theta) -> z^2 exact
{
    CheckAngular({.angular=AngularKind::GaussLegendre, .nAngular=11}, 1e-10, 1e-12);
}
TEST(Mesh_Angular, EulerMaclaren)   // both the sum and z^2 are approximations for this scheme
{
    CheckAngular({.angular=AngularKind::EulerMaclaren, .nAngular=23, .em_m=2}, 2e-3, 5e-3);
}

//================================================================================================
//  4. ProductMesh + Integrate(ScalarField).
//================================================================================================
static Mesh MakeProduct()
{
    MeshParams p;
    p.radial=RadialKind::MHL; p.nRadial=200; p.mhl_m=2; p.mhl_alpha=3.0;
    p.angular=AngularKind::Lebedev; p.nAngular=12;
    RadialMesh  rad=MakeRadial(p);
    AngularMesh ang=MakeAngular(p);
    return ProductMesh(rad,ang);
}

TEST(Mesh_Product, ExpIntegral)        // integral exp(-2r) d^3r = 4*pi * 0.25 = pi
{
    Mesh m=MakeProduct();
    EXPECT_NEAR(Integrate(m,ExpScalar()),Pi,1e-12);
}
TEST(Mesh_Product, GaussianIntegral)   // integral exp(-r^2) d^3r = pi^{3/2}
{
    Mesh m=MakeProduct();
    EXPECT_NEAR(Integrate(m,GaussScalar()),Pi32,1e-6);
}

//================================================================================================
//  5. Quadrature free functions against analytic <exp(-2r)|...|exp(-2r)>.
//================================================================================================
TEST(Mesh_Quadrature, OverlapAndNormalize)
{
    Mesh m=MakeProduct();
    ExpBasis a;
    auto S=Overlap(m,a);                       // <a|a> = integral exp(-4r) d^3r = pi/8
    ASSERT_EQ(S.rows(),1u);
    EXPECT_NEAR(S(0,0),Pi/8.0,1e-10);
    auto nrm=Normalize(m,a);
    EXPECT_NEAR(nrm[0],1.0/std::sqrt(Pi/8.0),1e-10);
}
TEST(Mesh_Quadrature, WeightedOverlap_OneOverR)
{
    Mesh m=MakeProduct();
    ExpBasis a;
    auto V=WeightedOverlap(m,a,OneOverR());    // <a|1/r|a> = integral exp(-4r)/r d^3r = pi/4
    EXPECT_NEAR(V(0,0),Pi/4.0,1e-7);           // 1/r is non-smooth at the origin -> not machine-precision
}
TEST(Mesh_Quadrature, KineticGrad2)
{
    Mesh m=MakeProduct();
    ExpBasis a;
    auto K=KineticGrad2(m,a);                  // <grad a|grad a> = integral 4 exp(-4r) d^3r = pi/2
    EXPECT_NEAR(K(0,0),Pi/2.0,1e-10);
}

// Projection of a scalar field onto the basis: integral f a_i d^3r.  With a=f=exp(-2r),
// p_0 = integral exp(-4r) d^3r = pi/8.
TEST(Mesh_Quadrature, ScalarProjection)
{
    Mesh m=MakeProduct();
    ExpBasis a;
    auto p=Overlap(m,a,ExpScalar());           // ExpScalar = exp(-2r); a = exp(-2r)
    ASSERT_EQ(p.size(),1u);
    EXPECT_NEAR(p[0],Pi/8.0,1e-10);
}

//================================================================================================
//  6. Two-basis (rectangular) overlap: <a|a> off the Hermitian path equals pi/8.
//================================================================================================
TEST(Mesh_Quadrature, RectangularOverlap)
{
    Mesh m=MakeProduct();
    ExpBasis a;
    auto S=Overlap(m,a,a);
    ASSERT_EQ(S.rows(),1u);
    ASSERT_EQ(S.columns(),1u);
    EXPECT_NEAR(S(0,0),Pi/8.0,1e-10);
}

//================================================================================================
//  7. Angular-rule DEGREE AUDIT.  Each rule claims an algebraic degree L (it must integrate every
//  spherical polynomial of degree <= L exactly).  The hand-coded Lebedev-style tables ("Gauss")
//  were historically only smoke-checked (sum W, z^2), so this audits EVERY even monomial
//  x^{2a} y^{2b} z^{2c} with 2(a+b+c) <= L against the exact
//      integral x^{2a} y^{2b} z^{2c} dOmega = 4 pi (2a-1)!!(2b-1)!!(2c-1)!!/(2(a+b+c)+1)!!
//  and every odd monomial (degree <= L) against 0.  Direction constants in the 24/30 tables are
//  only ~7-figure, so the pass tolerance is 1e-6 relative to 4 pi.
//================================================================================================
namespace
{
double DoubleFactorial(int n) {double f=1; for (int i=n; i>1; i-=2) f*=i; return f;}

// Max |quadrature - exact| over ALL monomials of degree <= L, scaled by 1/(4 pi).
double AngularMomentError(const AngularMesh& ang, int L)
{
    double worst=0;
    for (int a=0; a<=L; a++)
        for (int b=0; a+b<=L; b++)
            for (int c=0; a+b+c<=L; c++)
            {
                double q=0;
                for (size_t i=0; i<ang.size(); i++)
                {
                    const rvec3_t& d=ang.Dirs()[i];
                    q+=ang.W()[i]*std::pow(d.x,a)*std::pow(d.y,b)*std::pow(d.z,c);
                }
                double exact=0;
                if (a%2==0 && b%2==0 && c%2==0)
                    exact=FourPi*DoubleFactorial(a-1)*DoubleFactorial(b-1)*DoubleFactorial(c-1)
                         /DoubleFactorial(a+b+c+1);
                worst=std::max(worst, std::fabs(q-exact)/FourPi);
            }
    return worst;
}
} //anon

TEST(Mesh_Angular, LebedevTablesHaveClaimedDegree)
{
    // {numDir, claimed L} straight from the table comments in LebedevAngularMesh.C.
    const std::pair<int,int> rules[]={{1,0},{2,0},{6,1},{8,3},{12,5},{24,7},{30,8},{50,11}};
    for (auto [n,L] : rules)
    {
        qcMesh::MeshParams mp; mp.angular=AngularKind::Lebedev; mp.nAngular=n;
        double err=AngularMomentError(MakeAngular(mp),L);
        EXPECT_LT(err,1e-6) << "Lebedev table numDir=" << n << " claims degree " << L
                            << " but max moment error is " << err;
    }
}

TEST(Mesh_Angular, GaussLegendreHasClaimedDegree)
{
    for (int L : {5,11,17,23,29})
    {
        qcMesh::MeshParams mp; mp.angular=AngularKind::GaussLegendre; mp.nAngular=L;
        double err=AngularMomentError(MakeAngular(mp),L);
        EXPECT_LT(err,1e-12) << "GaussLegendre L=" << L << " max moment error " << err;
    }
}

// Euler-Maclaren is trapezoidal in a transformed theta: it has NO algebraic degree, only slow
// ALGEBRAIC convergence (measured: degree-3 moment error ~1e-3 at L=29, m=2 -- three orders worse
// than GaussLegendre at the same L).  Audit it for what it is: the error must DROP with L, and the
// best transform (m=2) must reach 1e-3 by L=29.  For high-accuracy angular work use GaussLegendre.
TEST(Mesh_Angular, EulerMaclarenLowMomentsConverge)
{
    for (int m : {1,2,3})
    {
        qcMesh::MeshParams mp; mp.angular=AngularKind::EulerMaclaren; mp.em_m=m;
        mp.nAngular=11; double e11=AngularMomentError(MakeAngular(mp),3);
        mp.nAngular=29; double e29=AngularMomentError(MakeAngular(mp),3);
        EXPECT_LT(e29,e11) << "EulerMaclaren m=" << m << " does not converge: "
                           << e11 << " -> " << e29;
        if (m==2) EXPECT_LT(e29,1e-3) << "EulerMaclaren m=2 L=29 degree-3 moment error " << e29;
    }
}
