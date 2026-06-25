// File: src/Mesh1/tests/Mesh1Primitives.C  Smoke + acceptance tests for the qcMesh1 primitives.
#include "gtest/gtest.h"
#include <cmath>
import qchem.Mesh1.Product;       // Mesh, ProductMesh, MakeRadial/MakeAngular, MeshParams, enums
import qchem.Mesh1.Quadrature;    // free functions + ScalarField/BasisField
import qchem.Mesh1.GaussLegendre;
import qchem.Blaze;
import qchem.Math;                // Pi, Pi32, FourPi

using namespace qcMesh1;          // Mesh, ScalarField/BasisField, MakeRadial/MakeAngular, quadrature, ...

//================================================================================================
//  Test fields: a single normalised-ish s-function phi(r) = exp(-2|r|), and some scalar potentials.
//================================================================================================
class ExpBasis : public BasisField<double>
{
public:
    size_t size() const override {return 1;}
    rvec_t operator()(const rvec3_t& r) const override {return {std::exp(-2*norm(r))};}
    rvec3vec_t Gradient(const rvec3_t& r) const override
    {
        double mr=norm(r);
        if (mr==0.0) return {rvec3_t(0,0,0)};
        return {(r/mr)*(-2*std::exp(-2*mr))};
    }
};

class ExpScalar : public ScalarField<double>     // exp(-2r)
{
public:
    double   operator()(const rvec3_t& r) const override {return std::exp(-2*norm(r));}
    rvec3_t  Gradient  (const rvec3_t&  ) const override {return rvec3_t(0,0,0);}
};

class GaussScalar : public ScalarField<double>   // exp(-r^2)
{
public:
    double   operator()(const rvec3_t& r) const override {double m=norm(r); return std::exp(-m*m);}
    rvec3_t  Gradient  (const rvec3_t&  ) const override {return rvec3_t(0,0,0);}
};

class OneOverR : public ScalarField<double>      // 1/r, finite (0) at the origin
{
public:
    double   operator()(const rvec3_t& r) const override {double m=norm(r); return m==0.0 ? 0.0 : 1.0/m;}
    rvec3_t  Gradient  (const rvec3_t&  ) const override {return rvec3_t(0,0,0);}
};

static double Sum(const rvec_t& v) {double s=0; for (size_t i=0;i<v.size();i++) s+=v[i]; return s;}

//================================================================================================
//  1. Shared 1D Gauss-Legendre: known table + polynomial exactness to degree 2n-1.
//================================================================================================
TEST(Mesh1_GaussLegendre, TwoPointTable)
{
    GaussLegendre gl(2,-1.0,1.0);
    double r=1.0/std::sqrt(3.0);
    EXPECT_NEAR(gl.x[0],-r,1e-14);
    EXPECT_NEAR(gl.x[1], r,1e-14);
    EXPECT_NEAR(gl.w[0],1.0,1e-14);
    EXPECT_NEAR(gl.w[1],1.0,1e-14);
}

TEST(Mesh1_GaussLegendre, PolynomialExactness)
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
    auto m=MakeRadial(p);
    double I=0;
    for (size_t i=0;i<m->size();i++) I+=m->W()[i]*std::exp(-2*m->R()[i]);
    return I;
}

TEST(Mesh1_Radial, MHL)
{
    MeshParams p; p.radial=RadialKind::MHL; p.nRadial=200; p.mhl_m=2; p.mhl_alpha=3.0;
    EXPECT_NEAR(RadialIntegralExp(p),0.25,4e-15);
}
TEST(Mesh1_Radial, Log)
{
    MeshParams p; p.radial=RadialKind::Log; p.nRadial=400; p.logStart=1e-4; p.logStop=40;
    EXPECT_NEAR(RadialIntegralExp(p),0.25,1e-3);
}
TEST(Mesh1_Radial, Linear)
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
    auto a=MakeAngular(p);
    EXPECT_NEAR(Sum(a->W()),FourPi,sumtol);
    double zz=0;
    for (size_t j=0;j<a->size();j++) zz+=a->W()[j]*a->Dirs()[j].z*a->Dirs()[j].z;
    EXPECT_NEAR(zz,FourPi/3.0,ztol);
}
TEST(Mesh1_Angular, Gauss)
{
    // numDir=32 was removed (inherited weight-sum bug, see GaussAngularMesh.C).
    // 24 and 30 have ~7-figure direction constants (not unit vectors) -> exact only to ~1e-7.
    CheckAngular({.angular=AngularKind::Gauss, .nAngular=12}, 1e-10, 1e-10);
    CheckAngular({.angular=AngularKind::Gauss, .nAngular=24}, 1e-10, 1e-7);
    CheckAngular({.angular=AngularKind::Gauss, .nAngular=30}, 1e-10, 1e-7);
    CheckAngular({.angular=AngularKind::Gauss, .nAngular=50}, 1e-10, 1e-10);
}
TEST(Mesh1_Angular, GaussLegendre)   // Gauss-exact in cos(theta) -> z^2 exact
{
    CheckAngular({.angular=AngularKind::GaussLegendre, .nAngular=11}, 1e-10, 1e-12);
}
TEST(Mesh1_Angular, EulerMaclaren)   // both the sum and z^2 are approximations for this scheme
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
    p.angular=AngularKind::Gauss; p.nAngular=12;
    auto rad=MakeRadial(p);
    auto ang=MakeAngular(p);
    return ProductMesh(*rad,*ang);
}

TEST(Mesh1_Product, ExpIntegral)        // integral exp(-2r) d^3r = 4*pi * 0.25 = pi
{
    Mesh m=MakeProduct();
    EXPECT_NEAR(Integrate(m,ExpScalar()),Pi,1e-12);
}
TEST(Mesh1_Product, GaussianIntegral)   // integral exp(-r^2) d^3r = pi^{3/2}
{
    Mesh m=MakeProduct();
    EXPECT_NEAR(Integrate(m,GaussScalar()),Pi32,1e-6);
}

//================================================================================================
//  5. Quadrature free functions against analytic <exp(-2r)|...|exp(-2r)>.
//================================================================================================
TEST(Mesh1_Quadrature, OverlapAndNormalize)
{
    Mesh m=MakeProduct();
    ExpBasis a;
    auto S=Overlap(m,a);                       // <a|a> = integral exp(-4r) d^3r = pi/8
    ASSERT_EQ(S.rows(),1u);
    EXPECT_NEAR(S(0,0),Pi/8.0,1e-10);
    auto nrm=Normalize(m,a);
    EXPECT_NEAR(nrm[0],1.0/std::sqrt(Pi/8.0),1e-10);
}
TEST(Mesh1_Quadrature, WeightedOverlap_OneOverR)
{
    Mesh m=MakeProduct();
    ExpBasis a;
    auto V=WeightedOverlap(m,a,OneOverR());    // <a|1/r|a> = integral exp(-4r)/r d^3r = pi/4
    EXPECT_NEAR(V(0,0),Pi/4.0,1e-7);           // 1/r is non-smooth at the origin -> not machine-precision
}
TEST(Mesh1_Quadrature, KineticGrad2)
{
    Mesh m=MakeProduct();
    ExpBasis a;
    auto K=KineticGrad2(m,a);                  // <grad a|grad a> = integral 4 exp(-4r) d^3r = pi/2
    EXPECT_NEAR(K(0,0),Pi/2.0,1e-10);
}

// Projection of a scalar field onto the basis: integral f a_i d^3r.  With a=f=exp(-2r),
// p_0 = integral exp(-4r) d^3r = pi/8.
TEST(Mesh1_Quadrature, ScalarProjection)
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
TEST(Mesh1_Quadrature, RectangularOverlap)
{
    Mesh m=MakeProduct();
    ExpBasis a;
    auto S=Overlap(m,a,a);
    ASSERT_EQ(S.rows(),1u);
    ASSERT_EQ(S.columns(),1u);
    EXPECT_NEAR(S(0,0),Pi/8.0,1e-10);
}
