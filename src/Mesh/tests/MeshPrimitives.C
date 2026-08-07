// File: src/Mesh/tests/MeshPrimitives.C  Smoke + acceptance tests for the qcMesh primitives.
#include "gtest/gtest.h"
#include <cmath>
#include <stdexcept>
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
static void CheckAngular(const AngularMesh& a, double sumtol, double ztol)
{
    EXPECT_NEAR(Sum(a.W()),FourPi,sumtol);
    double zz=0;
    for (size_t j=0;j<a.size();j++) zz+=a.W()[j]*a.Dirs()[j].z*a.Dirs()[j].z;
    EXPECT_NEAR(zz,FourPi/3.0,ztol);
}
TEST(Mesh_Angular, Gauss)
{
    // numDir=32 was removed (inherited weight-sum bug, see LebedevAngularMesh.C).
    // 24 and 30 have ~7-figure direction constants (not unit vectors) -> exact only to ~1e-7.
    CheckAngular(LebedevAngular(12), 1e-10, 1e-10);
    CheckAngular(LebedevAngular(24), 1e-10, 1e-7);
    CheckAngular(LebedevAngular(30), 1e-10, 1e-7);
    CheckAngular(LebedevAngular(50), 1e-10, 1e-10);
}
TEST(Mesh_Angular, GaussLegendre)   // Gauss-exact in cos(theta) -> z^2 exact
{
    CheckAngular(MakeAngular({.angular=AngularKind::GaussLegendre, .angularDegree=11}), 1e-10, 1e-12);
}

//================================================================================================
//  4. ProductMesh + Integrate(ScalarField).
//================================================================================================
static Mesh MakeProduct()
{
    MeshParams p;
    p.radial=RadialKind::MHL; p.nRadial=200; p.mhl_m=2; p.mhl_alpha=3.0;
    p.angular=AngularKind::Lebedev; p.angularDegree=5;   // the 12-direction rule
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
    // Driven from LebedevMenu() -- the ONE place the ladder is written down -- so a table and a test
    // can no longer disagree.  (They had: two claims here were understated before the degrees were
    // MEASURED, and this test passed anyway because it only checks "no worse than claimed".)
    for (const auto& r : qcMesh::LebedevMenu())
    {
        if (r.nDir>50) continue;                       // the hand-written rules; canonical ones below
        const int n=r.nDir, L=r.degree;
        double err=AngularMomentError(LebedevAngular(n),L);
        EXPECT_LT(err,1e-6) << "Lebedev table numDir=" << n << " claims degree " << L
                            << " but max moment error is " << err;
    }
}

TEST(Mesh_Angular, CanonicalLebedevOrdersHaveClaimedDegree)
{
    // The canonical Lebedev-Laikov orders imported 2026-08-02 (16-digit constants, so the moment
    // audit is tight); the negative-weight orders 74/230/266 are excluded from the menu by the
    // generator audit.  This IS the transcription audit: full monomial sweep to the claimed
    // degree + positivity + sum w = 4pi.
    for (const auto& r : qcMesh::LebedevMenu())
    {
        if (r.nDir<86) continue;                       // the canonical Lebedev-Laikov orders
        const int n=r.nDir, L=r.degree;
        qcMesh::AngularMesh ang=LebedevAngular(n);
        ASSERT_EQ(int(ang.size()), n);
        double wsum=0;
        for (size_t i=0; i<ang.size(); i++)
        {
            EXPECT_GT(ang.W()[i], 0.0);
            wsum+=ang.W()[i];
            const rvec3_t& d=ang.Dirs()[i];
            EXPECT_NEAR(d.x*d.x+d.y*d.y+d.z*d.z, 1.0, 1e-14);
        }
        EXPECT_NEAR(wsum, FourPi, 1e-11);
        double err=AngularMomentError(ang,L);
        EXPECT_LT(err,1e-10) << "canonical Lebedev numDir=" << n << " claims degree " << L
                             << " but max moment error is " << err;
    }
}

TEST(Mesh_Angular, GaussLegendreHasClaimedDegree)
{
    for (int L : {5,11,17,23,29})
    {
        qcMesh::MeshParams mp; mp.angular=AngularKind::GaussLegendre; mp.angularDegree=L;
        double err=AngularMomentError(MakeAngular(mp),L);
        EXPECT_LT(err,1e-12) << "GaussLegendre L=" << L << " max moment error " << err;
    }
}

//================================================================================================
//  9. MEASURED algebraic exactness degree of every angular rule (R2.15 groundwork).
//
//  R2.15 wants the angular knob to become a polynomial DEGREE rather than "a count for Lebedev, an L
//  for the others".  That makes each rule's degree LOAD-BEARING: it stops being a comment and becomes
//  the interface contract, so a wrong annotation silently hands the caller a different grid than asked
//  for.  This block MEASURES the degree instead of trusting the annotation.
//
//  A rule is exact to degree D iff it integrates every monomial x^a y^b z^c with a+b+c <= D exactly
//  (monomials span the polynomials).  The exact sphere integral is closed-form:
//      integral x^a y^b z^c dOmega = 0                                        if any exponent is ODD
//                                  = 4pi (a-1)!!(b-1)!!(c-1)!! / (a+b+c+1)!!  if all are EVEN
//  so no spherical-harmonic machinery is needed and the test has no dependencies beyond the mesh.
//
//  This file's own header records why such a check is worth having: a 32-direction rule was shipped
//  VERBATIM from the old library with sum W = 0.971*4pi and had to be removed.  A degree measurement
//  would have caught it immediately.
//================================================================================================
//! Exact \int x^a y^b z^c dOmega over the unit sphere.  (DoubleFactorial -- n!! with (-1)!!=0!!=1 --
//! is already defined above for the radial moment checks.)
static double ExactMonomial(int a, int b, int c)
{
    if (a%2 || b%2 || c%2) return 0.0;          // odd power => antisymmetric => vanishes
    return FourPi*DoubleFactorial(a-1)*DoubleFactorial(b-1)*DoubleFactorial(c-1)
                 /DoubleFactorial(a+b+c+1);
}
//! The largest D such that the rule integrates EVERY monomial of degree <= D to within \a rtol
//! (relative to 4pi).  Stops at \a dMax; returns -1 if it cannot even integrate the constant.
static int MeasureDegree(const AngularMesh& m, double rtol, int dMax=40)
{
    for (int d=0; d<=dMax; d++)
        for (int a=d; a>=0; a--)
            for (int b=d-a; b>=0; b--)
            {
                const int c=d-a-b;
                double q=0.0;
                for (size_t j=0;j<m.size();j++)
                {
                    const rvec3_t& u=m.Dirs()[j];
                    q += m.W()[j]*std::pow(u.x,a)*std::pow(u.y,b)*std::pow(u.z,c);
                }
                if (std::fabs(q-ExactMonomial(a,b,c)) > rtol*FourPi) return d-1;
            }
    return dMax;
}

// Every rule in LebedevMenu() must DELIVER at least the degree it claims.
//
// The contract is EXPECT_GE, not EQ: R2.15's resolution rule is "round UP to a rule of at least the
// requested degree", so delivering more than claimed is fine and only under-delivery is a defect.  It
// also has to be GE for a measurement reason -- a monomial scan can report MORE than the constructed
// degree, because monomials odd under the rule's octahedral symmetry vanish identically on both sides
// and so never discriminate.  The measured degree and the SPECIAL ORBITS are printed for every rule:
// that table is what a same-degree tie has to be broken on (6 and 8 are the lowest-order instance).
TEST(Mesh_AngularDegree, LebedevOrdersMeetTheirClaimedDegree)
{
    for (const auto& r : qcMesh::LebedevMenu())
    {
        const double rtol = (r.nDir==24||r.nDir==30) ? 1e-6      // ~7-figure direction constants
                          : (r.nDir>=86)            ? 1e-10 : 1e-12;
        AngularMesh a=LebedevAngular(r.nDir);
        ASSERT_EQ(a.size(), size_t(r.nDir)) << "Lebedev " << r.nDir << ": wrong direction count";
        const int d=MeasureDegree(a,rtol);
        const qcMesh::SpecialOrbits o=qcMesh::ClassifyOrbits(a);
        std::cout << "[angular degree] Lebedev " << r.nDir << " dirs: claims " << r.degree
                  << ", measures " << d << ", orbits:"
                  << (o.axes100?" <100>":"") << (o.edges110?" <110>":"") << (o.corners111?" <111>":"")
                  << std::endl;
        EXPECT_GE(d, r.degree)
            << "Lebedev " << r.nDir << " directions UNDER-delivers: measured degree " << d
            << " < claimed " << r.degree << ".  Degree is the R2.15 interface -- a rule that claims more "
               "than it delivers silently under-integrates every caller that asks for that degree.";
    }
}

// The same-degree tie the user flagged: 6 and 8 both deliver degree 3 and are NOT interchangeable --
// they occupy DIFFERENT high-symmetry directions.  Pinned because the resolver picks the cheaper one,
// and that choice must be explainable as "cheapest at this degree", never as "the other is redundant".
TEST(Mesh_AngularDegree, EqualDegreeRulesCanDifferInOrbitDirection)
{
    const qcMesh::SpecialOrbits o6=qcMesh::ClassifyOrbits(LebedevAngular(6));
    const qcMesh::SpecialOrbits o8=qcMesh::ClassifyOrbits(LebedevAngular(8));
    EXPECT_TRUE (o6.axes100)    << "the 6-point rule should sit on the <100> axes";
    EXPECT_FALSE(o6.corners111);
    EXPECT_TRUE (o8.corners111) << "the 8-point rule should sit on the <111> body diagonals";
    EXPECT_FALSE(o8.axes100);
    EXPECT_EQ(qcMesh::ResolveLebedev(3).nDir, 6) << "degree 3 should resolve to the CHEAPER of the two";
}

// The resolver rounds UP across the menu's principled gaps (degrees 9/13/25/27 are absent: those orders
// carry negative weights, or shipped a weight-sum bug).  Asking for a gap degree must never return LESS
// exactness than requested.
TEST(Mesh_AngularDegree, ResolveRoundsUpAcrossTheGaps)
{
    for (int want=0; want<=35; want++)
    {
        const qcMesh::LebedevRule r=qcMesh::ResolveLebedev(want);
        EXPECT_GE(r.degree, want) << "ResolveLebedev(" << want << ") returned degree " << r.degree;
    }
    EXPECT_EQ(qcMesh::ResolveLebedev(25).nDir, 302) << "degree 25 falls in a gap -> the degree-29 rule";
    EXPECT_EQ(qcMesh::ResolveLebedev(29).nDir, 302);
    EXPECT_THROW(qcMesh::ResolveLebedev(36), std::runtime_error) << "past the menu must THROW, not clamp";
}

// Gauss-Legendre(cos theta) x uniform(phi): algebraically exact, so its L IS a degree.
TEST(Mesh_AngularDegree, GaussLegendreLIsADegree)
{
    for (int L : {5, 11, 17, 29})
    {
        AngularMesh a=MakeAngular({.angular=AngularKind::GaussLegendre, .angularDegree=L});
        EXPECT_GE(MeasureDegree(a,1e-10), L)
            << "GaussLegendre L=" << L << " is not exact to degree L -- its knob is documented as a degree.";
    }
}

// The R2.15 headline, measured rather than asserted: Lebedev reaches a given degree with far fewer
// directions than the Gauss-Legendre product grid, which is why the default flip is worth making.
TEST(Mesh_AngularDegree, LebedevIsCheaperThanGaussLegendreAtEqualDegree)
{
    AngularMesh leb=MakeAngular({.angular=AngularKind::Lebedev,       .angularDegree=29});
    AngularMesh gl =MakeAngular({.angular=AngularKind::GaussLegendre, .angularDegree=29});
    const int dLeb=MeasureDegree(leb,1e-10), dGL=MeasureDegree(gl,1e-10);
    EXPECT_GE(dLeb, 29);
    EXPECT_GE(dGL , 29);
    EXPECT_LT(leb.size(), gl.size())
        << "Leb-302 (" << leb.size() << " dirs, degree " << dLeb << ") vs GL-29 ("
        << gl.size() << " dirs, degree " << dGL << ")";
    std::cout << "[angular degree] equal degree ~29: Lebedev " << leb.size() << " dirs vs GaussLegendre "
              << gl.size() << " dirs (" << (100*leb.size())/gl.size() << "%)" << std::endl;
}
