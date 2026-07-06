// File: src/Structure/tests/MolecularMeshTests.C  Acceptance tests for the Becke molecular mesh.
//
// The mirror of the old MoleculeMesh bug: a homonuclear dimer at separation -> 0 must integrate to
// the SINGLE-atom result, and a well-separated dimer must be ADDITIVE (= 2 x the atom).
#include "gtest/gtest.h"
#include <cmath>

import qchem.Structure;                 // Molecule, Atom, Structure
import qchem.UnitCell;                   // UnitCell (uniform periodic mesh)
import qchem.Mesh.Quadrature;          // qcMesh::Integrate, ScalarField
import qchem.Math;                      // Pi32, Pi
using namespace qchem;

// Note: this TU sees BOTH the old global `Mesh`/`MeshParams` (via qchem.Structure) and the new
// qcMesh::* -- so we qualify qcMesh:: throughout and never name a bare Mesh/MeshParams.

namespace
{
// Gaussian exp(-|r-c|^2); integral over R^3 = pi^{3/2}.
class GaussAt : public qcMesh::ScalarField<double>
{
    rvec3_t itsC;
public:
    explicit GaussAt(const rvec3_t& c) : itsC(c) {}
    double  operator()(const rvec3_t& r) const override {double m=norm(r-itsC); return std::exp(-m*m);}
    rvec3_t Gradient  (const rvec3_t&)   const override {return rvec3_t(0,0,0);}
};

// Sum of two unit Gaussians centred at a and b; integral = 2 pi^{3/2}.
class TwoGauss : public qcMesh::ScalarField<double>
{
    rvec3_t itsA, itsB;
public:
    TwoGauss(const rvec3_t& a, const rvec3_t& b) : itsA(a), itsB(b) {}
    double  operator()(const rvec3_t& r) const override
    { double ma=norm(r-itsA), mb=norm(r-itsB); return std::exp(-ma*ma)+std::exp(-mb*mb); }
    rvec3_t Gradient  (const rvec3_t&)   const override {return rvec3_t(0,0,0);}
};

qcMesh::MeshParams Params()
{
    qcMesh::MeshParams mp;
    mp.radial=qcMesh::RadialKind::MHL; mp.nRadial=40; mp.mhl_m=2; mp.mhl_alpha=2.0;
    mp.angular=qcMesh::AngularKind::Gauss; mp.nAngular=24;
    return mp;
}
} //anon

// The bug's mirror: two coincident atoms must integrate to exactly the single-atom result.
TEST(MolecularMesh, CoincidentDimerEqualsAtom)
{
    qcMesh::MeshParams mp=Params();
    rvec3_t o(0,0,0);

    Atom a(1,o);
    auto m1=a.CreateIntegrationMesh(mp);

    Molecule dimer;
    dimer.Insert(new Atom(1,o));
    dimer.Insert(new Atom(1,o));          // right on top of each other (R_ab = 0)
    auto m2=dimer.CreateIntegrationMesh(mp);

    GaussAt g(o);
    double I1=qcMesh::Integrate(m1,g);
    double I2=qcMesh::Integrate(m2,g);
    EXPECT_NEAR(I2,I1,1e-12);             // coincident dimer == single atom (no 0/0 corruption)
    EXPECT_NEAR(I1,Pi32,1e-4);            // and both reproduce integral exp(-r^2) d^3r = pi^{3/2}
}

// Additivity: a well-separated dimer integrates the sum of two atom-centred Gaussians to 2 pi^{3/2}.
TEST(MolecularMesh, SeparatedDimerIsAdditive)
{
    qcMesh::MeshParams mp=Params();
    rvec3_t RA(-8,0,0), RB(8,0,0);

    Molecule dimer;
    dimer.Insert(new Atom(1,RA));
    dimer.Insert(new Atom(1,RB));
    auto m=dimer.CreateIntegrationMesh(mp);

    TwoGauss f(RA,RB);
    EXPECT_NEAR(qcMesh::Integrate(m,f),2*Pi32,1e-3);
}

// A far second atom must not corrupt the integral of a function localised on the first.
TEST(MolecularMesh, FarAtomDoesNotCorrupt)
{
    qcMesh::MeshParams mp=Params();
    rvec3_t RA(-8,0,0), RB(8,0,0);

    Molecule dimer;
    dimer.Insert(new Atom(1,RA));
    dimer.Insert(new Atom(1,RB));
    auto m=dimer.CreateIntegrationMesh(mp);

    GaussAt g(RA);                        // localised on atom A only
    EXPECT_NEAR(qcMesh::Integrate(m,g),Pi32,1e-3);
}

// ---- The uniform periodic (UnitCell) mesh -----------------------------------------------------------
namespace
{
// f == 1 everywhere: its cell integral is exactly the cell volume Omega (validates the equal weights
// sum to Omega for ANY n).
class One : public qcMesh::ScalarField<double>
{
public:
    double  operator()(const rvec3_t&) const override {return 1.0;}
    rvec3_t Gradient  (const rvec3_t&) const override {return rvec3_t(0,0,0);}
};

// cos^2(2 pi x / a): a smooth cell-periodic field on a cubic cell of edge a.  Exact cell integral = Omega/2
// (cos^2 = 1/2 + 1/2 cos(4 pi x/a), and the cos term averages to zero over the cell).  The midpoint rule
// integrates it exactly for n >= 3 points per axis -- a strong check that the fractional-midpoint mapping
// and weights are right.
class CosSqX : public qcMesh::ScalarField<double>
{
    double itsA;
public:
    explicit CosSqX(double a) : itsA(a) {}
    double  operator()(const rvec3_t& r) const override { double c=std::cos(2*Pi*r.x/itsA); return c*c; }
    rvec3_t Gradient  (const rvec3_t&)   const override {return rvec3_t(0,0,0);}
};
} //anon

// Uniform mesh: the equal weights must sum to the cell volume (Omega = a^3 for a cubic cell).
TEST(LatticeMesh, UniformCellIntegratesConstantToVolume)
{
    const double a=5.0;
    UnitCell cell(a);
    qcMesh::MeshParams mp; mp.nUniform=8;
    auto m=cell.CreateIntegrationMesh(mp);

    EXPECT_EQ(m.size(), size_t(8*8*8));
    EXPECT_NEAR(qcMesh::Integrate(m,One()), a*a*a, 1e-9);   // sum of weights = Omega
}

// A smooth cell-periodic integrand: the midpoint rule is exact -> Omega/2.
TEST(LatticeMesh, UniformCellIntegratesPeriodicCosine)
{
    const double a=5.0;
    UnitCell cell(a);
    qcMesh::MeshParams mp; mp.nUniform=6;
    auto m=cell.CreateIntegrationMesh(mp);

    EXPECT_NEAR(qcMesh::Integrate(m,CosSqX(a)), 0.5*a*a*a, 1e-9);
}

// Item H: with a physical eCut set, nUniform is DERIVED from the Nyquist bound n = ceil(2 a sqrt(2 eCut)/pi)
// (a = longest edge, x2 for density bandwidth) and the manual nUniform is ignored.  The longest edge binds
// an isotropic n, so a cubic cell's three axes share it.
TEST(LatticeMesh, ECutDerivesNyquistDivisions)
{
    const double a=6.0, eCut=8.0;
    UnitCell cell(a);
    qcMesh::MeshParams mp; mp.nUniform=3; mp.eCut=eCut;   // nUniform deliberately too small -> must be ignored
    auto m=cell.CreateIntegrationMesh(mp);

    const int nExpect=int(std::ceil(2.0*a*std::sqrt(2.0*eCut)/Pi));
    EXPECT_EQ(m.size(), size_t(nExpect)*nExpect*nExpect);
    EXPECT_GT(nExpect, 3);                                 // and the manual nUniform=3 was NOT used
    EXPECT_NEAR(qcMesh::Integrate(m,One()), a*a*a, 1e-9);  // still a valid volume quadrature
}

// The derived grid is alias-free BY CONSTRUCTION for a field at the density bandwidth: cos^2(2pi x/a) has a
// component at G = 4pi/a (|G|^2/2 = 8pi^2/a^2 ~ 2.19 a.u. for a=6), well under the density bandwidth eCut
// resolves, so the midpoint rule integrates it exactly.
TEST(LatticeMesh, ECutGridIsAliasFreeForDensityBandwidth)
{
    const double a=6.0, eCut=8.0;
    UnitCell cell(a);
    qcMesh::MeshParams mp; mp.eCut=eCut;
    auto m=cell.CreateIntegrationMesh(mp);

    EXPECT_NEAR(qcMesh::Integrate(m,CosSqX(a)), 0.5*a*a*a, 1e-9);
}
