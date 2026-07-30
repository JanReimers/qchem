// File: src/Structure/tests/MolecularMeshTests.C  Acceptance tests for the Becke molecular mesh.
//
// The mirror of the old MoleculeMesh bug: a homonuclear dimer at separation -> 0 must integrate to
// the SINGLE-atom result, and a well-separated dimer must be ADDITIVE (= 2 x the atom).
#include "gtest/gtest.h"
#include <cmath>
#include <cstdio>

import qchem.Structure;                 // Molecule, Atom, Structure
import qchem.UnitCell;                   // UnitCell (uniform periodic mesh)
import qchem.Mesh.Quadrature;          // qcMesh::Integrate (over qcMath ScalarFunction)
import qchem.Math;                      // Pi32, Pi
using namespace qchem;

// Note: this TU sees BOTH the old global `Mesh`/`MeshParams` (via qchem.Structure) and the new
// qcMesh::* -- so we qualify qcMesh:: throughout and never name a bare Mesh/MeshParams.

namespace
{
// Gaussian exp(-|r-c|^2); integral over R^3 = pi^{3/2}.
class GaussAt : public ScalarFunction<double>
{
    rvec3_t itsC;
public:
    explicit GaussAt(const rvec3_t& c) : itsC(c) {}
    double  operator()(const rvec3_t& r) const override {double m=norm(r-itsC); return std::exp(-m*m);}
    rvec3_t Gradient  (const rvec3_t&)   const override {return rvec3_t(0,0,0);}
};

// Sum of two unit Gaussians centred at a and b; integral = 2 pi^{3/2}.
class TwoGauss : public ScalarFunction<double>
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
class One : public ScalarFunction<double>
{
public:
    double  operator()(const rvec3_t&) const override {return 1.0;}
    rvec3_t Gradient  (const rvec3_t&) const override {return rvec3_t(0,0,0);}
};

// cos^2(2 pi x / a): a smooth cell-periodic field on a cubic cell of edge a.  Exact cell integral = Omega/2
// (cos^2 = 1/2 + 1/2 cos(4 pi x/a), and the cos term averages to zero over the cell).  The midpoint rule
// integrates it exactly for n >= 3 points per axis -- a strong check that the fractional-midpoint mapping
// and weights are right.
class CosSqX : public ScalarFunction<double>
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

// ---- The periodic Becke (UnitCell) mesh -------------------------------------------------------------------
namespace
{
// Lattice-summed Gaussians on the cell's atom sites: f(r) = sum_{a,L} exp(-alpha |r - R_a - L|^2), a smooth
// cell-periodic field whose exact cell integral is natom (pi/alpha)^{3/2} (each site contributes one full
// Gaussian per cell).  The image sum is epsilon-converged by magnitude (exp(-alpha d^2) < tiny).
class LatticeGauss : public ScalarFunction<double>
{
    const UnitCell& itsCell;
    double          itsAlpha;
    std::vector<rvec3_t> itsSites;   // all images within the significant range of every cell point
public:
    LatticeGauss(const UnitCell& cell, double alpha) : itsCell(cell), itsAlpha(alpha)
    {
        // Significant range from any point in the home cell: exp(-alpha d^2) >= 1e-16.  The home cell
        // (fractional [0,1)^3) reaches out to its farthest corner, so cover that plus the Gaussian reach.
        double rCorner=0;
        for (int c=0; c<8; c++)
            rCorner=std::max(rCorner, norm(cell.ToCartesian(rvec3_t(c&1,(c>>1)&1,(c>>2)&1))));
        double range=std::sqrt(37.0/itsAlpha)+rCorner;
        for (auto n : cell.CellsInSphere(range))
        {
            rvec3_t L=cell.ToCartesian(rvec3_t(n.x,n.y,n.z));
            for (auto a : cell) itsSites.push_back(a->itsR+L);
        }
    }
    double  operator()(const rvec3_t& r) const override
    {
        double f=0;
        for (const auto& s : itsSites) {double d2=norm(r-s); f+=std::exp(-itsAlpha*d2*d2);}
        return f;
    }
    rvec3_t Gradient  (const rvec3_t&) const override {return rvec3_t(0,0,0);}
};

// The periodic-Becke test resolution == the GPW Becke-XC-gate resolution (BeckeXCParams in
// GPW_SCF_UT.C): MHL nr=40 radial, GaussLegendre L=29 angular.  GaussLegendre, NOT the Lebedev
// tables (they stop at L=11) and NOT EulerMaclaren (no algebraic degree; converges slowly) -- the
// fuzzy Voronoi switching shell is angular-quadrature limited (see Mesh_Angular audits).
qcMesh::MeshParams BeckeParams()
{
    qcMesh::MeshParams mp=Params();
    mp.cellKind=qcMesh::UnitCellKind::Becke;
    mp.angular=qcMesh::AngularKind::GaussLegendre;
    mp.nAngular=29;
    return mp;
}
} //anon

// The Si-gate FCC diamond cell + its periodic Becke mesh, built ONCE and shared by the acceptance
// tests below (the ~45 s partition build is the expensive part; the integrals are cheap).
struct FCCBecke
{
    FCCUnitCell cell;
    qcMesh::Mesh mesh;
    FCCBecke() : cell(10.26)
    {
        cell.AddAtom(14, rvec3_t(0,0,0));
        cell.AddAtom(14, rvec3_t(0.25,0.25,0.25));
        mesh=cell.CreateIntegrationMesh(BeckeParams());
    }
    static const FCCBecke& Get() {static FCCBecke f; return f;}
};

// Partition-of-unity sanity: the lattice translates of the weighted per-atom grids tile the crystal, so a
// constant integrates to the cell volume Omega.  (This is quadrature-limited, not partition-limited: the
// cell functions are smooth but not band-limited, so the tolerance is the mesh's, not machine's --
// measured -2.9e-4 at nr=40/GL-29 on the sparse single-atom cell, the partition-shell worst case.)
TEST(LatticeMesh, BeckeCellIntegratesConstantToVolume)
{
    const double a=6.0;
    UnitCell cell(a);
    cell.AddAtom(1, rvec3_t(0,0,0));
    auto m=cell.CreateIntegrationMesh(BeckeParams());

    EXPECT_NEAR(qcMesh::Integrate(m,One())/(a*a*a), 1.0, 1e-3);
}

// Two-atom (diamond FCC) cell: a smooth lattice-summed Gaussian must integrate to natom (pi/alpha)^{3/2}.
TEST(LatticeMesh, BeckeFCCIntegratesLatticeGaussian)
{
    const auto& f=FCCBecke::Get();
    const double alpha=1.0, exact=2*std::pow(Pi/alpha,1.5);
    EXPECT_NEAR(qcMesh::Integrate(f.mesh,LatticeGauss(f.cell,alpha))/exact, 1.0, 1e-4);
}

// The reason the Becke grid exists: a SHARP atom-centred field (alpha=100, the two-scale integrand that
// explodes the uniform grid) is captured by the dense radial region at ordinary point counts.  The uniform
// mesh at the same total point count misses it badly.
TEST(LatticeMesh, BeckeCapturesSharpCoreField)
{
    const auto& f=FCCBecke::Get();
    const double alpha=100.0, exact=2*std::pow(Pi/alpha,1.5);
    LatticeGauss sharp(f.cell,alpha);
    EXPECT_NEAR(qcMesh::Integrate(f.mesh,sharp)/exact, 1.0, 1e-4);

    qcMesh::MeshParams um; um.nUniform=int(std::ceil(std::cbrt(double(f.mesh.size()))));
    auto u=f.cell.CreateIntegrationMesh(um);         // uniform grid at the SAME total point count
    EXPECT_GT(std::fabs(qcMesh::Integrate(u,sharp)/exact-1.0), 1e-2);
}

// Wrapped emission: every Becke point lies in the home cell (fractional coordinates in [0,1)).
TEST(LatticeMesh, BeckePointsAreWrappedIntoHomeCell)
{
    const auto& f=FCCBecke::Get();
    for (size_t q=0; q<f.mesh.size(); q++)
    {
        rvec3_t fr=f.cell.ToFractional(f.mesh.Points()[q]);
        EXPECT_GE(fr.x,-1e-12); EXPECT_LT(fr.x,1.0+1e-12);
        EXPECT_GE(fr.y,-1e-12); EXPECT_LT(fr.y,1.0+1e-12);
        EXPECT_GE(fr.z,-1e-12); EXPECT_LT(fr.z,1.0+1e-12);
    }
}
