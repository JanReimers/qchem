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
    mp.angular=qcMesh::AngularKind::Lebedev; mp.angularDegree=7;   // the 24-direction rule (was nAngular=24)
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

// The periodic-Becke test resolution: MHL nr=40 radial, GaussLegendre L=29 angular -- PINNED to GL
// deliberately, so these tests stay bit-stable across angular-scheme default changes (the GPW gate
// recipe qcMesh::BeckeXCParams flipped its default to the Lebedev tables 2026-08-17, R2.15; the
// canonical tables now reach degree 35, so the old "they stop at L=11" reason is history).
// (EulerMaclaren, also degree-less, was retired 2026-08-07.)  The fuzzy Voronoi switching shell is
// angular-quadrature limited (see Mesh_Angular audits).
qcMesh::MeshParams BeckeParams()
{
    qcMesh::MeshParams mp=Params();
    mp.cellKind=qcMesh::UnitCellKind::Becke;
    mp.angular=qcMesh::AngularKind::GaussLegendre;
    mp.angularDegree=29;
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

// ---- SITE BLOCKS: the partition an atom-centred mesh already carries (doc/OpenWork.md Step 0a) --------
//
// A Becke mesh IS a union of per-centre grids whose weights already fold in that centre's partition
// function w_A(r) -- that is the scheme's defining identity, Integral f = Sum_A Integral w_A f.  Recording
// which block came from which centre costs one index per atom and turns "the atomic moment/charge" from a
// point sample of a field into a genuine integral.  These three gates pin the contract.

// (1) STRUCTURE: one block per atom, contiguous, covering the whole mesh.
TEST(LatticeMesh, BeckeSiteBlocksPartitionTheMesh)
{
    const auto& f=FCCBecke::Get();
    ASSERT_EQ(f.mesh.NSites(), 2u) << "one site block per atom of the FCC diamond cell";
    EXPECT_EQ(f.mesh.SiteBegin(0), 0u);
    EXPECT_EQ(f.mesh.SiteEnd(f.mesh.NSites()-1), f.mesh.size());
    for (size_t a=0; a+1<f.mesh.NSites(); a++)
        EXPECT_EQ(f.mesh.SiteEnd(a), f.mesh.SiteBegin(a+1)) << "blocks must be contiguous, site " << a;
    for (size_t a=0; a<f.mesh.NSites(); a++)
        EXPECT_GT(f.mesh.SiteEnd(a), f.mesh.SiteBegin(a)) << "site " << a << " kept no points";
}

// (2) THE DEFINING IDENTITY: Sum_A Integral w_A f == Integral f, for a non-trivial f.  This is what makes a
// per-site integral a decomposition of the total rather than an independent (and arbitrary) quantity.
TEST(LatticeMesh, BeckeSiteIntegralsSumToTheTotal)
{
    const auto& f=FCCBecke::Get();
    const double alpha=1.0;
    LatticeGauss g(f.cell,alpha);
    rvec_t fv(f.mesh.size());
    for (size_t q=0; q<f.mesh.size(); q++) fv[q]=g(f.mesh.Points()[q]);

    rvec_t site=qcMesh::SiteIntegrals(f.mesh, fv);
    ASSERT_EQ(site.size(), f.mesh.NSites());
    double sum=0.0; for (size_t a=0;a<site.size();a++) sum+=site[a];
    // Exact to roundoff: both sides are the SAME sum of w_q f_q, only regrouped by block.
    EXPECT_NEAR(sum, qcMesh::Integrate(f.mesh,g), 1e-12*std::fabs(sum));
}

// (3) EQUIVALENT SITES GET EQUAL SHARES.  The two Si of the diamond cell are crystallographically
// equivalent, so each must own half the cell volume and half of any lattice-symmetric field.  This is the
// site-resolved check the aggregate "[Becke grid] ... dropped N tail pts" line cannot make -- and a
// site-DEPENDENT mesh defect is exactly what bit the MnO campaign (the moment died on the corner atom and
// survived on the centre one, 2026-08-11).  Quadrature-limited, not partition-limited.
TEST(LatticeMesh, BeckeEquivalentSitesOwnEqualShares)
{
    const auto& f=FCCBecke::Get();
    rvec_t ones(f.mesh.size(), 1.0);
    rvec_t vol=qcMesh::SiteIntegrals(f.mesh, ones);          // Integral w_A d3r = site A's share of Omega
    ASSERT_EQ(vol.size(), 2u);
    const double tot=vol[0]+vol[1];
    // THE ASSERTION IS THE SHARE, NOT THE VOLUME.  Site equality is a PARTITION property and holds to
    // ~1e-9 here; the absolute Integral w_A d3r sits 0.43% under Omega/2 because a CONSTANT is the worst
    // case for an atom-centred grid (its radial meshes must reach infinity and the partition tails are
    // eps-cut) -- that is quadrature accuracy, a different property, already pinned by
    // BeckeCellIntegratesConstantToVolume.  Testing the ratio keeps this gate on the thing it names.
    EXPECT_NEAR(vol[0]/tot, 0.5, 1e-8) << "crystallographically equivalent sites must own equal shares";
    EXPECT_NEAR(vol[1]/tot, 0.5, 1e-8);
    EXPECT_NEAR(tot/f.cell.GetCellVolume(), 1.0, 1e-2);      // loose: the constant-integration worst case
}
