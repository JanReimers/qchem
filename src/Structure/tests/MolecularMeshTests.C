// File: src/Structure/tests/MolecularMeshTests.C  Acceptance tests for the Becke molecular mesh.
//
// The mirror of the old MoleculeMesh bug: a homonuclear dimer at separation -> 0 must integrate to
// the SINGLE-atom result, and a well-separated dimer must be ADDITIVE (= 2 x the atom).
#include "gtest/gtest.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <memory>
#include <stdexcept>

import qchem.Structure;                 // Molecule, Atom, Structure
import qchem.UnitCell;                   // UnitCell (uniform periodic mesh)
import qchem.Mesh.Quadrature;          // qcMesh::Integrate (over qcMath ScalarFunction)
import qchem.Reporting;                 // report::Scope: the mesh files its fields one level into the caller's scope
import qchem.LatticePoints;   // LatticePointsInSphere (V-STRUCTURE-ISP F2: images off the cell)
import qchem.Structure.MeshFactory;   // MakeMesh, StructureMeshMaker<S> (V-MESH-FACTORY)
import qchem.Lattice_3D.Imposition;   // ResolveImposition: a crystal run's SymmetrySource (the strain gate)
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
    auto m1=MakeMesh(a, mp);

    Molecule dimer;
    dimer.Insert(Site{1, o});
    dimer.Insert(Site{1, o});          // right on top of each other (R_ab = 0)
    auto m2=MakeMesh(dimer, mp);

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
    dimer.Insert(Site{1, RA});
    dimer.Insert(Site{1, RB});
    auto m=MakeMesh(dimer, mp);

    TwoGauss f(RA,RB);
    EXPECT_NEAR(qcMesh::Integrate(m,f),2*Pi32,1e-3);
}

// A far second atom must not corrupt the integral of a function localised on the first.
TEST(MolecularMesh, FarAtomDoesNotCorrupt)
{
    qcMesh::MeshParams mp=Params();
    rvec3_t RA(-8,0,0), RB(8,0,0);

    Molecule dimer;
    dimer.Insert(Site{1, RA});
    dimer.Insert(Site{1, RB});
    auto m=MakeMesh(dimer, mp);

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
    auto m=MakeMesh(cell, mp);

    EXPECT_EQ(m.size(), size_t(8*8*8));
    EXPECT_NEAR(qcMesh::Integrate(m,One()), a*a*a, 1e-9);   // sum of weights = Omega
}

// A smooth cell-periodic integrand: the midpoint rule is exact -> Omega/2.
TEST(LatticeMesh, UniformCellIntegratesPeriodicCosine)
{
    const double a=5.0;
    UnitCell cell(a);
    qcMesh::MeshParams mp; mp.nUniform=6;
    auto m=MakeMesh(cell, mp);

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
    auto m=MakeMesh(cell, mp);

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
    auto m=MakeMesh(cell, mp);

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
        for (auto n : LatticePointsInSphere(cell, range))
        {
            rvec3_t L=cell.ToCartesian(rvec3_t(n.x,n.y,n.z));
            for (const Site& a : cell.Sites()) itsSites.push_back(a.R+L);
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
    mp.cellKind=qcMesh::CellKind::RecipeStated(qcMesh::UnitCellKind::Becke);
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
        mesh=MakeMesh(cell, BeckeParams());
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
    auto m=MakeMesh(cell, BeckeParams());

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
    auto u=MakeMesh(f.cell, um);         // uniform grid at the SAME total point count
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

// THE MESH REPORTS ITSELF, ONE LEVEL IN (OpenWork §3, user ruling 2026-11-02): the CALLER opens a report::Scope naming
// itself and the builder writes its fields there, choosing no key -- so two callers' meshes land in two places, never one
// overwriting the other, and a scope opened under another nests (purely relative).
TEST(MolecularMesh, TheMeshFilesItsFieldsOneLevelIntoTheCallersScope)
{
    Molecule n2;
    n2.Insert(Site{7, rvec3_t(-1.03,0,0)});
    n2.Insert(Site{7, rvec3_t( 1.04,0,0)});
    namespace rpt=qchem::report;
    rpt::Begin("MolecularMeshReport", "2026-11-02T00:00:00");
    const qcMesh::Mesh a=[&]{ rpt::Scope s("callerA"); return MakeMesh(n2, qcMesh::MeshParams{.nRadial=20, .angularDegree=5}); }();
    const qcMesh::Mesh b=[&]{ rpt::Scope o("outer"); rpt::Scope s("callerB");
                              return MakeMesh(n2, qcMesh::MeshParams{.nRadial=30, .angularDegree=11}); }();
    const auto doc=rpt::CurrentRunReport();       // COPY before End (End moves the doc out)
    rpt::End();
    EXPECT_EQ(doc["callerA"]["points"].get<long>(), long(a.size()));
    EXPECT_EQ(doc["callerA"]["nRadial"].get<int>(), 20);
    EXPECT_EQ(doc["callerA"]["natom"].get<long>(), 2);
    EXPECT_EQ(doc["outer"]["callerB"]["points"].get<long>(), long(b.size()));
    EXPECT_EQ(doc["outer"]["callerB"]["nRadial"].get<int>(), 30);
    EXPECT_FALSE(doc.contains("points")) << "the builder chose a key of its own, or wrote at the root";
}

// THE MOLECULAR BECKE MESH's WEIGHT RESPONSE (lane B 2026-11-22): for a field f held per point (here a smooth function of the point's
// owner offset u, so it rides with its atom), Pullback(m, f) is d(Σ_g f_g w_g)/dR_D, the points riding their owners.  Against the central
// difference of the rebuilt mesh at R_D ± h, on a bent triatomic (water's shape, unequal Z).  The partition is smooth: h = 1e-4, 1e-8.
TEST(MolecularMesh, TheWeightResponseIsTheFiniteDifferenceOfTheRebuiltMesh)
{
    qcMesh::MeshParams mp=Params();
    const std::vector<int>     Z{8, 1, 1};
    const std::vector<rvec3_t> R0{rvec3_t(0.0, 0.0, 0.2214), rvec3_t(0.0, 1.4309, -0.8857), rvec3_t(0.0, -1.4309, -0.8857)};
    auto build=[&](const std::vector<rvec3_t>& R){ Molecule m; for (size_t a=0;a<R.size();a++) m.Insert(Site{Z[a], R[a]}); return MakeMesh(m, mp); };
    const qcMesh::Mesh m0=build(R0);
    ASSERT_NE(m0.WeightResponse(), nullptr) << "the molecular mesh installs its weight response";
    ASSERT_EQ(m0.OwnerOffsets().size(), m0.size()) << "every point records its owner offset";
    ASSERT_EQ(m0.NSites(), R0.size());
    auto field=[](const qcMesh::Mesh& m)
    {
        rvec_t f(m.size());
        for (size_t g=0;g<f.size();g++) { const rvec3_t& u=m.OwnerOffsets()[g]; f[g]=std::exp(-0.3*(u*u))*(1.0+0.4*u.x-0.2*u.y+0.7*u.z); }
        return f;
    };
    const rvec_t f0=field(m0);
    const auto r=m0.WeightResponse()->Pullback(m0, f0);
    ASSERT_EQ(r.gradient.size(), R0.size());
    auto E=[&](size_t a, int k, double s)
    {
        std::vector<rvec3_t> R=R0; (k==0 ? R[a].x : k==1 ? R[a].y : R[a].z)+=s;
        const qcMesh::Mesh m=build(R);
        EXPECT_EQ(m.size(), m0.size()) << "the moved mesh kept the same points";
        const rvec_t f=field(m);
        double e=0.0; for (size_t g=0;g<f.size();g++) e+=f[g]*m.Weights()[g];
        return e;
    };
    const double h=1e-4;
    double big=0.0, worst=0.0;
    rvec3_t sum(0,0,0);
    for (size_t a=0;a<R0.size();a++)
    {
        sum+=r.gradient[a];
        for (int k=0;k<3;k++)
        {
            const double fd=(E(a,k,h)-E(a,k,-h))/(2*h);
            const double an=k==0 ? r.gradient[a].x : k==1 ? r.gradient[a].y : r.gradient[a].z;
            std::printf("  [Becke weight response] atom %zu %c  analytic % .10e  FD % .10e  diff %.2e\n", a, "xyz"[k], an, fd, an-fd);
            big=std::max(big, std::fabs(fd)); worst=std::max(worst, std::fabs(an-fd));
        }
    }
    EXPECT_GT(big, 1e-3) << "a vanishing response cannot show it is wired";
    EXPECT_LT(worst, 1e-8*std::max(1.0,big)) << "the weight response is not the finite difference of the rebuilt mesh";
    EXPECT_LT(norm(sum), 1e-10) << "translation invariance";
}

// THE STRUCTURE MESH FACTORY (V-MESH-FACTORY §8 step 2).  A StructureMeshMaker<S> handed to a library must build EXACTLY what the
// factory overload builds for the same structure and recipe -- same points, same order, same weights, same per-site blocks -- or a
// library that receives a maker instead of building through the structure would move every number it integrates.
namespace
{
void ExpectSameMesh(const qcMesh::Mesh& a, const qcMesh::Mesh& b)
{
    ASSERT_EQ(a.size(), b.size());
    ASSERT_EQ(a.NSites(), b.NSites());
    for (size_t s=0; s<a.NSites(); s++) EXPECT_EQ(a.SiteBegin(s), b.SiteBegin(s));
    for (size_t q=0; q<a.size(); q++)
    {
        EXPECT_EQ(a.Weights()[q], b.Weights()[q]);
        EXPECT_EQ(norm(a.Points()[q]-b.Points()[q]), 0.0);
    }
}
} //anon

TEST(MeshFactory, AMakerBuildsWhatTheFactoryBuildsForEveryStructureKind)
{
    const qcMesh::MeshParams mp=Params();
    auto atom=std::make_shared<const Atom>(8, rvec3_t(0.1,0.2,0.3));
    auto mol=std::make_shared<Molecule>();
    mol->Insert(Site{8, rvec3_t(0,0,0)}); mol->Insert(Site{1, rvec3_t(0,1.4,1.1)}); mol->Insert(Site{1, rvec3_t(0,-1.4,1.1)});
    auto cell=std::make_shared<UnitCell>(6.0);
    cell->AddAtom(1, rvec3_t(0,0,0));
    cell->AddAtom(1, rvec3_t(0.5,0.5,0.5));
    qcMesh::MeshParams uni; uni.nUniform=6;

    ExpectSameMesh(StructureMeshMaker<Atom>(atom).Make(mp), MakeMesh(*atom, mp));
    ExpectSameMesh(StructureMeshMaker<Molecule>(mol).Make(mp), MakeMesh(*mol, mp));
    ExpectSameMesh(StructureMeshMaker<UnitCell>(cell).Make(uni), MakeMesh(*cell, uni));
    ExpectSameMesh(StructureMeshMaker<UnitCell>(cell).Make(BeckeParams()), MakeMesh(*cell, BeckeParams()));
    // Every ATOM-CENTRED mesh carries one site block per centre (the factory checks it on every overload, NDEBUG-proof).
    EXPECT_EQ(MakeMesh(*atom, mp).NSites(), size_t(1));
    EXPECT_EQ(MakeMesh(*mol, mp).NSites(), size_t(3));
    EXPECT_EQ(MakeMesh(*cell, BeckeParams()).NSites(), size_t(2));
}

TEST(MeshFactory, AMakerDemandsItsStructure)
{
    EXPECT_THROW(StructureMeshMaker<Molecule>(nullptr), std::invalid_argument);
}

// THE STRAIN GATE (controller 2026-11-22, the user's "Does the grid builder hold state information?"): a symmetry source supplies OPS only, and
// the MAKER owns the geometry.  An imposition resolved at the REFERENCE cell, handed to a maker over a STRAINED cell (an isotropic strain keeps
// the group, so the session keeps the imposition's ops), must build the strained cell's invariant mesh -- never the reference cell's.
TEST(MeshFactory, AStrainedCellsMakerBuildsOverItsOwnCellUnderTheReferencesOps)
{
    const auto bcc=[](double a){ auto c=std::make_shared<UnitCell>(a); c->AddAtom(1, rvec3_t(0,0,0)); c->AddAtom(1, rvec3_t(0.5,0.5,0.5)); return c; };
    const auto ref=bcc(6.0), strained=bcc(6.3);
    const CrystalImposition imp=ResolveImposition(Lattice_3D(*ref, ivec3_t(1,1,1)), {.imposeSymmetry=true});
    const std::vector<Symmetry::SymOp> ops=imp.MeshOps();
    ASSERT_FALSE(ops.empty());
    ASSERT_TRUE(imp.AveragesDensityOverOrbits());
    qcMesh::MeshParams uni; uni.nUniform=6;
    const qcMesh::FoldedMesh q=StructureMeshMaker<UnitCell>(strained).Make(uni, imp);
    ASSERT_TRUE(q.GetMesh());
    EXPECT_TRUE(q.HasFold()) << "an imposed run folds its mesh";
    ExpectSameMesh(*q.GetMesh(), MakeMesh(*strained, uni, ops));
    const qcMesh::Mesh atRef=MakeMesh(*ref, uni, ops);
    ASSERT_EQ(q.GetMesh()->size(), atRef.size());
    double moved=0.0;
    for (size_t i=0;i<atRef.size();i++) moved=std::max(moved, norm(q.GetMesh()->Points()[i]-atRef.Points()[i]));
    EXPECT_GT(moved, 0.1) << "the strained maker built the REFERENCE cell's mesh";
}

// The symmetric overload's per-structure hooks: an atom's invariant mesh is its radial-only set (an angular one refused); a molecule handed
// the identity alone gets its plain mesh, folded trivially (V-FOLD-NEUTRAL: every point its own orbit); and with NO ops neither hook is
// reached -- the plain mesh, unfolded.
TEST(MeshFactory, TheSymmetricOverloadsHooksPerStructure)
{
    auto atom=std::make_shared<const Atom>(8, rvec3_t(0,0,0));
    auto mol=std::make_shared<Molecule>();
    mol->Insert(Site{8, rvec3_t(0,0,0)}); mol->Insert(Site{1, rvec3_t(0,1.4,1.1)}); mol->Insert(Site{1, rvec3_t(0,-1.4,1.1)});
    struct OneOp : qcMesh::SymmetrySource
    {
        std::vector<Symmetry::SymOp> MeshOps() const override {return {Symmetry::SymOp{Matrix3D<double>(1,0,0, 0,1,0, 0,0,1), rvec3_t(0,0,0)}};}
        bool AveragesDensityOverOrbits() const override {return true;}
    } oneOp;
    qcMesh::MeshParams radial=Params(); radial.angular=qcMesh::AngularKind::None; radial.angularDegree=0;
    const qcMesh::FoldedMesh a=StructureMeshMaker<Atom>(atom).Make(radial, oneOp);
    ExpectSameMesh(*a.GetMesh(), MakeMesh(*atom, radial));
    EXPECT_FALSE(a.HasFold()) << "one direction per shell: the trivial fold";
    EXPECT_THROW(StructureMeshMaker<Atom>(atom).Make(Params(), oneOp), std::invalid_argument) << "an angular set on an atom claims a symmetry it lacks";
    const qcMesh::FoldedMesh e=StructureMeshMaker<Molecule>(mol).Make(Params(), oneOp);
    // The identity's frame is the lab frame: the plain mesh's points and weights, less the zero-weight points a radial node at r = 0 stacks
    // on each nucleus (the symmetric build drops them: coincident points would make the fold ambiguous), and every point its own orbit.
    const qcMesh::Mesh plain=MakeMesh(*mol, Params());
    size_t j=0;
    for (size_t i=0; i<plain.size(); i++)
    {
        if (plain.Weights()[i]==0.0) continue;
        ASSERT_LT(j, e.size());
        EXPECT_EQ(e.GetMesh()->Weights()[j], plain.Weights()[i]);
        EXPECT_EQ(norm(e.GetMesh()->Points()[j]-plain.Points()[i]), 0.0);
        j++;
    }
    EXPECT_EQ(j, e.size());
    EXPECT_TRUE(e.HasFold());
    EXPECT_EQ(e.NumOrbits(), e.size()) << "under the identity alone every point is its own orbit";
    const qcMesh::FoldedMesh m=StructureMeshMaker<Molecule>(mol).Make(Params(), qcMesh::NoSymmetry());
    ExpectSameMesh(*m.GetMesh(), MakeMesh(*mol, Params()));
    EXPECT_FALSE(m.HasFold());
}

// ---- V-FOLD-NEUTRAL step 1: a molecule's mesh invariant under its point group (the D2h-type abelian subgroup) ---------------------------
namespace
{
//! The test recipe at the smallest Lebedev rule that carries a D2h frame (degree 9, 38 directions).
qcMesh::MeshParams SymParams() {qcMesh::MeshParams mp=Params(); mp.angularDegree=9; return mp;}
Matrix3D<double> Diag(double a, double b, double c) {return Matrix3D<double>(a,0,0, 0,b,0, 0,0,c);}
//! A rotation by \a angle about the unit axis \a n (Rodrigues).
Matrix3D<double> Rotation(rvec3_t n, double angle)
{
    n=n/norm(n);
    const double c=std::cos(angle), s=std::sin(angle), t=1.0-c;
    return Matrix3D<double>(t*n.x*n.x+c,     t*n.x*n.y-s*n.z, t*n.x*n.z+s*n.y,
                            t*n.x*n.y+s*n.z, t*n.y*n.y+c,     t*n.y*n.z-s*n.x,
                            t*n.x*n.z-s*n.y, t*n.y*n.z+s*n.x, t*n.z*n.z+c);
}
const Matrix3D<double> kQ=Rotation(rvec3_t(0.3,-0.7,0.5), 0.9);   // an arbitrary deck orientation

//! A molecule given by its atoms in its own symmetry frame, placed in the deck by \a Q and moved by \a shift.
std::shared_ptr<const Molecule> Placed(const std::vector<Site>& frameSites, const Matrix3D<double>& Q, const rvec3_t& shift=rvec3_t(0,0,0))
{
    auto m=std::make_shared<Molecule>();
    for (Site s : frameSites) { s.R=Q*s.R+shift; m->Insert(s); }
    return m;
}
//! The ops \a W (frame) as the deck sees them: Q W Q^T, tau = 0.
std::vector<Symmetry::SymOp> PlacedOps(const std::vector<Matrix3D<double>>& W, const Matrix3D<double>& Q)
{
    std::vector<Symmetry::SymOp> ops;
    for (const auto& w : W) ops.push_back({Q*w*Transpose(Q), rvec3_t(0,0,0)});
    return ops;
}
struct Ops : qcMesh::SymmetrySource
{
    std::vector<Symmetry::SymOp> ops; bool averages=true;
    std::vector<Symmetry::SymOp> MeshOps() const override {return ops;}
    bool AveragesDensityOverOrbits() const override {return averages;}
};
const std::vector<Site> kWater={{8, rvec3_t(0,0,0.12)}, {1, rvec3_t(0,1.43,-0.98)}, {1, rvec3_t(0,-1.43,-0.98)}};       // C2v: C2 along z, yz plane
const std::vector<Matrix3D<double>> kC2v={Diag(1,1,1), Diag(-1,-1,1), Diag(-1,1,1), Diag(1,-1,1)};
std::vector<Matrix3D<double>> D2h()
{
    std::vector<Matrix3D<double>> w;
    for (double a : {1.0,-1.0}) for (double b : {1.0,-1.0}) for (double c : {1.0,-1.0}) w.push_back(Diag(a,b,c));
    return w;
}
std::vector<Site> SF6()
{
    const double d=2.96;
    return {{16, rvec3_t(0,0,0)}, {9, rvec3_t(d,0,0)}, {9, rvec3_t(-d,0,0)}, {9, rvec3_t(0,d,0)}, {9, rvec3_t(0,-d,0)},
            {9, rvec3_t(0,0,d)}, {9, rvec3_t(0,0,-d)}};
}
// NH3: C3v; its D2h-type subgroup is Cs, the mirror x -> -x (an N-H bond in the yz plane).
std::vector<Site> NH3()
{
    std::vector<Site> s{{7, rvec3_t(0,0,0.21)}};
    for (int k=0;k<3;k++)
    {
        const double phi=M_PI/2+2*M_PI*k/3, r=1.77;
        s.push_back({1, rvec3_t(r*std::cos(phi), r*std::sin(phi), -0.49)});
    }
    return s;
}
//! An anisotropic test integrand that moves with the molecule: Gaussians on the atoms times a polynomial in the frame coordinates.
double Field(const rvec3_t& r, const Molecule& mol, const Matrix3D<double>& Q, const rvec3_t& shift)
{
    const rvec3_t u=Transpose(Q)*(r-shift);   // frame coordinates
    double g=0.0;
    for (const Site& s : mol.Sites()) { const rvec3_t d=r-s.R; g+=std::exp(-0.8*(d*d)); }
    return g*(1.0+0.3*u.x+0.2*u.y*u.z+0.1*u.x*u.x);
}
double Integral(const qcMesh::Mesh& m, const Molecule& mol, const Matrix3D<double>& Q, const rvec3_t& shift=rvec3_t(0,0,0))
{
    double sum=0.0;
    for (size_t i=0;i<m.size();i++) sum+=m.Weights()[i]*Field(m.Points()[i], mol, Q, shift);
    return sum;
}
//! Every orbit of \a q carries one weight, and the star-average is idempotent.
void ExpectAnOrbitConsistentFold(const qcMesh::FoldedMesh& q)
{
    ASSERT_TRUE(q.HasFold());
    rvec_t w=q.GetMesh()->Weights(), w2=w;
    q.Symmetrize(w2);
    for (size_t i=0;i<w.size();i++) EXPECT_NEAR(w2[i], w[i], 1e-10*std::abs(w[i])) << "the weights differ across the orbit of point " << i;
    rvec_t f(w.size());
    for (size_t i=0;i<f.size();i++) f[i]=std::sin(1.0+q.GetMesh()->Points()[i].x)*std::cos(q.GetMesh()->Points()[i].y);
    q.Symmetrize(f);
    rvec_t g=f;
    q.Symmetrize(g);
    for (size_t i=0;i<f.size();i++) EXPECT_NEAR(g[i], f[i], 1e-14) << "the star-average is a projector";
}
} //anon

// The deck's orientation no longer reaches the answer: water placed in an arbitrary orientation (and off the origin) integrates a
// co-moving field to the same number as water in its symmetry frame -- the plain lab-frame mesh does not.
TEST(MolecularPointGroupMesh, TheDecksOrientationDoesNotReachTheIntegral)
{
    const qcMesh::MeshParams mp=SymParams();
    const Matrix3D<double> I=Diag(1,1,1);
    auto aligned=Placed(kWater, I);
    const rvec3_t shift(0.4,-1.1,2.3);
    auto rotated=Placed(kWater, kQ, shift);
    Ops a, r; a.ops=PlacedOps(kC2v, I); r.ops=PlacedOps(kC2v, kQ);
    const qcMesh::FoldedMesh qa=StructureMeshMaker<Molecule>(aligned).Make(mp, a);
    const qcMesh::FoldedMesh qr=StructureMeshMaker<Molecule>(rotated).Make(mp, r);
    const double Ea=Integral(*qa.GetMesh(), *aligned, I), Er=Integral(*qr.GetMesh(), *rotated, kQ, shift);
    // An aligned deck keeps the lab frame: the plain mesh's point set (the mirror image H carries its partner's points, in its own order).
    EXPECT_NEAR(Integral(MakeMesh(*aligned, mp), *aligned, I), Ea, 1e-13*std::abs(Ea));
    ASSERT_EQ(qr.size(), qa.size());
    EXPECT_NEAR(Er, Ea, 1e-12*std::abs(Ea));
    const double Eplain=Integral(MakeMesh(*rotated, mp), *rotated, kQ, shift);
    EXPECT_GT(std::abs(Eplain-Ea), 1e-10) << "the lab-frame mesh on the rotated deck: the artifact this removes (else the gate above is vacuous)";
    ExpectAnOrbitConsistentFold(qa);
    ExpectAnOrbitConsistentFold(qr);
    EXPECT_LT(qr.NumOrbits(), qr.size()/2) << "C2v folds the mesh";
}

// SF6 (Oh) under its D2h subgroup, in an arbitrary orientation: the 8 ops fold the mesh, the S atom at the centre carries all of them.
TEST(MolecularPointGroupMesh, SF6UnderD2hInAnyOrientation)
{
    qcMesh::MeshParams mp=Params(); mp.angularDegree=17;
    auto mol=Placed(SF6(), kQ, rvec3_t(1.0,2.0,3.0));
    Ops s; s.ops=PlacedOps(D2h(), kQ);
    const qcMesh::FoldedMesh q=StructureMeshMaker<Molecule>(mol).Make(mp, s);
    ExpectAnOrbitConsistentFold(q);
    EXPECT_LT(q.NumOrbits(), q.size()/5) << "D2h's 8 ops: most orbits are full (short ones lie on the planes)";
    // An orientation, not a new rule: the plain mesh of SF6 in its own frame (where the lab frame IS the group's), less its zero weights.
    const qcMesh::Mesh plain=MakeMesh(*Placed(SF6(), Diag(1,1,1)), mp);
    size_t nonzero=0; for (size_t i=0;i<plain.size();i++) nonzero+=(plain.Weights()[i]!=0.0);
    EXPECT_EQ(q.size(), nonzero);
}

// NH3 (C3v) under its D2h-type subgroup Cs, rotated; and the C3 it also carries is REFUSED by name (the full group is the follow-up).
TEST(MolecularPointGroupMesh, NH3UnderCsAndItsC3Refused)
{
    const qcMesh::MeshParams mp=SymParams();
    auto mol=Placed(NH3(), kQ);
    Ops cs; cs.ops=PlacedOps({Diag(1,1,1), Diag(-1,1,1)}, kQ);
    ExpectAnOrbitConsistentFold(StructureMeshMaker<Molecule>(mol).Make(mp, cs));
    const Matrix3D<double> C3=Rotation(rvec3_t(0,0,1), 2*M_PI/3);
    Ops c3v; c3v.ops=PlacedOps({Diag(1,1,1), C3, C3*C3}, kQ);
    try { StructureMeshMaker<Molecule>(mol).Make(mp, c3v); ADD_FAILURE() << "a C3 was accepted"; }
    catch (const std::invalid_argument& e) { EXPECT_NE(std::string(e.what()).find("D2h type"), std::string::npos) << e.what(); }
}

// A free run gets the same invariant mesh, unfolded; a geometry that carries the ops only roughly, an op with a translation or a spin
// action, and a non-closed set are refused.
TEST(MolecularPointGroupMesh, FreeRunsAndRefusals)
{
    const qcMesh::MeshParams mp=SymParams();
    auto mol=Placed(kWater, kQ);
    Ops fold, free; fold.ops=free.ops=PlacedOps(kC2v, kQ); free.averages=false;
    const qcMesh::FoldedMesh qf=StructureMeshMaker<Molecule>(mol).Make(mp, fold), qn=StructureMeshMaker<Molecule>(mol).Make(mp, free);
    EXPECT_FALSE(qn.HasFold());
    ExpectSameMesh(*qn.GetMesh(), *qf.GetMesh());

    std::vector<Site> bent=kWater; bent[1].R.y+=1e-6;
    EXPECT_THROW(StructureMeshMaker<Molecule>(Placed(bent, kQ)).Make(mp, fold), std::invalid_argument) << "the geometry does not carry the mirror";
    Ops shifted=fold; shifted.ops[1].tau=rvec3_t(0.1,0,0);
    EXPECT_THROW(StructureMeshMaker<Molecule>(mol).Make(mp, shifted), std::invalid_argument) << "a point-group op has no translation";
    Ops magnetic=fold; magnetic.ops[1].sigma=Symmetry::SpinAction::Flip;
    EXPECT_THROW(StructureMeshMaker<Molecule>(mol).Make(mp, magnetic), std::invalid_argument) << "no molecular Shubnikov imposition";
    Ops open; open.ops=PlacedOps({Diag(1,1,1), Diag(-1,-1,1), Diag(-1,1,1)}, kQ);
    EXPECT_ANY_THROW(StructureMeshMaker<Molecule>(mol).Make(mp, open)) << "not a group";
}

// Which angular rules carry the frame's reflections: every Lebedev rule from degree 9 up (the molecular entry rung is 29); the 24- and
// 30-direction rules (degrees 7, 8) do not, and are REFUSED by name with the degree to state instead.
TEST(MolecularPointGroupMesh, EveryLebedevRuleFromDegree9CarriesTheFrame)
{
    auto mol=Placed({{8, rvec3_t(0,0,0)}}, Diag(1,1,1));
    Ops s; s.ops=PlacedOps(D2h(), Diag(1,1,1));
    for (int d : {9,11,15,17,19,21,23,29,31,35,53,65})
    {
        qcMesh::MeshParams mp=Params(); mp.angularDegree=d;
        EXPECT_NO_THROW(StructureMeshMaker<Molecule>(mol).Make(mp, s)) << "degree " << d;
    }
    for (int d : {7,8})
    {
        qcMesh::MeshParams mp=Params(); mp.angularDegree=d;
        try { StructureMeshMaker<Molecule>(mol).Make(mp, s); ADD_FAILURE() << "degree " << d << " accepted"; }
        catch (const std::invalid_argument& e) { EXPECT_NE(std::string(e.what()).find("degree 9"), std::string::npos) << e.what(); }
    }
}
