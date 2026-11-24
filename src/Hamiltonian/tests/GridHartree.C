// File: Hamiltonian/tests/GridHartree.C  The grid Hartree route: the direct Coulomb sum against a closed form.
//
// A Gaussian charge has a closed-form potential, q erf(sqrt(a) r)/r, so the sum over a mesh can be judged point by point
// without an SCF.  Tolerances are measured: the sum is a coarse rule, and these pin what it delivers, not an ideal.
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <utility>
#include <vector>
#include <map>
#include <memory>
#include <tuple>
#include <numbers>
import qchem.Mesh.CoulombOperator;
import qchem.Mesh.Internal.RadialPoissonCoulomb;   // the solve's adjoint, stage by stage (a unit test may cheat)
import qchem.Hamiltonian.Factory;
import qchem.Blaze;
import qchem.Mesh;
import qchem.Mesh.Radial;
import qchem.Mesh.Angular;
import qchem.Structure;
import qchem.StructureData;
import qchem.Types;
import qchem.Structure.MeshFactory;   // MakeMesh: the geometry's mesh, by overload on its concrete type (V-MESH-FACTORY step 4)

using namespace qchem;

namespace {
const qcMesh::MeshParams theMesh = {.radial=qcMesh::RadialKind::MHL, .nRadial=80, .mhl_m=3, .mhl_alpha=2.0,
                                    .angular=qcMesh::AngularKind::Lebedev, .angularDegree=29, .beckeOrder=2};
qcMesh::Mesh OneAtomMesh(const qcMesh::MeshParams& mp=theMesh)
{
    Molecule m;
    m.Insert(Site{8, Vector3D<double>(0,0,0)});
    return MakeMesh(m, mp);
}
std::shared_ptr<const qcMesh::SiteMultipoles> Sites(const qcMesh::Mesh& mesh, const qcMesh::MeshParams& mp, int lMax, const rvec3vec_t& centres={rvec3_t(0,0,0)})
{
    return std::make_shared<const qcMesh::SiteMultipoles>(mesh, centres, qcMesh::MakeRadial(mp), qcMesh::MakeAngular(mp), lMax);
}
double Rho(double a, double r2) {return std::pow(a/std::numbers::pi,1.5)*std::exp(-a*r2);}   // unit charge
}

// The charge-weighted potential is the self-energy, which the erf form gives in closed form: E = q^2 sqrt(2a/pi).
TEST(DirectSumCoulomb, GaussianChargeSelfEnergy)
{
    const qcMesh::Mesh mesh=OneAtomMesh();
    const double a=1.3;
    rvec_t q(mesh.size());
    for (size_t g=0; g<mesh.size(); g++)
    {
        const auto& r=mesh.Points()[g];
        q[g]=mesh.Weights()[g]*Rho(a,r.x*r.x+r.y*r.y+r.z*r.z);
    }
    const rvec_t v=qcMesh::MakeDirectSumCoulomb(mesh)->Potential(q);
    double selfEnergy=0.0, charge=0.0;
    for (size_t g=0; g<mesh.size(); g++) {selfEnergy+=q[g]*v[g]; charge+=q[g];}
    EXPECT_NEAR(charge, 1.0, 1e-8);
    EXPECT_NEAR(selfEnergy, std::sqrt(2*a/std::numbers::pi), 1e-2);
}

// The radial Poisson solve on the same charge, against the closed form at every point: a one-site density is all l=0.  The error
// is the radial rule's, so it must fall with the radial count.
TEST(RadialPoissonCoulomb, GaussianChargePotential)
{
    auto worstAt=[](int nRadial)
    {
        qcMesh::MeshParams mp=theMesh; mp.nRadial=nRadial;
        const qcMesh::Mesh mesh=OneAtomMesh(mp);
        const double a=1.3;
        rvec_t q(mesh.size());
        for (size_t g=0; g<mesh.size(); g++)
        {
            const auto& r=mesh.Points()[g];
            q[g]=mesh.Weights()[g]*Rho(a,r.x*r.x+r.y*r.y+r.z*r.z);
        }
        const rvec_t v=qcMesh::MakeRadialPoissonCoulomb(Sites(mesh,mp,4))->Potential(q);
        double worst=0.0;
        for (size_t g=0; g<mesh.size(); g++)
        {
            const double r=norm(mesh.Points()[g]);
            if (r==0.0) continue;   // the rule's zero-weight node
            if (!std::isfinite(v[g])) return 1e300;
            worst=std::max(worst,std::fabs(v[g]-std::erf(std::sqrt(a)*r)/r));
        }
        return worst;
    };
    const double w80=worstAt(80), w160=worstAt(160), w320=worstAt(320);
    EXPECT_LT(w80,5e-5);
    EXPECT_LT(w160,w80/4);
    EXPECT_LT(w320,w160/4);
}

// A Gaussian off the centre needs every l: the expansion about the origin must converge on the closed form at the points.
TEST(RadialPoissonCoulomb, OffCentreGaussianNeedsHigherL)
{
    const qcMesh::Mesh mesh=OneAtomMesh();
    const double a=2.0; const rvec3_t c(0.3,-0.2,0.4);
    rvec_t q(mesh.size());
    for (size_t g=0; g<mesh.size(); g++)
    {
        const rvec3_t d=mesh.Points()[g]-c;
        q[g]=mesh.Weights()[g]*Rho(a,d*d);
    }
    auto worstAt=[&](int l)
    {
        const rvec_t v=qcMesh::MakeRadialPoissonCoulomb(Sites(mesh,theMesh,l))->Potential(q);
        double worst=0.0;
        for (size_t g=0; g<mesh.size(); g+=13)
        {
            const double r=norm(mesh.Points()[g]-c);
            if (!std::isfinite(v[g])) return 1e300;
            worst=std::max(worst,std::fabs(v[g]-std::erf(std::sqrt(a)*r)/r));
        }
        return worst;
    };
    const double w2=worstAt(2), w10=worstAt(10);
    EXPECT_LT(w10,w2/20);
    EXPECT_LT(w10,1e-4);
}

// Linear in the charges, self term included.
TEST(DirectSumCoulomb, Linear)
{
    const qcMesh::Mesh mesh=OneAtomMesh();
    rvec_t q(mesh.size());
    for (size_t g=0; g<mesh.size(); g++) q[g]=mesh.Weights()[g]*std::exp(-double(g%7));
    const auto op=qcMesh::MakeDirectSumCoulomb(mesh);
    const rvec_t v=op->Potential(q), v2=op->Potential(2.0*q);
    for (size_t g=0; g<mesh.size(); g+=997) EXPECT_NEAR(v2[g], 2.0*v[g], 1e-10*std::fabs(v[g]));
    EXPECT_EQ(op->HartreeEnergy(7.0,q,v), 0.5*op->SelfRepulsion(q,v)) << "on points both legs are one number, so Tr(DV) is not needed";
}

// The energy the solve delivers, not its pointwise error.  The radial map r(t) inside each interval only shows here: a quadratic in log r
// passed every pointwise test and still put 2.7e-7 of E_H into TeF₆, 95 μHa, at any rung.
TEST(RadialPoissonCoulomb, GaussianSelfEnergyIsNotBiasedByTheRadialMap)
{
    qcMesh::MeshParams mp=theMesh; mp.nRadial=300;
    const qcMesh::Mesh mesh=OneAtomMesh(mp);
    const double a=1.3;
    rvec_t q(mesh.size());
    for (size_t g=0; g<mesh.size(); g++)
    {
        const auto& r=mesh.Points()[g];
        q[g]=mesh.Weights()[g]*Rho(a,r.x*r.x+r.y*r.y+r.z*r.z);
    }
    const rvec_t v=qcMesh::MakeRadialPoissonCoulomb(Sites(mesh,mp,0))->Potential(q);
    double selfEnergy=0.0;
    for (size_t g=0; g<mesh.size(); g++) selfEnergy+=q[g]*v[g];
    const double exact=std::sqrt(2*a/std::numbers::pi);
    EXPECT_NEAR(selfEnergy/exact, 1.0, 1e-8);
}

// The channels of a smooth two-centre field, carried back to the points.  The partition of unity is what makes the sites' shares
// add up, so this is the projection and the synthesis checked together.  The error falls with l and is measured, not ideal.
TEST(SiteMultipoles, ProjectionThenSynthesisReturnsTheField)
{
    Molecule m;
    const rvec3_t A(0,0,0), B(0,0,2.0);
    m.Insert(Site{8, Vector3D<double>(0,0,0)});
    m.Insert(Site{1, Vector3D<double>(0,0,2.0)});
    const qcMesh::Mesh mesh=MakeMesh(m, theMesh);
    rvec_t f(mesh.size()), q(mesh.size());
    for (size_t g=0; g<mesh.size(); g++)
    {
        const rvec3_t a=mesh.Points()[g]-A, b=mesh.Points()[g]-B;
        f[g]=Rho(0.8,a*a)+0.5*Rho(1.1,b*b);
        q[g]=mesh.Weights()[g]*f[g];
    }
    auto worstAt=[&](int l)
    {
        const auto sites=Sites(mesh,theMesh,l,{A,B});
        auto Q=sites->Moments(q);
        const auto& W=sites->ShellWeights();
        for (auto& site : Q)
            for (size_t i=0; i<W.size(); i++)
                for (size_t c=0; c<sites->NumLM(); c++) site[i*sites->NumLM()+c]/=W[i];
        const rvec_t back=sites->Synthesize(Q);
        double worst=0.0;
        for (size_t g=0; g<mesh.size(); g++)
            if (std::isfinite(back[g])) worst=std::max(worst,std::fabs(back[g]-f[g]));
            else return 1e300;
        return worst;
    };
    const double l4=worstAt(4), l10=worstAt(10);
    EXPECT_LT(l10,l4/3.0);
    EXPECT_LT(l10,5e-3);
}

// The multipole energy's second leg on two overlapping Gaussian charges, against the closed form.  This rung failed by a large margin when
// the leg synthesized rho-tilde onto the mesh; on each site's own shells its charge is the mesh's, so it lands with the first-order one.
TEST(SiteMultipoles, OwnShellSelfEnergyMatchesTheClosedForm)
{
    Molecule m;
    const double R=1.5, a=1.3, b=0.9;
    const rvec3_t A(0,0,0), B(0,0,R);
    m.Insert(Site{8, Vector3D<double>(0,0,0)});
    m.Insert(Site{1, Vector3D<double>(0,0,R)});
    const qcMesh::Mesh mesh=MakeMesh(m, theMesh);
    rvec_t q(mesh.size());
    for (size_t g=0; g<mesh.size(); g++)
    {
        const rvec3_t x=mesh.Points()[g]-A, y=mesh.Points()[g]-B;
        q[g]=mesh.Weights()[g]*(Rho(a,x*x)+Rho(b,y*y));
    }
    const auto sites=Sites(mesh,theMesh,10,{A,B});
    const auto Q=sites->Moments(q);
    const rvec_t metric=[&]
    {
        rvec_t d(sites->NumChannels());
        size_t k=0;
        for (size_t s=0; s<sites->NumSites(); s++) for (double w : sites->ShellWeights()) for (size_t c=0; c<sites->NumLM(); c++) d[k++]=w;
        return d;
    }();
    const rvec_t qf=sites->Flatten(Q);
    const auto op=qcMesh::MakeMultipoleCoulomb(sites,metric);
    const rvec_t v=op->Potential(qf);
    const double p=a*b/(a+b);
    const double exact=std::sqrt(2*a/std::numbers::pi)+std::sqrt(2*b/std::numbers::pi)+2*std::erf(std::sqrt(p)*R)/R;
    EXPECT_NEAR(op->SelfRepulsion(qf,v), exact, 2e-5);
    EXPECT_EQ(op->HartreeEnergy(7.0,qf,v), 7.0-0.5*op->SelfRepulsion(qf,v)) << "a multipole's energy has two legs: Tr(DV) is the density's";
}

// H = dE/dD needs v = K q with K symmetric, which forces rely on.  The direct sum is, by construction.  The radial solve expands each
// source about its own site, so its raw K_gh and K_hg differ across sites; the operator answers ½(K + Kᵀ) (user ruling 2026-11-23), so
// what is pinned is symmetry to roundoff, on two different charge sets.
namespace {
struct TwoChargeSets {qcMesh::Mesh mesh; rvec3vec_t centres; rvec_t q, p;};
TwoChargeSets ChargesOn(const std::vector<std::pair<int,rvec3_t>>& atoms, const qcMesh::MeshParams& mp)
{
    Molecule m;
    rvec3vec_t centres(atoms.size());
    for (size_t a=0; a<atoms.size(); a++)
    {
        m.Insert(Site{atoms[a].first, Vector3D<double>(atoms[a].second.x,atoms[a].second.y,atoms[a].second.z)});
        centres[a]=atoms[a].second;
    }
    qcMesh::Mesh mesh=MakeMesh(m, mp);
    rvec_t q(mesh.size()), p(mesh.size());
    for (size_t g=0; g<mesh.size(); g++)
    {
        double a=0.0, b=0.0;
        for (size_t s=0; s<centres.size(); s++)
        {
            const rvec3_t d=mesh.Points()[g]-centres[s];
            a+=Rho(1.1,d*d)+0.3*Rho(6.0,d*d);
            b+=(1.0+0.7*d.x-0.4*d.z)*Rho(0.7,d*d);
        }
        q[g]=mesh.Weights()[g]*a;
        p[g]=mesh.Weights()[g]*b;
    }
    return {std::move(mesh),std::move(centres),std::move(q),std::move(p)};
}
double Asymmetry(const qcMesh::CoulombOperator& K, const rvec_t& q, const rvec_t& p)
{
    const rvec_t Kq=K.Potential(q), Kp=K.Potential(p);
    double pKq=0.0, qKp=0.0;
    for (size_t g=0; g<q.size(); g++) {pKq+=p[g]*Kq[g]; qKp+=q[g]*Kp[g];}
    return pKq-qKp;
}
}
TEST(CoulombOperatorSymmetry, DirectSumIsSymmetricAndThePoissonAsymmetryIsSmall)
{
    const std::vector<std::pair<int,rvec3_t>> water={{8,rvec3_t(0,0,0)},{1,rvec3_t(0,1.431,1.107)},{1,rvec3_t(0,-1.431,1.107)}};
    qcMesh::MeshParams mp=theMesh; mp.nRadial=40; mp.angularDegree=17;
    const auto c=ChargesOn(water,mp);
    EXPECT_NEAR(Asymmetry(*qcMesh::MakeDirectSumCoulomb(c.mesh),c.q,c.p), 0.0, 1e-9);
    const auto f=ChargesOn(water,theMesh);
    EXPECT_NEAR(Asymmetry(*qcMesh::MakeRadialPoissonCoulomb(Sites(f.mesh,theMesh,14,f.centres)),f.q,f.p), 0.0, 1e-12)
        << "the raw K measured -2.6e-6 here; the symmetrized one is exact to roundoff";
}

// Kᵀ is built stage by stage, each the TRANSPOSE of its forward stage: <w, S V> = <Sᵀ w, V> for Synthesize, Moments and the radial solve,
// on random-signed inputs (every channel and point live).  A wrong adjoint is a wrong Fock, so each stage is pinned alone.
namespace {
double Dot(const qcMesh::SiteMultipoles::channels_t& A, const qcMesh::SiteMultipoles::channels_t& B)
{
    double s=0.0;
    for (size_t a=0; a<A.size(); a++) for (size_t c=0; c<A[a].size(); c++) s+=A[a][c]*B[a][c];
    return s;
}
double Dot(const rvec_t& a, const rvec_t& b) {double s=0.0; for (size_t g=0; g<a.size(); g++) s+=a[g]*b[g]; return s;}
double Hash(size_t k) {return std::sin(1.7*double(k)+0.3*std::cos(0.37*double(k)));}   // deterministic, signed, unstructured
}
TEST(SiteMultipoles, EachStageOfTheAdjointIsTheTransposeOfItsForwardStage)
{
    const std::vector<std::pair<int,rvec3_t>> atoms={{8,rvec3_t(0,0,0)},{1,rvec3_t(0.3,1.431,1.107)},{9,rvec3_t(-0.2,-1.6,0.9)}};
    qcMesh::MeshParams mp=theMesh; mp.nRadial=40; mp.angularDegree=17;
    const auto c=ChargesOn(atoms,mp);
    const auto sites=Sites(c.mesh,mp,8,c.centres);
    const size_t n=c.mesh.size(), nCh=sites->Radii().size()*sites->NumLM();
    rvec_t w(n), q(n);
    for (size_t g=0; g<n; g++) {w[g]=Hash(g); q[g]=Hash(g+n);}
    qcMesh::SiteMultipoles::channels_t V(sites->NumSites(),std::vector<double>(nCh)), G=V;
    for (size_t a=0; a<V.size(); a++) for (size_t k=0; k<nCh; k++) {V[a][k]=Hash(7*a*nCh+k); G[a][k]=Hash(3*a*nCh+k+11);}

    const double sy=Dot(w,sites->Synthesize(V)), sa=Dot(sites->SynthesizeAdjoint(w),V);
    EXPECT_NEAR(sy, sa, 1e-12*std::fabs(sy)) << "Synthesize";
    const double my=Dot(sites->Moments(q),G), ma=Dot(q,sites->MomentsAdjoint(G));
    EXPECT_NEAR(my, ma, 1e-12*std::fabs(my)) << "Moments";
    const auto solver=qcMesh::RadialPoissonSolver(sites);
    const double py=Dot(G,solver.Solve(V)), pa=Dot(solver.SolveAdjoint(G),V);
    EXPECT_NEAR(py, pa, 1e-12*std::fabs(py)) << "the radial solve";
}

// A multipole fit is the one combination refused on a periodic structure: a correctness limit, stated as one.  The direct sum runs
// over points, so it has no meaning for channels on any structure, and past 2 l = degree the channels are not orthonormal.
TEST(CoulombRouteAdmission, MultipoleIsRefusedOnAPeriodicStructureOnTheDirectSumAndPastItsBound)
{
    using namespace Hamiltonian;
    Molecule m;
    m.Insert(Site{8, Vector3D<double>(0,0,0)});
    EXPECT_TRUE(AdmitCoulomb(m, {.fit=CoulombFit::Multipole}, theMesh));
    EXPECT_TRUE(AdmitVxc(m, VxcFit::Multipole, theMesh));
    EXPECT_FALSE(AdmitCoulomb(m, {.fit=CoulombFit::Multipole, .solver=CoulombSolver::DirectSum}, theMesh));
    EXPECT_TRUE (AdmitCoulomb(m, {.fit=CoulombFit::Multipole, .lMax=14}, theMesh)) << "2 l = 28 <= degree 29";
    EXPECT_FALSE(AdmitCoulomb(m, {.fit=CoulombFit::Multipole, .lMax=15}, theMesh)) << "past it the channels alias";
    EXPECT_FALSE(AdmitCoulomb(m, {.lMax=15}, theMesh, VxcFit::Multipole)) << "the Vxc fit has the same channels";
    EXPECT_TRUE (AdmitCoulomb(m, {.fit=CoulombFit::Delta, .lMax=20}, theMesh)) << "a Poisson solve on points expands nothing in them";

    UnitCell cell(8.0); cell.AddAtom(14, {0.5,0.5,0.5});
    const auto c=AdmitCoulomb(cell, {.fit=CoulombFit::Multipole}, theMesh);
    ASSERT_FALSE(c);
    EXPECT_NE(c.Error().find("conditionally convergent"), std::string::npos) << c.Error();
    const auto v=AdmitVxc(cell, VxcFit::Multipole, theMesh);
    ASSERT_FALSE(v);
    EXPECT_NE(v.Error().find("conditionally convergent"), std::string::npos) << v.Error();
    EXPECT_TRUE(AdmitCoulomb(cell, {.fit=CoulombFit::Delta}, theMesh)) << "the refusal is the multipole's alone";
}

// THE OPERATOR LEG OF THE δ HARTREE FORCE (lane S 2026-11-23): MakeSelfRepulsionPullback is d(qᵀK(R)q)/dR_A at FIXED q, the points riding
// their owners and the mesh's sites (shells, directions) moving with them.  Against the central difference of the REBUILT mesh and sites at
// R_A ± h, the same q handed to each.  A point is matched by its IDENTITY (its owner and its offset u from it), not its index: the builder
// drops a point whose Becke weight underflows to zero, so a move can change the count by a few.  Such a point's charge is w q ≈ 0 in both
// meshes, so a point present in only one carries zero.  An off-symmetric triatomic (unequal Z, no mirror) at the δ default lMax 16 on
// 80/L29, so every component is live.  Translation invariance: the gradient sums to zero.
TEST(RadialPoissonCoulomb, TheSelfRepulsionPullbackIsTheFiniteDifferenceOfTheRebuiltOperator)
{
    const int lMax=16;
    const std::vector<int> Z{8, 1, 9};
    const std::vector<rvec3_t> R0{rvec3_t(0.1,-0.2,0.0), rvec3_t(0.3,1.5,1.0), rvec3_t(-0.4,-1.2,1.9)};
    auto build=[&](const std::vector<rvec3_t>& R)
    {
        Molecule m;
        for (size_t a=0;a<R.size();a++) m.Insert(Site{Z[a], Vector3D<double>(R[a].x,R[a].y,R[a].z)});
        return MakeMesh(m, theMesh);
    };
    const qcMesh::Mesh m0=build(R0);
    rvec3vec_t c0(R0.size()); for (size_t a=0;a<R0.size();a++) c0[a]=R0[a];
    rvec_t q(m0.size());
    for (size_t g=0; g<m0.size(); g++)
    {
        double a=0.0;
        for (size_t s=0; s<R0.size(); s++) {const rvec3_t d=m0.Points()[g]-R0[s]; a+=(1.0+0.3*d.x-0.2*d.y)*Rho(1.1,d*d)+0.3*Rho(5.0,d*d);}
        q[g]=m0.Weights()[g]*a;
    }
    auto key=[](const qcMesh::Mesh& m, size_t a, size_t g)
    {
        const rvec3_t& u=m.OwnerOffsets()[g];
        return std::tuple<size_t,long long,long long,long long>(a, std::llround(u.x*1e8), std::llround(u.y*1e8), std::llround(u.z*1e8));
    };
    ASSERT_EQ(m0.OwnerOffsets().size(), m0.size());
    std::map<std::tuple<size_t,long long,long long,long long>,double> held;
    for (size_t a=0;a<m0.NSites();a++) for (size_t g=m0.SiteBegin(a);g<m0.SiteEnd(a);g++)
    {
        const auto [it,fresh]=held.emplace(key(m0,a,g), q[g]);
        if (!fresh) ASSERT_TRUE(q[g]==0.0 && it->second==0.0) << "only the rule's zero-weight r = 0 node repeats (once per direction)";
    }
    auto heldOn=[&](const qcMesh::Mesh& m)
    {
        rvec_t qm(m.size());
        for (size_t a=0;a<m.NSites();a++) for (size_t g=m.SiteBegin(a);g<m.SiteEnd(a);g++)
        {
            const auto it=held.find(key(m,a,g));
            qm[g]= it==held.end() ? 0.0 : it->second;
        }
        return qm;
    };
    const auto K0=qcMesh::MakeRadialPoissonCoulomb(Sites(m0,theMesh,lMax,c0));
    ASSERT_TRUE(K0->CanPullbackSelfRepulsion());
    EXPECT_FALSE(qcMesh::MakeDirectSumCoulomb(m0)->CanPullbackSelfRepulsion()) << "the direct sum's self term reads the weights: no pullback";
    const rvec3vec_t G=K0->MakeSelfRepulsionPullback(q);
    ASSERT_EQ(G.size(), R0.size());
    auto E=[&](size_t a, int k, double s)
    {
        std::vector<rvec3_t> R=R0; (k==0 ? R[a].x : k==1 ? R[a].y : R[a].z)+=s;
        const qcMesh::Mesh m=build(R);
        EXPECT_NEAR(double(m.size()), double(m0.size()), 10.0) << "the moved mesh kept its points, short of a few zero weights";
        rvec3vec_t c(R.size()); for (size_t b=0;b<R.size();b++) c[b]=R[b];
        const auto K=qcMesh::MakeRadialPoissonCoulomb(Sites(m,theMesh,lMax,c));
        const rvec_t qm=heldOn(m);
        return K->SelfRepulsion(qm, K->Potential(qm));
    };
    const double h=1e-4;
    double big=0.0, worst=0.0;
    rvec3_t sum(0,0,0);
    for (size_t a=0;a<R0.size();a++)
    {
        sum+=G[a];
        for (int k=0;k<3;k++)
        {
            const double fd=(E(a,k,h)-E(a,k,-h))/(2*h), an=(k==0 ? G[a].x : k==1 ? G[a].y : G[a].z);
            big=std::max(big,std::fabs(an)); worst=std::max(worst,std::fabs(an-fd));
            std::printf("[op pullback] atom %zu %c  analytic %+.10e  FD %+.10e  diff %+.2e\n", a, "xyz"[k], an, fd, an-fd);
        }
    }
    EXPECT_GT(big, 1e-2) << "the gate is vacuous";
    EXPECT_LT(worst, 1e-8*big) << "max |analytic - FD| " << worst << " of " << big << " (measured 8e-10 of 0.97)";
    EXPECT_LT(norm(sum), 1e-10*big) << "the cross-site pullback sums to zero";
}
