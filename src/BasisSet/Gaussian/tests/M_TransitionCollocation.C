// File: BasisSet/Gaussian/tests/M_TransitionCollocation.C  The (k+q, k) transition collocation and its adjoint
// (doc/LinearResponsePlan.md §3d B1, R3 step 2).
//
// THREE GATES, each an independent route to the same object:
//   1. q = 0 with a HERMITIAN δD, it IS the ground-state collocation -- CollocateDensity under the geometry-only
//      screener (the transition entry points screen on geometry alone, by construction).  Run on an uncontracted
//      fixture, on a CONTRACTED basis (dzvp: the ground-state walk vs the transition's primitive-pair cubes), and
//      through the SPHERICAL view (the T congruence).
//   2. at q != 0 the scatter and the gather are EXACT adjoints: Σ conj(δD_ij) h_ij = Σ_L w_L Σ_g conj(u_L) v_L.
//   3. the periodic part u against EXPLICIT Bloch sums: u(r) = e^{-iq.r} Σ_ij δD_ij Φ^{k+q}_i(r) conj(Φ^k_j(r)),
//      with Φ from BlochPointValues (the lattice sum at points -- no collocation machinery at all), at a non-TRIM
//      k and q on a single-level ladder (so every pair lands on the one grid and u is pointwise).
#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <cstdio>
#include <functional>
#include <memory>
#include <vector>

import qchem.BasisSet.Gaussian.Point.Factory;     // Factory(BasisSetData, cell, engine, angular)
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD;   // NR_Evaluator::ContractCubeOverride (D-CUBE0)
import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;     // LatticeCollocation, TransitionCollocation
import qchem.BasisSet.Gaussian.Lattice.LatticeScreener;  // GeometryOnlyScreener, CollocationEps
import qchem.BasisSet.Gaussian.Lattice.SphericalLatticeView;   // MakeSphericalLatticeView (the GPW_SPHERICAL block)
import qchem.BasisSet;
import qchem.BasisSet.Orbital_1E_IBS;
import qchem.UnitCell;
import qchem.Structure;
import qchem.Types;
import qchem.Blaze;

using namespace qchem;
using namespace qchem::BasisSet::Gaussian;

namespace {

constexpr double kTwoPi=6.283185307179586476925286766559;

struct Periodic
{
    std::shared_ptr<const BasisSet::Real_BS> bs;
    const LatticeCollocation*    lc=nullptr;
    const TransitionCollocation* tc=nullptr;
    size_t n=0;
};
//! \a spherical: the SPHERICAL LATTICE VIEW over the Cartesian block (what a GPW_SPHERICAL run uses -- the NiO
//! recipe), made the way M_Spherical makes it; the Factory's own Spherical angular is the molecular family.
Periodic Make(BasisSetData data, const UnitCell& cell, bool spherical)
{
    Periodic p;
    std::shared_ptr<const BasisSet::Real_BS> cart(Factory(data, &cell, Engine::MnD, Angular::Cartesian));
    p.bs = spherical ? std::shared_ptr<const BasisSet::Real_BS>(PG_Spherical::MakeSphericalLatticeView(cart)) : cart;
    for (auto ibs : const_cast<BasisSet::Real_BS&>(*p.bs).Iterate<BasisSet::Real_OIBS>())
    {
        p.lc=dynamic_cast<const LatticeCollocation*>(ibs);
        p.tc=dynamic_cast<const TransitionCollocation*>(ibs);
        p.n =ibs->GetNumFunctions();
        break;
    }
    return p;
}
UnitCell SiCell()
{
    UnitCell c(Matrix3D<double>(0.0,5.13,5.13, 5.13,0.0,5.13, 5.13,5.13,0.0));   // FCC primitive, a = 10.26
    c.AddAtom(14, {0.0 ,0.0 ,0.0 });
    c.AddAtom(14, {0.25,0.25,0.25});
    return c;
}
UnitCell MnOCell()   // rocksalt MnO, FCC primitive: Mn d shells, so the SPHERICAL view's T really drops contaminants
{
    UnitCell c(Matrix3D<double>(0.0,4.20,4.20, 4.20,0.0,4.20, 4.20,4.20,0.0));
    c.AddAtom(25, {0.0,0.0,0.0});
    c.AddAtom( 8, {0.5,0.5,0.5});
    return c;
}
UnitCell WaterBox()   // a contracted (dzvp) basis in a periodic box -- the collocation does not care about conditioning
{
    UnitCell c(Matrix3D<double>(9.0,0.0,0.0, 0.0,9.5,0.0, 0.0,0.0,10.0));
    c.AddAtom(8, {0.50,0.50,0.50});
    c.AddAtom(1, {0.60,0.58,0.50});
    c.AddAtom(1, {0.42,0.59,0.50});
    return c;
}
LatticeSum1E::cellphase_t PhaseOf(const rvec3_t& k)
{ return [k](const ivec3_t& n){ const double a=kTwoPi*(k.x*n.x+k.y*n.y+k.z*n.z); return dcmplx(std::cos(a),std::sin(a)); }; }

mat_t<dcmplx> RandomMatrix(size_t n, unsigned long long seed, bool hermitian)
{
    unsigned long long r=seed;
    auto next=[&]{ r^=r<<13; r^=r>>7; r^=r<<17; return double(r%2000)/1000.0-1.0; };
    mat_t<dcmplx> M(n,n);
    for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) M(i,j)=0.05*dcmplx(next(),next());
    if (hermitian)
        for (size_t i=0;i<n;i++) { M(i,i)=dcmplx(M(i,i).real(),0.0); for (size_t j=i+1;j<n;j++) M(j,i)=std::conj(M(i,j)); }
    return M;
}

// The two-level ladder the gates run on (the second level a factor 4 coarser in Ecut, as the production ladder).
struct Ladder { std::vector<ivec3_t> N; std::vector<double> ecut; };
Ladder Two() { return {{ivec3_t(24,24,24), ivec3_t(12,12,12)}, {30.0, 7.5}}; }
Ladder One() { return {{ivec3_t(20,20,20)}, {30.0}}; }

//! GATE 1: q = 0, Hermitian δD == the ground-state collocation (geometry-only screen), at a non-TRIM k.
void QZeroIsTheGroundState(const Periodic& p, const UnitCell& cell, const char* what)
{
    ASSERT_TRUE(p.lc && p.tc) << what << ": the block does not realise both collocation faces";
    const Ladder L=Two();
    const rvec3_t k(0.25,0.0,0.125);
    const mat_t<dcmplx> dD=RandomMatrix(p.n, 88172645463325252ULL, true);
    chmat_t D(p.n);
    for (size_t i=0;i<p.n;i++) for (size_t j=i;j<p.n;j++) D(i,j)=dD(i,j);
    const GeometryOnlyScreener screen(CollocationEps());
    const std::vector<rvec_t>  rho=p.lc->CollocateDensity(D, PhaseOf(k), cell, L.N, L.ecut, screen);
    const std::vector<cvec_t>  u  =p.tc->CollocateTransition(dD, PhaseOf(k), rvec3_t(0,0,0), cell, L.N, L.ecut);
    ASSERT_EQ(rho.size(), u.size());
    double d=0, sc=0, im=0;
    for (size_t l=0;l<u.size();l++)
        for (size_t g=0;g<u[l].size();g++)
        {
            d =std::max(d , std::abs(u[l][g].real()-rho[l][g]));
            im=std::max(im, std::abs(u[l][g].imag()));
            sc=std::max(sc, std::abs(rho[l][g]));
        }
    std::printf("  [transition q=0] %-28s max|u-rho| %.2e  max|Im u| %.2e  (max|rho| %.2e)\n", what, d, im, sc);
    ASSERT_GT(sc, 1e-6) << what;
    EXPECT_LT(d , 1e-10*sc) << what << ": q=0 transition collocation != ground-state collocation";
    EXPECT_LT(im, 1e-10*sc) << what << ": a Hermitian δD at q = 0 must collocate a REAL density";
}

} // namespace

TEST(M_TransitionCollocation, QZero_EqualsTheGroundState_Uncontracted)
{
    const UnitCell cell=SiCell();
    QZeroIsTheGroundState(Make(BasisSetData::SIPP_SR, cell, false), cell, "Si SIPP_SR cartesian");
}
TEST(M_TransitionCollocation, QZero_EqualsTheGroundState_Spherical)
{
    const UnitCell cell=MnOCell();
    const Periodic c=Make(BasisSetData::VALENCE_LOWQ_VA, cell, false), v=Make(BasisSetData::VALENCE_LOWQ_VA, cell, true);
    ASSERT_LT(v.n, c.n) << "the spherical view must drop Cartesian contaminants, or this gate does not test T";
    QZeroIsTheGroundState(v, cell, "MnO VA spherical view");
}
TEST(M_TransitionCollocation, QZero_EqualsTheGroundState_Contracted)
{
    const UnitCell cell=WaterBox();
    QZeroIsTheGroundState(Make(BasisSetData::DZVP, cell, false), cell, "H2O dzvp (contracted)");
}

//! GATE 2: at q != 0 the scatter and the gather are EXACT adjoints, on the two-level ladder.
TEST(M_TransitionCollocation, ScatterGather_AreExactAdjoints)
{
    const UnitCell cell=MnOCell();
    for (bool spherical : {false, true})
    {
        const Periodic p=Make(BasisSetData::VALENCE_LOWQ_VA, cell, spherical);
        ASSERT_TRUE(p.tc);
        const Ladder L=Two();
        const rvec3_t k(0.25,0.0,0.125), q(0.5,0.25,0.0);
        const mat_t<dcmplx> dD=RandomMatrix(p.n, 1181783497276652981ULL, false);
        const std::vector<cvec_t> u=p.tc->CollocateTransition(dD, PhaseOf(k), q, cell, L.N, L.ecut);
        std::vector<cvec_t> v(u.size());
        unsigned long long r=0x9E3779B97F4A7C15ULL;
        for (size_t l=0;l<u.size();l++)
        {
            v[l]=cvec_t(u[l].size());
            for (size_t g=0;g<v[l].size();g++)
            { r^=r<<13; r^=r>>7; r^=r<<17; const double a=double(r%2000)/1000.0-1.0;
              r^=r<<13; r^=r>>7; r^=r<<17; const double b=double(r%2000)/1000.0-1.0; v[l][g]=dcmplx(a,b); }
        }
        const mat_t<dcmplx> h=p.tc->IntegrateTransition(v, PhaseOf(k), q, cell, L.N, L.ecut);
        dcmplx lhs(0.0), rhs(0.0);
        for (size_t l=0;l<u.size();l++)
        {
            const double w=cell.GetCellVolume()/double(u[l].size());
            for (size_t g=0;g<u[l].size();g++) lhs+=w*std::conj(u[l][g])*v[l][g];
        }
        for (size_t i=0;i<p.n;i++) for (size_t j=0;j<p.n;j++) rhs+=std::conj(dD(i,j))*h(i,j);
        const double sc=std::abs(lhs);
        std::printf("  [transition adjoint] %s  |<u,v> - Tr(dD^H h)| = %.3e  (rel %.2e)\n",
                    spherical ? "spherical view" : "cartesian", std::abs(lhs-rhs), std::abs(lhs-rhs)/sc);
        ASSERT_GT(sc, 1e-8) << "the pairing must be non-trivial";
        EXPECT_LT(std::abs(lhs-rhs), 1e-12*sc) << "the transition scatter and gather are not adjoint";
    }
}

//! GATE 3: u against EXPLICIT Bloch sums at a non-TRIM k and q (a single-level ladder: u is then pointwise).
TEST(M_TransitionCollocation, PeriodicPart_EqualsExplicitBlochSums)
{
    const UnitCell cell=SiCell();
    const Periodic p=Make(BasisSetData::SIPP_SR, cell, false);
    ASSERT_TRUE(p.lc && p.tc);
    const Ladder L=One();
    const rvec3_t k(0.25,0.0,0.125), q(0.5,0.25,0.0), kq=k+q;
    const mat_t<dcmplx> dD=RandomMatrix(p.n, 7640891576956012809ULL, false);
    const std::vector<cvec_t> u=p.tc->CollocateTransition(dD, PhaseOf(k), q, cell, L.N, L.ecut);
    const ivec3_t N=L.N[0];
    blazem::VecBuilder<rvec3_t> ptb; std::vector<size_t> idx;
    for (int a : {0, 3, 7, 11, 17}) for (int b : {1, 9, 14}) for (int c : {2, 10, 19})
    {
        push_back(ptb, cell.ToCartesian(rvec3_t(double(a)/N.x, double(b)/N.y, double(c)/N.z)));
        idx.push_back((size_t(a)*N.y+b)*N.z+c);
    }
    const rvec3vec_t pts=ptb.take();
    mat_t<dcmplx> Pk, Pkq;
    p.lc->BlochPointValues(pts, PhaseOf(k),  cell, Pk);
    p.lc->BlochPointValues(pts, PhaseOf(kq), cell, Pkq);
    double d=0, sc=0;
    for (size_t t=0;t<pts.size();t++)
    {
        dcmplx drho(0.0);
        for (size_t i=0;i<p.n;i++) for (size_t j=0;j<p.n;j++) drho+=dD(i,j)*Pkq(t,i)*std::conj(Pk(t,j));
        const rvec3_t f=cell.ToFractional(pts[t]);
        const double a=-kTwoPi*(q.x*f.x+q.y*f.y+q.z*f.z);
        const dcmplx ref=drho*dcmplx(std::cos(a),std::sin(a));
        d =std::max(d , std::abs(u[0][idx[t]]-ref));
        sc=std::max(sc, std::abs(ref));
    }
    std::printf("  [transition vs Bloch sums] %zu points  max|u - e^{-iqr} drho| %.2e  (max|u| %.2e)\n", pts.size(), d, sc);
    ASSERT_GT(sc, 1e-6);
    EXPECT_LT(d, 1e-9*sc) << "the transition collocation disagrees with explicit Bloch sums";
}


//! D-CUBE0: the production collocate (CollocateDensity) and gather (IntegratePotential) must give the SAME numbers
//! on the default separable-contraction route and on the reference box walk (GPW_CONTRACT_CUBE=0).  Until this test
//! the walk was only ever checked kernel-by-kernel in M_PG_BoxWalk; nothing exercised the production routing, so a
//! regression in the walk arm (the oracle a future change is judged against) could sit unseen.  Both routes are eps-
//! converged answers to one sum, so they agree at the screening tier, not bitwise (no per-component screen on the cube).
TEST(M_TransitionCollocation, WalkAndContractionRoutesAgree)
{
    using Evaluators::PG_Cart_MnD::NR_Evaluator;
    const UnitCell cell=SiCell();
    const Ladder L=Two();
    const rvec3_t k(0.25,0.0,0.125);
    const GeometryOnlyScreener screen(CollocationEps());

    std::vector<rvec_t> rho[2]; chmat_t h[2]; size_t n=0;
    for (int route=0; route<2; route++)                       // 0 = contraction (default), 1 = reference walk
    {
        NR_Evaluator::ContractCubeOverride()=(route==0);
        const Periodic p=Make(BasisSetData::SIPP_SR, cell, false);   // a FRESH evaluator per route: no memo can leak across
        ASSERT_TRUE(p.lc);
        n=p.n;
        const mat_t<dcmplx> Dm=RandomMatrix(n, 88172645463325252ULL, true);
        chmat_t D(n);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)=Dm(i,j);
        rho[route]=p.lc->CollocateDensity(D, PhaseOf(k), cell, L.N, L.ecut, screen);
        std::vector<rvec_t> V(rho[route].size());             // an arbitrary smooth-ish field: the density itself
        for (size_t l=0;l<V.size();l++) V[l]=rho[route][l];
        h[route]=p.lc->IntegratePotential(V, PhaseOf(k), cell, L.N, L.ecut, screen);
    }
    NR_Evaluator::ContractCubeOverride()=std::nullopt;

    double dRho=0, sRho=0;
    for (size_t l=0;l<rho[0].size();l++)
        for (size_t g=0;g<rho[0][l].size();g++)
        { dRho=std::max(dRho, std::abs(rho[0][l][g]-rho[1][l][g])); sRho=std::max(sRho, std::abs(rho[0][l][g])); }
    double dH=0, sH=0;
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++)
    { dH=std::max(dH, std::abs(dcmplx(h[0](i,j))-dcmplx(h[1](i,j)))); sH=std::max(sH, std::abs(dcmplx(h[0](i,j)))); }
    std::printf("  [walk vs contraction] max|drho|/max|rho| = %.2e   max|dh|/max|h| = %.2e\n", dRho/sRho, dH/sH);
    ASSERT_GT(sRho, 1e-6);
    ASSERT_GT(sH, 1e-8);
    EXPECT_LT(dRho, 1e-10*sRho) << "the two collocation routes disagree on the density";
    EXPECT_LT(dH,   1e-10*sH)  << "the two gather routes disagree on the potential matrix";
}

// ---- D-ENV step 5/6a: the typed GPW tolerances; the environment is NOT a way to set them ----
#include <cstdlib>
import qchem.BasisSet.Gaussian.Lattice.GPWTolerances;
import qchem.Environment;
TEST(GPWTolerances, DefaultsAreTodaysBehaviourAndTheRetiredEnvironmentIsIgnoredAndReported)
{
    using namespace qchem::BasisSet::Gaussian;
    GPWTolerances t;
    EXPECT_EQ(t.vlocEps, 1e-5);
    EXPECT_EQ(t.localPPRelCutoff, 30.0);
    EXPECT_NEAR(t.relFieldSharp, 1.0/3.0, 1e-15);
    EXPECT_TRUE(t.mgridEcuts.empty());
    EXPECT_EQ(t.screenEps, 1e-10);  EXPECT_EQ(t.densityEps, 1e-10);  EXPECT_EQ(t.relCutoff, 0.0);
    EXPECT_NEAR(t.fieldSharp, 2.0/3.0, 1e-15);
    EXPECT_TRUE(t.Describe().empty()) << "defaults print nothing on the banner";
    EXPECT_TRUE(t == GPWTolerances{});
    t.vlocEps=1e-7; t.mgridEcuts={53.33,17.78};
    EXPECT_NE(t.Describe().find("vlocEps"), std::string::npos);
    EXPECT_NE(t.Describe().find("mgridEcuts=53.33,17.78"), std::string::npos);
    EXPECT_FALSE(t == GPWTolerances{});

    // 6a: the old variables are retired -- reported with the deck key that replaces them, and they change NOTHING (no ApplyEnvOverrides exists)
    unsetenv("GPW_SCREEN_EPS"); unsetenv("QCHEM_BECKE_NR");
    EXPECT_TRUE(qchem::RetiredEnvironmentSet().empty());
    setenv("GPW_SCREEN_EPS","1e-6",1); setenv("QCHEM_BECKE_NR","12",1);
    const auto set=qchem::RetiredEnvironmentSet();
    ASSERT_EQ(set.size(),2u);
    bool sawScreen=false, sawNR=false;
    for (const auto& r : set)
    {
        if (r.name=="GPW_SCREEN_EPS")  { sawScreen=true; EXPECT_EQ(r.deckKey,"solid.tolerances.screenEps"); }
        if (r.name=="QCHEM_BECKE_NR")  { sawNR=true;     EXPECT_EQ(r.deckKey,"solid.xcMesh.nRadial"); }
    }
    EXPECT_TRUE(sawScreen && sawNR);
    EXPECT_EQ(GPWTolerances{}.screenEps, 1e-10) << "a retired variable must not reach the typed default";
    unsetenv("GPW_SCREEN_EPS"); unsetenv("QCHEM_BECKE_NR");
}

// Option B: a tolerance handed to the molecular evaluator through ApplyTolerances REACHES its pair loops (a looser analytic screen
// moves the lattice-summed overlap, a little), and restoring the default restores the bits -- the tolerance-dependent caches are dropped
// and rebuilt, none of it leaks across settings.
TEST(GPWTolerances, ApplyTolerancesReachesThePairLoopsAndIsReversible)
{
    const UnitCell cell=SiCell();
    const Periodic p=Make(BasisSetData::SIPP_SR, cell, false);
    const Periodic_Gaussian_IBS* pg=nullptr;
    for (auto ibs : const_cast<BasisSet::Real_BS&>(*p.bs).Iterate<BasisSet::Real_OIBS>()) { pg=dynamic_cast<const Periodic_Gaussian_IBS*>(ibs); break; }
    ASSERT_NE(pg, nullptr);
    const auto phase=PhaseOf(rvec3_t(0.25,0.0,0.125));
    auto diff=[&](const chmat_t& a, const chmat_t& b)
    { double d=0; for (size_t i=0;i<p.n;i++) for (size_t j=i;j<p.n;j++) d=std::max(d,std::abs(dcmplx(a(i,j))-dcmplx(b(i,j)))); return d; };
    const chmat_t S0=pg->MakeOverlap(phase,cell);
    GPWTolerances loose; loose.screenEps=1e-3;
    pg->ApplyTolerances(loose);
    const chmat_t S1=pg->MakeOverlap(phase,cell);
    EXPECT_GT(diff(S0,S1), 1e-9) << "a 1e-3 screen must drop terms the 1e-10 screen kept";
    EXPECT_LT(diff(S0,S1), 1e-1);
    pg->ApplyTolerances(GPWTolerances{});
    EXPECT_EQ(diff(S0,pg->MakeOverlap(phase,cell)), 0.0) << "restoring the default restores the bits";
}

// The guard: once a basis has built collocation work, a DIFFERENT tolerance is a broken invariant (two evaluators sharing one basis),
// so it throws; the same tolerance is a no-op.
TEST(GPWTolerances, ApplyTolerancesThrowsIfChangedAfterCollocationWork)
{
    const UnitCell cell=SiCell();
    const Periodic p=Make(BasisSetData::SIPP_SR, cell, false);
    const Periodic_Gaussian_IBS* pg=nullptr;
    for (auto ibs : const_cast<BasisSet::Real_BS&>(*p.bs).Iterate<BasisSet::Real_OIBS>()) { pg=dynamic_cast<const Periodic_Gaussian_IBS*>(ibs); break; }
    ASSERT_NE(pg, nullptr);
    GPWTolerances loose; loose.screenEps=1e-6;
    pg->ApplyTolerances(loose);                    // before any work: fine
    const Ladder L=One();
    chmat_t D(p.n);
    for (size_t i=0;i<p.n;i++) for (size_t j=i;j<p.n;j++) D(i,j)=(i==j);
    const GeometryOnlyScreener screen(loose.densityEps);
    (void)p.lc->CollocateDensity(D, PhaseOf(rvec3_t(0,0,0)), cell, L.N, L.ecut, screen);   // builds the task list
    EXPECT_NO_THROW(pg->ApplyTolerances(loose));   // equal: no-op
    EXPECT_THROW(pg->ApplyTolerances(GPWTolerances{}), std::runtime_error);
}
