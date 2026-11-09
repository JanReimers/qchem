// File: BasisSet/Gaussian/tests/M_CollocationPullback.C  THE PULLBACK OF THE GPW GATHER (doc/GeometryOptimizationPlan.md §13.2; lane B
// 2026-11-09): LatticeCollocation::IntegratePotentialPullback, the jet twin of the reference box walk, against the central finite
// difference of the gather it differentiates, Re Tr(D h[V]) with h = IntegratePotential(V, ...) at FIXED level fields V and FIXED D,
// the basis rebuilt on displaced atoms.
//
// THREE fixtures, each a different arm of the kernel: Si SIPP_SR (uncontracted Cartesian, FCC, off-origin atoms), H3 dzvp in a box
// (CONTRACTED shells: the radial's r²-derivative sums over primitives), MnO through the SPHERICAL lattice view (d shells, D expanded
// through T).  All at a NON-TRIM k (complex phases: a sign or conjugation slip cannot hide), on the two-level ladder.  The FD runs on
// the PRODUCTION gather (the contracted cube), the pullback on the walk: the two agree to the walk-vs-cube tier (1e-10 relative,
// M_TransitionCollocation.WalkAndContractionRoutesAgree), well inside this gate.
//
// And TRANSLATION: under a constant field every atom moved together leaves Tr(D h) unchanged up to the raster's egg-box, so the
// pullbacks sum to (almost) zero while each one is not.
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <functional>
#include <memory>
#include <vector>

import qchem.BasisSet.Gaussian.Point.Factory;     // Factory(BasisSetData, cell, engine, angular)
import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;     // LatticeCollocation
import qchem.BasisSet.Gaussian.Lattice.LatticeScreener;  // GeometryOnlyScreener, CollocationEps
import qchem.BasisSet.Gaussian.Lattice.SphericalLatticeView;   // MakeSphericalLatticeView
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

struct Atom { int Z; rvec3_t f; };
struct Fixture { const char* what; Matrix3D<double> A; std::vector<Atom> atoms; BasisSetData data; bool spherical; };

UnitCell CellOf(const Fixture& x, size_t moved=size_t(-1), const rvec3_t& dr=rvec3_t(0,0,0))
{
    UnitCell c(x.A);
    for (size_t a=0;a<x.atoms.size();a++)
    {
        rvec3_t f=x.atoms[a].f;
        if (a==moved || moved==size_t(-2)) f=c.ToFractional(c.ToCartesian(f)+dr);   // -2: every atom (the translation)
        c.AddAtom(x.atoms[a].Z, f);
    }
    return c;
}
struct Block { std::shared_ptr<const BasisSet::Real_BS> bs; const LatticeCollocation* lc=nullptr; size_t n=0; };
Block Make(const Fixture& x, const UnitCell& cell)
{
    Block b;
    std::shared_ptr<const BasisSet::Real_BS> cart(Factory(x.data, &cell, Engine::MnD, Angular::Cartesian));
    b.bs = x.spherical ? std::shared_ptr<const BasisSet::Real_BS>(PG_Spherical::MakeSphericalLatticeView(cart)) : cart;
    for (auto ibs : const_cast<BasisSet::Real_BS&>(*b.bs).Iterate<BasisSet::Real_OIBS>())
    {
        b.lc=dynamic_cast<const LatticeCollocation*>(ibs);
        b.n =ibs->GetNumFunctions();
        break;
    }
    return b;
}
LatticeSum1E::cellphase_t PhaseOf(const rvec3_t& k)
{ return [k](const ivec3_t& n){ const double a=kTwoPi*(k.x*n.x+k.y*n.y+k.z*n.z); return dcmplx(std::cos(a),std::sin(a)); }; }

chmat_t RandomHermitian(size_t n, unsigned long long seed)
{
    unsigned long long r=seed;
    auto next=[&]{ r^=r<<13; r^=r>>7; r^=r<<17; return double(r%2000)/1000.0-1.0; };
    chmat_t D(n);
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)= i==j ? dcmplx(0.1*next(),0.0) : 0.05*dcmplx(next(),next());
    return D;
}
struct Ladder { std::vector<ivec3_t> N; std::vector<double> ecut; };
Ladder Two() { return {{ivec3_t(24,24,24), ivec3_t(12,12,12)}, {30.0, 7.5}}; }
//! A fixed smooth-ish field per level (a few plane waves on each level's own raster), or a constant.
std::vector<rvec_t> Field(const Ladder& L, bool constant)
{
    std::vector<rvec_t> V(L.N.size());
    for (size_t l=0;l<L.N.size();l++)
    {
        const ivec3_t N=L.N[l];
        V[l]=rvec_t(size_t(N.x)*N.y*N.z);
        for (int x=0;x<N.x;x++) for (int y=0;y<N.y;y++) for (int z=0;z<N.z;z++)
        {
            const double u=kTwoPi*x/N.x, v=kTwoPi*y/N.y, w=kTwoPi*z/N.z;
            V[l][(size_t(x)*N.y+y)*N.z+z] = constant ? 0.7 : 0.4+0.3*std::cos(u+0.3)-0.2*std::sin(2*v-w)+0.15*std::cos(u+v+w);
        }
    }
    return V;
}
//! Re Tr(D h) over the full Hermitian matrices, from the upper triangles.
double ReTr(const chmat_t& D, const chmat_t& h)
{
    double e=0.0;
    for (size_t i=0;i<D.rows();i++)
    {
        e+=std::real(dcmplx(D(i,i)))*std::real(dcmplx(h(i,i)));
        for (size_t j=i+1;j<D.rows();j++) e+=2.0*std::real(std::conj(dcmplx(D(i,j)))*dcmplx(h(i,j)));
    }
    return e;
}

const rvec3_t kNonTRIM(0.25,0.0,0.125);

void PullbackIsTheFiniteDifference(const Fixture& x)
{
    const Ladder L=Two();
    const GeometryOnlyScreener screen(CollocationEps());
    const std::vector<rvec_t> V=Field(L,false);
    const UnitCell cell=CellOf(x);
    const Block b=Make(x, cell);
    ASSERT_TRUE(b.lc) << x.what << ": the block realises no collocation face";
    const chmat_t D=RandomHermitian(b.n, 0x5DEECE66DULL+b.n);
    const rvec3vec_t G=b.lc->IntegratePotentialPullback(V, D, PhaseOf(kNonTRIM), cell, L.N, L.ecut, screen, -1.0, &cell);
    ASSERT_EQ(G.size(), x.atoms.size());
    auto E=[&](size_t a, const rvec3_t& dr)
    {
        const UnitCell c=CellOf(x, a, dr);
        const Block bb=Make(x, c);   // a FRESH basis: no memo carries across geometries
        return ReTr(D, bb.lc->IntegratePotential(V, PhaseOf(kNonTRIM), c, L.N, L.ecut, screen));
    };
    // RICHARDSON on the central difference (h, h/2): the d shells (alpha 36) and the contracted H core (alpha 51) are narrow enough that
    // the plain O(h^2) truncation would sit at the gate.
    const double h=1e-3;
    auto central=[&](size_t a, int c, double s)
    { rvec3_t e(0,0,0); (c==0 ? e.x : c==1 ? e.y : e.z)=s; return (E(a,e)-E(a,-1.0*e))/(2*s); };
    double big=0.0, worst=0.0;
    for (size_t a=0;a<x.atoms.size();a++)
        for (int c=0;c<3;c++)
        {
            const double fd=(4.0*central(a,c,0.5*h)-central(a,c,h))/3.0;
            const double an= c==0 ? G[a].x : c==1 ? G[a].y : G[a].z;
            std::printf("  [gather pullback] %-22s atom %zu %c  analytic % .10e  FD % .10e  diff %.2e\n", x.what, a, "xyz"[c], an, fd, an-fd);
            big=std::max(big,std::fabs(fd)); worst=std::max(worst,std::fabs(an-fd));
        }
    ASSERT_GT(big, 1e-4) << x.what << ": a vanishing pullback cannot show the walk is wired";
    EXPECT_LT(worst, 1e-6*big) << x.what << ": the gather pullback is not the finite difference of the gather";
}

Fixture Si()   // FCC primitive, a = 10.26, the atoms OFF the origin (no corner atom: the displaced cell must not re-wrap one)
{ return {"Si SIPP_SR cartesian", Matrix3D<double>(0.0,5.13,5.13, 5.13,0.0,5.13, 5.13,5.13,0.0), {{14,{0.1,0.05,0.12}},{14,{0.35,0.3,0.37}}},
          BasisSetData::SIPP_SR, false}; }
Fixture H3()   // CONTRACTED dzvp shells (H 1s: 4 primitives, alpha <= 51) in an orthorhombic box.  Not O: its all-electron 1s (alpha 5223,
               // width 0.01 bohr) is beyond any finite-difference step, and beyond the grid's resolution too.
{ return {"H3 dzvp (contracted)", Matrix3D<double>(9.0,0.0,0.0, 0.0,9.5,0.0, 0.0,0.0,10.0),
          {{1,{0.50,0.50,0.50}},{1,{0.60,0.52,0.50}},{1,{0.47,0.59,0.53}}}, BasisSetData::DZVP, false}; }
Fixture MnO()   // rocksalt, d shells through the spherical view's T
{ return {"MnO VA spherical view", Matrix3D<double>(0.0,4.20,4.20, 4.20,0.0,4.20, 4.20,4.20,0.0), {{25,{0.1,0.1,0.1}},{8,{0.6,0.58,0.62}}},
          BasisSetData::VALENCE_LOWQ_VA, true}; }

} // namespace

TEST(M_CollocationPullback, TheGatherPullbackIsItsFiniteDifference_Uncontracted) {PullbackIsTheFiniteDifference(Si());}
TEST(M_CollocationPullback, TheGatherPullbackIsItsFiniteDifference_Contracted)   {PullbackIsTheFiniteDifference(H3());}
TEST(M_CollocationPullback, TheGatherPullbackIsItsFiniteDifference_Spherical)    {PullbackIsTheFiniteDifference(MnO());}

TEST(M_CollocationPullback, UnderAConstantFieldThePullbacksSumToTheEggBoxOnly)
{
    const Fixture x=Si();
    const Ladder L=Two();
    const GeometryOnlyScreener screen(CollocationEps());
    const UnitCell cell=CellOf(x);
    const Block b=Make(x, cell);
    ASSERT_TRUE(b.lc);
    const chmat_t D=RandomHermitian(b.n, 0x2545F4914F6CDD1DULL);
    const rvec3vec_t G=b.lc->IntegratePotentialPullback(Field(L,true), D, PhaseOf(kNonTRIM), cell, L.N, L.ecut, screen, -1.0, &cell);
    rvec3_t sum(0,0,0); double each=0.0;
    for (const rvec3_t& g : G) { sum=sum+g; each=std::max({each, std::fabs(g.x), std::fabs(g.y), std::fabs(g.z)}); }
    const double s=std::max({std::fabs(sum.x), std::fabs(sum.y), std::fabs(sum.z)});
    std::printf("  [gather pullback, constant V] max|G_A| %.3e   |Σ_A G_A| %.3e\n", each, s);
    ASSERT_GT(each, 1e-4) << "each atom's pullback must be non-trivial (the overlap moves), or this gate tests nothing";
    EXPECT_LT(s, 1e-6*each) << "a rigid translation under a constant field must cost only the raster's egg-box";
}

TEST(M_CollocationPullback, RefusesAMisSizedDensityAndAMissingStructure)
{
    const Fixture x=Si();
    const Ladder L=Two();
    const GeometryOnlyScreener screen(CollocationEps());
    const UnitCell cell=CellOf(x);
    const Block b=Make(x, cell);
    ASSERT_TRUE(b.lc);
    const std::vector<rvec_t> V=Field(L,false);
    EXPECT_THROW((void)b.lc->IntegratePotentialPullback(V, chmat_t(b.n+1), PhaseOf(kNonTRIM), cell, L.N, L.ecut, screen, -1.0, &cell),
                 std::invalid_argument);
    EXPECT_THROW((void)b.lc->IntegratePotentialPullback(V, chmat_t(b.n), PhaseOf(kNonTRIM), cell, L.N, L.ecut, screen, -1.0, nullptr),
                 std::invalid_argument);
}
