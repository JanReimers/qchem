// File: Calculation/tests/WyckoffMap.C  The level-1 solid map (free Wyckoff parameters, the cell fixed) on ANALYTIC MODEL SURFACES only
// (user ruling 2026-11-07, doc/GeometryOptimizationPlan.md §12.3: no real optimization before analytic forces).  The model is a smooth,
// short-ranged periodic pair potential with exact forces, its lattice sum carried until the Gaussians are below roundoff, so what is
// tested is the map and the optimizer alone:
//  - rutile (P4₂/mnm, non-symmorphic, 16 ops): Dimension() = 1 (O 4f's u; Ti 2a is fixed), the optimizer finds the u a direct 1-d
//    minimisation of E(u) finds (an independent route that never touches the map), and every geometry it visits keeps all 16 ops --
//    checked by the map's Asymmetry and, independently, by re-detecting the space group -- even under asymmetric force noise;
//  - wurtzite (P6₃mc, POLAR): two free z's, one of them the floating origin, held by Project; the optimizer finds u = z_O − z_Zn;
//  - the refusals, and a P1 cell's 3N parameters with its three translations held.
#include "gtest/gtest.h"
#include <cmath>
#include <functional>
#include <optional>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>
import qchem.WyckoffMap;
import qchem.Blaze;

using namespace qchem;
using Symmetry::Lattice_3D::SpaceGroup;
using Symmetry::Lattice_3D::AtomSite;

namespace
{
//! A pair potential φ(r) = Σ_k c_k exp(−r²/s_k²) per species pair: smooth, and below 1e-17 of its scale past ~7 widths.
struct Gaussians { std::vector<double> c, s; };
double Phi (const Gaussians& g, double r) {double v=0; for (size_t k=0;k<g.c.size();++k) v+=g.c[k]*std::exp(-r*r/(g.s[k]*g.s[k])); return v;}
double DPhi(const Gaussians& g, double r) {double v=0; for (size_t k=0;k<g.c.size();++k) v+=-2*r/(g.s[k]*g.s[k])*g.c[k]*std::exp(-r*r/(g.s[k]*g.s[k])); return v;}

//! The periodic model surface: E = ½ Σ_{a,b,L}' φ_{Z_a Z_b}(|R_b + A L − R_a|), F_a = Σ_{b,L}' φ'(r) d/r.  \a noise adds a FIXED
//! asymmetric pattern to the forces only (the model of a grid's force noise; not the gradient of E).
class PeriodicPairSurface : public virtual PotentialEnergySurface
{
public:
    using Pair_t=std::function<const Gaussians&(int, int)>;
    PeriodicPairSurface(rmat3d_t A, std::vector<rvec3_t> R0, std::vector<int> Z, Pair_t pair, ivec3_t images, double noise=0.0)
        : itsA(A), itsR0(std::move(R0)), itsZ(std::move(Z)), itsPair(std::move(pair)), itsImages(images), itsNoise(noise) {}
    PESGeometry Reference() const override {return PESGeometry{itsA, itsR0, itsZ};}
    Outcome<PESPoint, PESFailure> ConvergeAt(const std::optional<rmat3d_t>& cell, const std::vector<rvec3_t>& R, const PESRequest& want) override
    {
        EXPECT_TRUE(cell.has_value());
        visited.push_back(R);
        PESPoint p; p.iterations=1; p.energy=Energy(R);
        if (want.forces)
        {
            std::vector<rvec3_t> F(R.size(), rvec3_t(0,0,0));
            Sum(R, [&](size_t a, size_t, const rvec3_t& d, double r, const Gaussians& g){ F[a]+=(DPhi(g, r)/r)*d; });
            for (size_t a=0;a<F.size();++a) F[a]+=itsNoise*rvec3_t(std::sin(1.0+a), std::cos(2.0+3*a), std::sin(0.5+7*a));
            p.forces=F;
        }
        return Outcome<PESPoint, PESFailure>::Ok(p);
    }
    double Energy(const std::vector<rvec3_t>& R) const
    {
        double E=0;
        Sum(R, [&](size_t, size_t, const rvec3_t&, double r, const Gaussians& g){ E+=0.5*Phi(g, r); });
        return E;
    }
    std::vector<std::vector<rvec3_t>> visited;
private:
    template <class F> void Sum(const std::vector<rvec3_t>& R, F f) const
    {
        for (size_t a=0;a<R.size();++a) for (size_t b=0;b<R.size();++b)
        for (int i=-itsImages.x;i<=itsImages.x;++i) for (int j=-itsImages.y;j<=itsImages.y;++j) for (int k=-itsImages.z;k<=itsImages.z;++k)
        {
            if (a==b && i==0 && j==0 && k==0) continue;
            const rvec3_t d=R[b]+itsA*rvec3_t(i,j,k)-R[a];
            f(a, b, d, norm(d), itsPair(itsZ[a], itsZ[b]));
        }
    }
    rmat3d_t             itsA;
    std::vector<rvec3_t> itsR0;
    std::vector<int>     itsZ;
    Pair_t               itsPair;
    ivec3_t              itsImages;
    double               itsNoise;
};

//! Cation-anion attraction with a repulsive core, like-ion repulsion: widths <= 4 bohr, so a 34 bohr image sum is exact to roundoff.
const Gaussians& IonPair(int Za, int Zb)
{
    static const Gaussians cation_anion{{6.0, -1.2}, {1.9, 4.0}}, anion_anion{{1.5}, {3.0}}, cation_cation{{2.0}, {3.2}};
    const bool ca=Za>8, cb=Zb>8;
    return ca!=cb ? cation_anion : (ca ? cation_cation : anion_anion);
}

rmat3d_t Columns(const rvec3_t& a1, const rvec3_t& a2, const rvec3_t& a3)
{
    return rmat3d_t(a1.x, a2.x, a3.x, a1.y, a2.y, a3.y, a1.z, a2.z, a3.z);
}

// RUTILE: a = 8.68, c = 5.59 bohr; Ti 2a (0,0,0), (½,½,½); O 4f (u,u,0), (1−u,1−u,0), (½+u,½−u,½), (½−u,½+u,½).
const rmat3d_t RutileA=Columns(rvec3_t(8.68,0,0), rvec3_t(0,8.68,0), rvec3_t(0,0,5.59));
std::vector<rvec3_t> RutileFrac(double u)
{
    return {rvec3_t(0,0,0), rvec3_t(0.5,0.5,0.5), rvec3_t(u,u,0), rvec3_t(1-u,1-u,0), rvec3_t(0.5+u,0.5-u,0.5), rvec3_t(0.5-u,0.5+u,0.5)};
}
const std::vector<int> RutileZ={22,22,8,8,8,8};
std::vector<rvec3_t> Cartesian(const rmat3d_t& A, const std::vector<rvec3_t>& f) {std::vector<rvec3_t> R; for (const auto& x : f) R.push_back(A*x); return R;}

SpaceGroup Detect(const rmat3d_t& A, const std::vector<rvec3_t>& R, const std::vector<int>& Z)
{
    const rmat3d_t Ainv=Invert(A);
    std::vector<AtomSite> basis;
    for (size_t a=0;a<R.size();++a) basis.push_back({Z[a], Ainv*R[a]});
    return SpaceGroup::Detect(A, basis);
}

//! Golden-section minimum of a 1-d function on [lo, hi] (the independent route: it never touches the map).
double GoldenMin(const std::function<double(double)>& E, double lo, double hi, double tol=1e-10)
{
    const double r=(std::sqrt(5.0)-1)/2;
    double x1=hi-r*(hi-lo), x2=lo+r*(hi-lo), f1=E(x1), f2=E(x2);
    while (hi-lo>tol)
    {
        if (f1<f2) {hi=x2; x2=x1; f2=f1; x1=hi-r*(hi-lo); f1=E(x1);}
        else       {lo=x1; x1=x2; f1=f2; x2=lo+r*(hi-lo); f2=E(x2);}
    }
    return 0.5*(lo+hi);
}

std::string History(const std::vector<OptimizationStep>& h)
{
    std::ostringstream o; o.precision(15);
    for (const auto& s : h) o << "\n  step " << s.index << " E " << s.energy << " max|g| " << s.maxGradient << " |s| " << s.stepLength
                              << " trust " << s.trust << (s.accepted ? "" : " REJECTED");
    return o.str();
}

OptimizerOptions Tight()
{
    // maxStep 1e-4, not tighter: with |g| at 1e-9 the next step's energy change is below the energy's roundoff, and the optimizer
    // rejects a roundoff rise as an energy rise until the trust radius collapses (OpenWork §3 row, 2026-11-09).
    OptimizerOptions o; o.maxGradient=1e-8; o.rmsGradient=1e-8; o.maxStep=1e-4; o.energyChange=1e-14; o.initialCurvature=1.0;
    return o;
}
}

TEST(WyckoffMap, RutileHasOneFreeParameterAndTheOptimizerFindsItsMinimum)
{
    const auto R0=Cartesian(RutileA, RutileFrac(0.28));
    PeriodicPairSurface S(RutileA, R0, RutileZ, IonPair, ivec3_t(4,4,6));
    const SpaceGroup G=Detect(RutileA, R0, RutileZ);
    ASSERT_EQ(G.Order(), 16u) << "P4₂/mnm";
    ASSERT_FALSE(G.isSymmorphic());
    const WyckoffMap map(S.Reference(), G);
    EXPECT_EQ(map.Dimension(), 1u) << "Ti 2a has no free parameter, O 4f has u";
    EXPECT_EQ(map.Orbits(), 2u);
    EXPECT_EQ(map.HeldModes(), 0u) << "P4₂/mnm is not polar";
    EXPECT_LT(map.Asymmetry(map.Start()), 1e-12);

    // The independent route: minimise E(u) directly over the rutile parameterisation, no map involved.
    const double uStar=GoldenMin([&](double u){return S.Energy(Cartesian(RutileA, RutileFrac(u)));}, 0.2, 0.45);
    ASSERT_GT(uStar, 0.21); ASSERT_LT(uStar, 0.44) << "the model's minimum is interior";

    auto out=GeometryOptimizer(S, map, Tight()).Run();
    ASSERT_TRUE(out) << out.Error().details << History(out.Error().history);
    const double u=(Invert(RutileA)*out->geometry.R[2]).x;
    EXPECT_NEAR(u, uStar, 1e-7);
    EXPECT_NEAR((Invert(RutileA)*out->geometry.R[2]).y, u, 1e-12) << "O stays on (u,u,0)";
    EXPECT_NEAR(out->energy, S.Energy(Cartesian(RutileA, RutileFrac(uStar))), 1e-12);
    EXPECT_LT(out->Steps(), 15u);
    ASSERT_GT(S.visited.size(), 2u);
    for (const auto& R : S.visited)
    {
        EXPECT_EQ(Detect(RutileA, R, RutileZ).Order(), 16u) << "every visited geometry keeps the whole group";
        EXPECT_NEAR(R[0].x, 0.0, 1e-13) << "Ti 2a never moves";
    }
}

TEST(WyckoffMap, AsymmetricForceNoiseNeverLeavesTheGroup)
{
    const auto R0=Cartesian(RutileA, RutileFrac(0.28));
    PeriodicPairSurface S(RutileA, R0, RutileZ, IonPair, ivec3_t(4,4,6), /*noise*/1e-3);
    const WyckoffMap map(S.Reference(), Detect(RutileA, R0, RutileZ));
    // The noise carries a symmetric part, which moves the minimum: the claim here is the GROUP, not the optimum.
    OptimizerOptions o=Tight(); o.maxGradient=1e-6; o.rmsGradient=1e-6; o.maxStep=1e-5;
    auto out=GeometryOptimizer(S, map, o).Run();
    ASSERT_TRUE(out) << out.Error().details << History(out.Error().history);
    for (const auto& R : S.visited)
    {
        PESGeometry g{RutileA, R, RutileZ};
        const WyckoffMap at(g, Detect(RutileA, R0, RutileZ));   // THROWS if R is off the group by more than 1e-6 bohr
        EXPECT_LT(at.Asymmetry(at.Start()), 1e-12);
        EXPECT_EQ(Detect(RutileA, R, RutileZ).Order(), 16u);
    }
    for (const auto& h : out->history) EXPECT_GE(h.maxGradient, 0.0);
}

TEST(WyckoffMap, WurtziteIsPolarAndItsFloatingOriginIsHeld)
{
    const double a=6.14, c=9.83, u0=0.36;
    const rmat3d_t A=Columns(rvec3_t(a,0,0), rvec3_t(-a/2, a*std::sqrt(3.0)/2, 0), rvec3_t(0,0,c));
    auto frac=[](double zZn, double u){return std::vector<rvec3_t>{rvec3_t(1.0/3,2.0/3,zZn), rvec3_t(2.0/3,1.0/3,0.5+zZn),
                                                                     rvec3_t(1.0/3,2.0/3,zZn+u), rvec3_t(2.0/3,1.0/3,0.5+zZn+u)};};
    const std::vector<int> Z={30,30,8,8};
    const auto R0=Cartesian(A, frac(0.0, u0));
    PeriodicPairSurface S(A, R0, Z, IonPair, ivec3_t(6,6,4));
    const SpaceGroup G=Detect(A, R0, Z);
    ASSERT_EQ(G.Order(), 12u) << "P6₃mc";
    const WyckoffMap map(S.Reference(), G);
    EXPECT_EQ(map.Dimension(), 2u) << "Zn 2b z and O 2b z";
    EXPECT_EQ(map.HeldModes(), 1u) << "the c axis is polar: the common z shift is a zero mode";

    const double uStar=GoldenMin([&](double u){return S.Energy(Cartesian(A, frac(0.0, u)));}, 0.25, 0.5);
    ASSERT_GT(uStar, 0.26); ASSERT_LT(uStar, 0.49);
    auto out=GeometryOptimizer(S, map, Tight()).Run();
    ASSERT_TRUE(out) << out.Error().details << History(out.Error().history);
    const rmat3d_t Ainv=Invert(A);
    const double zZn=(Ainv*out->geometry.R[0]).z, zO=(Ainv*out->geometry.R[2]).z;
    EXPECT_NEAR(zO-zZn, uStar, 1e-7);
    // The held mode: the centroid along c never moves (the step is J times a vector orthogonal to the common shift).
    double z0=0, z1=0;
    for (size_t k=0;k<4;++k) {z0+=R0[k].z; z1+=out->geometry.R[k].z;}
    EXPECT_NEAR(z1, z0, 1e-10);
    for (const auto& R : S.visited) EXPECT_EQ(Detect(A, R, Z).Order(), 12u);
}

TEST(WyckoffMap, P1HasThreeNParametersAndHoldsTheTranslations)
{
    const rmat3d_t A=Columns(rvec3_t(7.0,0.3,0.1), rvec3_t(0.4,6.5,-0.2), rvec3_t(0.2,0.5,8.1));
    const std::vector<rvec3_t> R={A*rvec3_t(0.1,0.2,0.3), A*rvec3_t(0.55,0.6,0.71)};
    Symmetry::SymOp E{rmat3d_t(1,0,0, 0,1,0, 0,0,1), rvec3_t(0,0,0)};
    const WyckoffMap map(PESGeometry{A, R, {11,9}}, std::vector<Symmetry::SymOp>{E});
    EXPECT_EQ(map.Dimension(), 6u);
    EXPECT_EQ(map.HeldModes(), 3u);
    const rvec_t q=map.Start();
    EXPECT_NEAR(norm(map.At(q).R[1]-R[1]), 0.0, 1e-12) << "Start() is the reference";
}

TEST(WyckoffMap, Refusals)
{
    const auto R0=Cartesian(RutileA, RutileFrac(0.3));
    const SpaceGroup G=Detect(RutileA, R0, RutileZ);
    EXPECT_THROW(WyckoffMap(PESGeometry{std::nullopt, R0, RutileZ}, G), std::invalid_argument) << "no cell: a molecule maps Cartesian";
    auto broken=R0; broken[2].x+=0.01;
    EXPECT_THROW(WyckoffMap(PESGeometry{RutileA, broken, RutileZ}, G), std::invalid_argument) << "the reference is off the group";
    auto swapped=RutileZ; swapped[0]=8; swapped[2]=22;
    EXPECT_THROW(WyckoffMap(PESGeometry{RutileA, R0, swapped}, G), std::invalid_argument) << "an op maps an atom onto another species";
    // Not a group: the 4₂ screw alone, without its powers.
    std::vector<Symmetry::SymOp> open;
    for (const auto& op : G.Ops()) if (op.tau.z>0.25) {open.push_back({op.W, op.tau}); break;}
    EXPECT_THROW(WyckoffMap(PESGeometry{RutileA, R0, RutileZ}, open), std::exception);
    const WyckoffMap map(PESGeometry{RutileA, R0, RutileZ}, G);
    EXPECT_THROW((void)map.At(rvec_t(2, 0.0)), std::invalid_argument);
    EXPECT_THROW((void)map.Gradient(map.Start(), {rvec3_t(0,0,0)}), std::invalid_argument);
}
