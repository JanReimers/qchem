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
#include <algorithm>
#include <cmath>
#include <functional>
#include <optional>
#include <ostream>
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

//! The periodic model surface: E = ½ Σ_{a,b,L}' φ_{Z_a Z_b}(|d|), d = R_b + A L − R_a; F_a = Σ_{b,L}' φ'(r) d/r; and the STRESS, exact:
//! under R → (1+ε)R, A → (1+ε)A every d → (1+ε)d, so ∂E/∂ε = ½ Σ' (φ'(r)/r) d dᵀ (the pair virial) and σ = that / Ω.  The cell is
//! the one ConvergeAt is given.  \a noise adds a FIXED asymmetric pattern to the forces only (the model of a grid's force noise).
class PeriodicPairSurface : public virtual PotentialEnergySurface
{
public:
    using Pair_t=std::function<const Gaussians&(int, int)>;
    PeriodicPairSurface(rmat3d_t A, std::vector<rvec3_t> R0, std::vector<int> Z, Pair_t pair, ivec3_t images, double noise=0.0)
        : itsA(A), itsR0(std::move(R0)), itsZ(std::move(Z)), itsPair(std::move(pair)), itsImages(images), itsNoise(noise) {}
    //! The images sized PER CELL: every lattice vector within \a cutoff (bohr) of the home cell, for a surface whose cell moves.
    PeriodicPairSurface(rmat3d_t A, std::vector<rvec3_t> R0, std::vector<int> Z, Pair_t pair, double cutoff)
        : itsA(A), itsR0(std::move(R0)), itsZ(std::move(Z)), itsPair(std::move(pair)), itsImages(0,0,0), itsNoise(0.0), itsCutoff(cutoff) {}
    PESGeometry Reference() const override {return PESGeometry{itsA, itsR0, itsZ};}
    Outcome<PESPoint, PESFailure> ConvergeAt(const std::optional<rmat3d_t>& cell, const std::vector<rvec3_t>& R, const PESRequest& want) override
    {
        EXPECT_TRUE(cell.has_value());
        const rmat3d_t A=cell ? *cell : itsA;
        visited.push_back(R);
        visitedCells.push_back(A);
        PESPoint p; p.iterations=1; p.energy=Energy(A, R);
        if (want.forces)
        {
            std::vector<rvec3_t> F(R.size(), rvec3_t(0,0,0));
            Sum(A, R, [&](size_t a, size_t, const rvec3_t& d, double r, const Gaussians& g){ F[a]+=(DPhi(g, r)/r)*d; });
            for (size_t a=0;a<F.size();++a) F[a]+=itsNoise*rvec3_t(std::sin(1.0+a), std::cos(2.0+3*a), std::sin(0.5+7*a));
            p.forces=F;
        }
        if (want.stress) p.stress=Stress(A, R);
        return Outcome<PESPoint, PESFailure>::Ok(p);
    }
    double Energy(const std::vector<rvec3_t>& R) const {return Energy(itsA, R);}
    double Energy(const rmat3d_t& A, const std::vector<rvec3_t>& R) const
    {
        double E=0;
        Sum(A, R, [&](size_t, size_t, const rvec3_t&, double r, const Gaussians& g){ E+=0.5*Phi(g, r); });
        return E;
    }
    rmat3d_t Stress(const rmat3d_t& A, const std::vector<rvec3_t>& R) const
    {
        rmat3d_t W(0,0,0, 0,0,0, 0,0,0);
        Sum(A, R, [&](size_t, size_t, const rvec3_t& d, double r, const Gaussians& g)
        {
            const double c=0.5*DPhi(g, r)/r;
            for (size_t i=1;i<=3;++i) for (size_t j=1;j<=3;++j) W(i,j)+=c*Comp(d,i-1)*Comp(d,j-1);
        });
        return W*(1.0/std::fabs(Determinant(A)));
    }
    std::vector<std::vector<rvec3_t>> visited;
    std::vector<rmat3d_t>             visitedCells;
private:
    static double Comp(const rvec3_t& v, size_t k) {return k==0 ? v.x : k==1 ? v.y : v.z;}
    template <class F> void Sum(const rmat3d_t& A, const std::vector<rvec3_t>& R, F f) const
    {
        ivec3_t n=itsImages;
        if (itsCutoff>0.0)   // the cell's three heights Ω / |a_j x a_k|: n_i images reach the cutoff along a_i's direction
        {
            const double Omega=std::fabs(Determinant(A));
            int m[3];
            for (size_t i=0;i<3;++i) m[i]=int(std::ceil(itsCutoff/(Omega/norm(Cross(A.GetCol((i+1)%3+1), A.GetCol((i+2)%3+1)))))) + 1;
            n=ivec3_t(m[0], m[1], m[2]);
        }
        for (size_t a=0;a<R.size();++a) for (size_t b=0;b<R.size();++b)
        for (int i=-n.x;i<=n.x;++i) for (int j=-n.y;j<=n.y;++j) for (int k=-n.z;k<=n.z;++k)
        {
            if (a==b && i==0 && j==0 && k==0) continue;
            const rvec3_t d=R[b]+A*rvec3_t(i,j,k)-R[a];
            f(a, b, d, norm(d), itsPair(itsZ[a], itsZ[b]));
        }
    }
    rmat3d_t             itsA;
    std::vector<rvec3_t> itsR0;
    std::vector<int>     itsZ;
    Pair_t               itsPair;
    ivec3_t              itsImages;
    double               itsNoise;
    double               itsCutoff=0.0;
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

std::string History(const std::vector<Optimize::Iteration>& h)
{
    std::ostringstream o; o.precision(15);
    for (const auto& s : h) o << "\n  step " << s.index << " E " << s.value << " max|g| " << s.maxGradient << " |s| " << s.stepLength
                              << " trust " << s.bound << (s.accepted ? "" : " REJECTED");
    return o.str();
}

//! The test's own name for the stepper: a local enum, so gtest finds \c PrintTo by ADL here and ctest names the case by it, not by
//! its bytes.
enum class Algo {BFGS_RFO, L_BFGS};
std::string Name(Algo a) {return a==Algo::BFGS_RFO ? "BFGS_RFO" : "L_BFGS";}
void PrintTo(Algo a, std::ostream* os) {*os << Name(a);}
OptimizerOptions::Algorithm Of(Algo a) {return a==Algo::BFGS_RFO ? OptimizerOptions::Algorithm::BFGS_RFO : OptimizerOptions::Algorithm::L_BFGS;}

OptimizerOptions Tight()
{
    // maxStep 1e-6.  This test once needed 1e-4: with |g| at 1e-9 the next step's energy change was below the energy's roundoff, and
    // the optimizer rejected a roundoff rise as a rise until its radius collapsed (lane H, 2026-11-09).  The noise rules of
    // V-OPTIMIZE-DIP §4 (qcOptimize's BFGS-RFO) fixed it; this test is their gate on a real map.
    OptimizerOptions o; o.maxGradient=1e-8; o.rmsGradient=1e-8; o.maxStep=1e-6; o.energyChange=1e-14; o.initialCurvature=1.0;
    return o;
}
}

//! THE LSP PAIR (V-OPTIMIZE-DIP §8 step 3): the solid map's optimizer gates run, unchanged, under each stepper.
class WyckoffMapBoth : public ::testing::TestWithParam<Algo>
{
protected:
    OptimizerOptions With(OptimizerOptions o) const {o.algorithm=Of(GetParam()); return o;}
};

TEST_P(WyckoffMapBoth, RutileHasOneFreeParameterAndTheOptimizerFindsItsMinimum)
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

    auto out=GeometryOptimizer(S, map, With(Tight())).Run();
    ASSERT_TRUE(out) << out.Error().details << History(out.Error().history);
    const double u=(Invert(RutileA)*out->geometry.R[2]).x;
    EXPECT_NEAR(u, uStar, 1e-7);
    EXPECT_NEAR((Invert(RutileA)*out->geometry.R[2]).y, u, 1e-12) << "O stays on (u,u,0)";
    EXPECT_NEAR(out->value, S.Energy(Cartesian(RutileA, RutileFrac(uStar))), 1e-12);
    EXPECT_LT(out->Steps(), 15u);
    ASSERT_GT(S.visited.size(), 2u);
    for (const auto& R : S.visited)
    {
        EXPECT_EQ(Detect(RutileA, R, RutileZ).Order(), 16u) << "every visited geometry keeps the whole group";
        EXPECT_NEAR(R[0].x, 0.0, 1e-13) << "Ti 2a never moves";
    }
}

TEST_P(WyckoffMapBoth, AsymmetricForceNoiseNeverLeavesTheGroup)
{
    const auto R0=Cartesian(RutileA, RutileFrac(0.28));
    PeriodicPairSurface S(RutileA, R0, RutileZ, IonPair, ivec3_t(4,4,6), /*noise*/1e-3);
    const WyckoffMap map(S.Reference(), Detect(RutileA, R0, RutileZ));
    // The noise carries a symmetric part, which moves the minimum: the claim here is the GROUP, not the optimum.
    OptimizerOptions o=Tight(); o.maxGradient=1e-6; o.rmsGradient=1e-6; o.maxStep=1e-5;
    auto out=GeometryOptimizer(S, map, With(o)).Run();
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

TEST_P(WyckoffMapBoth, WurtziteIsPolarAndItsFloatingOriginIsHeld)
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
    auto out=GeometryOptimizer(S, map, With(Tight())).Run();
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
    EXPECT_THROW((void)map.Gradient(map.Start(), PESPoint{.energy=0.0, .forces=std::vector<rvec3_t>{rvec3_t(0,0,0)}}), std::invalid_argument);
    EXPECT_THROW((void)map.Gradient(map.Start(), PESPoint{}), std::invalid_argument) << "no forces";
}

//====================================================================================================================================
// THE CELL HALF, step 1 (doc/Records/CellHalf_Design.md §4 step 1, ruled 2026-11-10): CellWyckoffMap on the same model surfaces, now
// with their exact stress.  The claims: the model's stress IS the strain derivative (FD, every component, and nothing under a
// rotation); the free strains ARE the group's (rutile 2, wurtzite 2, rock salt 1, P1 6); the map's gradient IS the FD of E through
// the map (a triclinic strain, where F and σ do not commute); and the optimizer finds the minimum a direct minimisation of E(a, c, u)
// finds, under both steppers, every visited cell keeping the whole group.
namespace
{
//! Nelder-Mead on a few variables: the independent route for the cell (it never touches a map, a stress or a force).
std::vector<double> NelderMead(const std::function<double(const std::vector<double>&)>& E, std::vector<double> x0, const std::vector<double>& step,
                               size_t maxEval=20000)
{
    const size_t n=x0.size();
    std::vector<std::vector<double>> P(n+1, x0);
    for (size_t i=0;i<n;++i) P[i+1][i]+=step[i];
    std::vector<double> f(n+1);
    for (size_t i=0;i<=n;++i) f[i]=E(P[i]);
    for (size_t eval=n+1; eval<maxEval;)
    {
        std::vector<size_t> o(n+1); for (size_t i=0;i<=n;++i) o[i]=i;
        std::sort(o.begin(), o.end(), [&](size_t a, size_t b){return f[a]<f[b];});
        const size_t lo=o[0], hi=o[n], nh=o[n-1];
        double spread=0; for (size_t i=0;i<=n;++i) for (size_t k=0;k<n;++k) spread=std::max(spread, std::fabs(P[i][k]-P[lo][k]));
        if (spread<1e-9) break;
        std::vector<double> c(n, 0.0);
        for (size_t i=0;i<=n;++i) if (i!=hi) for (size_t k=0;k<n;++k) c[k]+=P[i][k]/double(n);
        auto along=[&](double t){std::vector<double> x(n); for (size_t k=0;k<n;++k) x[k]=c[k]+t*(P[hi][k]-c[k]); return x;};
        auto xr=along(-1.0); double fr=E(xr); ++eval;
        if (fr<f[lo]) { auto xe=along(-2.0); double fe=E(xe); ++eval; if (fe<fr) {P[hi]=xe; f[hi]=fe;} else {P[hi]=xr; f[hi]=fr;} }
        else if (fr<f[nh]) {P[hi]=xr; f[hi]=fr;}
        else
        {
            auto xc=along(fr<f[hi] ? -0.5 : 0.5); double fc=E(xc); ++eval;
            if (fc<std::min(fr, f[hi])) {P[hi]=xc; f[hi]=fc;}
            else for (size_t i=0;i<=n;++i) if (i!=lo) { for (size_t k=0;k<n;++k) P[i][k]=P[lo][k]+0.5*(P[i][k]-P[lo][k]); f[i]=E(P[i]); ++eval; }
        }
    }
    size_t best=0; for (size_t i=1;i<=n;++i) if (f[i]<f[best]) best=i;
    return P[best];
}

//! The CELL tests' pair set: IonPair's cation-anion term with stiffer like-ion repulsion.  IonPair is no rutile once the cell may move
//! (O collapses to u = ½, a = 5.5, c = 7.7); this one keeps rutile a STABLE minimum over all 24 P1 coordinates (a = 10.32, c = 7.16,
//! u = 0.301; the numerical P1 Hessian's eigenvalues are the 3 translations' zeros and the rest >= 0.022 Ha/bohr²; scipy scan,
//! lane H 2026-11-10).  Widths <= 4 bohr: the per-cell 32 bohr image cutoff holds it to roundoff.
const Gaussians& CellIonPair(int Za, int Zb)
{
    static const Gaussians cation_anion{{6.0, -1.2}, {1.9, 4.0}}, anion_anion{{3.0}, {3.0}}, cation_cation{{4.0}, {3.8}};
    const bool ca=Za>8, cb=Zb>8;
    return ca!=cb ? cation_anion : (ca ? cation_cation : anion_anion);
}
constexpr double kCutoff=32.0;
const rmat3d_t CellRutileA=Columns(rvec3_t(9.8,0,0), rvec3_t(0,9.8,0), rvec3_t(0,0,6.8));   //!< ~5% off the model's own cell

rmat3d_t RutileCell(double a, double c) {return Columns(rvec3_t(a,0,0), rvec3_t(0,a,0), rvec3_t(0,0,c));}
rmat3d_t WurtziteCell(double a, double c) {return Columns(rvec3_t(a,0,0), rvec3_t(-a/2, a*std::sqrt(3.0)/2, 0), rvec3_t(0,0,c));}
std::vector<rvec3_t> WurtziteFrac(double zZn, double u)
{
    return {rvec3_t(1.0/3,2.0/3,zZn), rvec3_t(2.0/3,1.0/3,0.5+zZn), rvec3_t(1.0/3,2.0/3,zZn+u), rvec3_t(2.0/3,1.0/3,0.5+zZn+u)};
}
const std::vector<int> WurtziteZ={30,30,8,8};
double Length(const rmat3d_t& A, size_t col) {return norm(A.GetCol(col));}
std::vector<Symmetry::SymOp> P1() {return {Symmetry::SymOp{rmat3d_t(1,0,0, 0,1,0, 0,0,1), rvec3_t(0,0,0)}};}
double Trace3(const rmat3d_t& a) {return a(1,1)+a(2,2)+a(3,3);}

//! The rutile model's minimum over (a, c, u), found directly (cached: three tests share it).
const std::vector<double>& RutileDirectMinimum()
{
    static const std::vector<double> x=[]
    {
        PeriodicPairSurface S(CellRutileA, Cartesian(CellRutileA, RutileFrac(0.28)), RutileZ, CellIonPair, kCutoff);
        auto E=[&](const std::vector<double>& p){const rmat3d_t A=RutileCell(p[0], p[1]); return S.Energy(A, Cartesian(A, RutileFrac(p[2])));};
        auto x=NelderMead(E, {9.8, 6.8, 0.28}, {0.2, 0.2, 0.02});
        return NelderMead(E, x, {0.01, 0.01, 0.001});   // a restart, so a collapsed simplex cannot stall it
    }();
    return x;
}
}

TEST(CellWyckoffMap, TheModelStressIsTheStrainDerivativeAndARotationCostsNothing)
{
    // A triclinic cell, so every component of σ is independent.
    const rmat3d_t A=Columns(rvec3_t(7.0,0.3,0.1), rvec3_t(0.4,6.5,-0.2), rvec3_t(0.2,0.5,8.1));
    const std::vector<rvec3_t> R={A*rvec3_t(0.1,0.2,0.3), A*rvec3_t(0.55,0.6,0.71)};
    PeriodicPairSurface S(A, R, {11,9}, IonPair, ivec3_t(5,5,5));
    const rmat3d_t sigma=S.Stress(A, R);
    const double Omega=std::fabs(Determinant(A)), h=1e-5;
    auto Estrained=[&](const rmat3d_t& D, double t)
    {
        const rmat3d_t F=rmat3d_t(1,0,0, 0,1,0, 0,0,1)+t*D;
        std::vector<rvec3_t> RF; for (const auto& r : R) RF.push_back(F*r);
        return S.Energy(F*A, RF);
    };
    double scale=0; for (size_t i=1;i<=3;++i) for (size_t j=1;j<=3;++j) scale=std::max(scale, std::fabs(Omega*sigma(i,j)));
    ASSERT_GT(scale, 1e-3) << "a model far from its minimum has a sizable stress";
    for (size_t i=1;i<=3;++i) for (size_t j=1;j<=3;++j)
    {
        EXPECT_NEAR(sigma(i,j), sigma(j,i), 1e-15*scale) << "σ is symmetric";
        rmat3d_t D(0,0,0, 0,0,0, 0,0,0); D(i,j)=1.0;
        const double fd=(Estrained(D, h)-Estrained(D, -h))/(2*h);
        EXPECT_NEAR(fd, Omega*sigma(i,j), 1e-7*scale) << "∂E/∂ε_" << i << j;
    }
    for (size_t i=1;i<=3;++i) for (size_t j=i+1;j<=3;++j)
    {
        rmat3d_t W(0,0,0, 0,0,0, 0,0,0); W(i,j)=1.0; W(j,i)=-1.0;
        EXPECT_NEAR((Estrained(W, h)-Estrained(W, -h))/(2*h), 0.0, 1e-9*scale) << "an antisymmetric ε is a rotation: no energy, no stress";
    }
}

TEST(CellWyckoffMap, TheFreeStrainsAreTheGroups)
{
    struct Case {std::string name; rmat3d_t A; std::vector<rvec3_t> R; std::vector<int> Z; size_t cellDim, intDim;};
    const rmat3d_t W=WurtziteCell(6.14, 9.83);
    const double a=10.6;
    const rmat3d_t fcc=Columns(rvec3_t(0,a/2,a/2), rvec3_t(a/2,0,a/2), rvec3_t(a/2,a/2,0));
    const std::vector<Case> cases={
        {"rutile",   RutileA, Cartesian(RutileA, RutileFrac(0.3)), RutileZ, 2, 1},
        {"wurtzite", W, Cartesian(W, WurtziteFrac(0.0, 0.37)), WurtziteZ, 2, 2},
        {"rocksalt", fcc, {rvec3_t(0,0,0), fcc*rvec3_t(0.5,0.5,0.5)}, {11,17}, 1, 0}};
    for (const auto& c : cases)
    {
        const SpaceGroup G=Detect(c.A, c.R, c.Z);
        const CellWyckoffMap map(PESGeometry{c.A, c.R, c.Z}, G);
        EXPECT_EQ(map.CellDimension(), c.cellDim) << c.name;
        EXPECT_EQ(map.InternalDimension(), c.intDim) << c.name;
        const rmat3d_t Ainv=Invert(c.A);
        const auto& E=map.StrainBasis();
        for (size_t k=0;k<E.size();++k)
        {
            for (size_t l=0;l<E.size();++l)
            {
                double f=0; for (size_t i=1;i<=3;++i) for (size_t j=1;j<=3;++j) f+=E[k](i,j)*E[l](i,j);
                EXPECT_NEAR(f, k==l ? 1.0 : 0.0, 1e-12) << c.name << ": Frobenius-orthonormal";
            }
            for (size_t i=1;i<=3;++i) for (size_t j=1;j<=3;++j) EXPECT_NEAR(E[k](i,j), E[k](j,i), 1e-14) << c.name << ": symmetric";
            for (const auto& op : G.Ops())
            {
                const rmat3d_t C=c.A*op.W*Ainv, D=C*E[k]*Transpose(C);
                for (size_t i=1;i<=3;++i) for (size_t j=1;j<=3;++j) EXPECT_NEAR(D(i,j), E[k](i,j), 1e-12) << c.name << ": invariant";
            }
        }
        // Every free strain keeps every op: the strained crystal re-detects the whole group.
        rvec_t q=map.Start(); for (size_t k=0;k<map.CellDimension();++k) q[map.InternalDimension()+k]+=0.3*(k+1);
        const PESGeometry g=map.At(q);
        EXPECT_EQ(Detect(*g.cell, g.R, g.Z).Order(), G.Order()) << c.name;
    }
    const CellWyckoffMap p1(PESGeometry{RutileA, Cartesian(RutileA, RutileFrac(0.3)), RutileZ}, P1());
    EXPECT_EQ(p1.CellDimension(), 6u) << "triclinic: all six strains";
    EXPECT_EQ(p1.InternalDimension(), 18u);
    EXPECT_EQ(p1.HeldModes(), 3u);
    EXPECT_EQ(p1.PointRequest().stress, true);
    EXPECT_NEAR(Length(*p1.At(p1.Start()).cell, 1), Length(RutileA, 1), 1e-14) << "Start() is the reference cell";
}

TEST(CellWyckoffMap, TheGradientIsTheFiniteDifferenceOfTheEnergyThroughTheMap)
{
    const auto R0=Cartesian(RutileA, RutileFrac(0.3));
    PeriodicPairSurface S(RutileA, R0, RutileZ, IonPair, kCutoff);
    for (const bool triclinic : {true, false})
    {
        const CellWyckoffMap map=triclinic ? CellWyckoffMap(S.Reference(), P1()) : CellWyckoffMap(S.Reference(), Detect(RutileA, R0, RutileZ));
        // Away from the reference in every coordinate: under a general (triclinic) strain F and σ do not commute.
        rvec_t q=map.Start();
        for (size_t i=0;i<q.size();++i) q[i]+=(i<map.InternalDimension() ? 0.05 : 0.25)*std::sin(1.0+2.0*i);
        const PESGeometry g=map.At(q);
        auto p=S.ConvergeAt(g.cell, g.R, map.PointRequest());
        ASSERT_TRUE(p);
        const rvec_t grad=map.Gradient(q, *p);
        const double h=1e-5;
        for (size_t i=0;i<q.size();++i)
        {
            rvec_t qp=q, qm=q; qp[i]+=h; qm[i]-=h;
            const PESGeometry gp=map.At(qp), gm=map.At(qm);
            const double fd=(S.Energy(*gp.cell, gp.R)-S.Energy(*gm.cell, gm.R))/(2*h);
            EXPECT_NEAR(grad[i], fd, 1e-8) << (triclinic ? "P1" : "P4₂/mnm") << " coordinate " << i
                                           << (i<map.InternalDimension() ? " (internal)" : " (cell)");
        }
        // The model Hessian is symmetric and positive on the cell block.
        const size_t n=map.Dimension();
        rmat_t H(n, n);
        for (size_t j=0;j<n;++j) { rvec_t e(n, 0.0); e[j]=1.0; const rvec_t c=map.ApplyModelHessian(q, e, 0.5); for (size_t i=0;i<n;++i) H(i,j)=c[i]; }
        for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) EXPECT_NEAR(H(i,j), H(j,i), 1e-15);
        for (size_t k=map.InternalDimension();k<n;++k) EXPECT_GT(H(k,k), 0.0);
    }
    EXPECT_THROW((void)CellWyckoffMap(S.Reference(), P1()).Gradient(rvec_t(24, 0.0), PESPoint{.energy=0.0, .forces=std::vector<rvec3_t>(6, rvec3_t(0,0,0))}),
                 std::invalid_argument) << "a cell map needs the stress";
}

//! Tight() with the model's energy noise DECLARED: a per-cell lattice sum of ~10⁵ O(1) Gaussians, mostly cancelling, carries ~1e-13 Ha
//! of roundoff, above the 16 ε|E| ≈ 1e-14 the optimizer assumes by default (rutile, L-BFGS: a 2.2e-14 rise at max|g| 1e-9 read as a
//! real one).  Declaring it is what the cost's ValueNoise is for (V-OPTIMIZE-DIP §4 rule 1); the gradient thresholds are unchanged.
OptimizerOptions CellTight() {OptimizerOptions o=Tight(); o.energyNoise=1e-13; return o;}

class CellWyckoffMapBoth : public ::testing::TestWithParam<Algo>
{
protected:
    OptimizerOptions With(OptimizerOptions o) const {o.algorithm=Of(GetParam()); return o;}
};

TEST_P(CellWyckoffMapBoth, RutileRelaxesACAndUToTheDirectMinimum)
{
    const auto R0=Cartesian(CellRutileA, RutileFrac(0.28));
    PeriodicPairSurface S(CellRutileA, R0, RutileZ, CellIonPair, kCutoff);
    const SpaceGroup G=Detect(CellRutileA, R0, RutileZ);
    const CellWyckoffMap map(S.Reference(), G);
    ASSERT_EQ(map.Dimension(), 3u) << "u, and the strains a and c";

    const auto& x=RutileDirectMinimum();
    ASSERT_GT(std::fabs(x[0]-9.8)+std::fabs(x[1]-6.8), 0.3) << "the model's cell is not the reference: the cell has work to do";
    const rmat3d_t Astar=RutileCell(x[0], x[1]);
    const double Estar=S.Energy(Astar, Cartesian(Astar, RutileFrac(x[2])));

    auto out=GeometryOptimizer(S, map, With(CellTight())).Run();
    ASSERT_TRUE(out) << out.Error().details << History(out.Error().history);
    const rmat3d_t A=*out->geometry.cell;
    EXPECT_NEAR(Length(A, 1), x[0], 1e-5);
    EXPECT_NEAR(Length(A, 2), x[0], 1e-5);
    EXPECT_NEAR(Length(A, 3), x[1], 1e-5);
    EXPECT_NEAR(A(1,2), 0.0, 1e-12) << "the cell stays tetragonal";
    EXPECT_NEAR((Invert(A)*out->geometry.R[2]).x, x[2], 1e-6);
    EXPECT_LT(out->value, Estar+1e-12) << "at least as low as the direct minimum";
    EXPECT_NEAR(out->value, Estar, 1e-10);
    const rmat3d_t sigma=S.Stress(A, out->geometry.R);
    for (size_t i=1;i<=3;++i) for (size_t j=1;j<=3;++j) EXPECT_LT(std::fabs(sigma(i,j)), 1e-9) << "σ_" << i << j << " at the minimum";
    EXPECT_LT(out->Steps(), 40u);
    ASSERT_EQ(S.visited.size(), S.visitedCells.size());
    for (size_t v=0;v<S.visited.size();++v) EXPECT_EQ(Detect(S.visitedCells[v], S.visited[v], RutileZ).Order(), 16u) << "visited cell " << v;
}

TEST_P(CellWyckoffMapBoth, WurtziteRelaxesACAndUWithItsOriginHeld)
{
    const rmat3d_t A0=WurtziteCell(6.14, 9.83);
    const auto R0=Cartesian(A0, WurtziteFrac(0.0, 0.36));
    PeriodicPairSurface S(A0, R0, WurtziteZ, CellIonPair, kCutoff);
    const CellWyckoffMap map(S.Reference(), Detect(A0, R0, WurtziteZ));
    ASSERT_EQ(map.Dimension(), 4u);
    ASSERT_EQ(map.HeldModes(), 1u);
    auto E=[&](const std::vector<double>& p){const rmat3d_t A=WurtziteCell(p[0], p[1]); return S.Energy(A, Cartesian(A, WurtziteFrac(0.0, p[2])));};
    auto x=NelderMead(E, {6.14, 9.83, 0.36}, {0.2, 0.2, 0.02});
    x=NelderMead(E, x, {0.01, 0.01, 0.001});

    auto out=GeometryOptimizer(S, map, With(CellTight())).Run();
    ASSERT_TRUE(out) << out.Error().details << History(out.Error().history);
    const rmat3d_t A=*out->geometry.cell, Ainv=Invert(A);
    EXPECT_NEAR(Length(A, 1), x[0], 1e-5);
    EXPECT_NEAR(Length(A, 3), x[1], 1e-5);
    EXPECT_NEAR((Ainv*out->geometry.R[2]).z-(Ainv*out->geometry.R[0]).z, x[2], 1e-6);
    EXPECT_NEAR(out->value, E(x), 1e-10);
    // The held mode: the FRACTIONAL centroid along c never moves.
    double z0=0, z1=0;
    for (size_t k=0;k<4;++k) {z0+=(Invert(A0)*R0[k]).z; z1+=(Ainv*out->geometry.R[k]).z;}
    EXPECT_NEAR(z1, z0, 1e-10);
    for (size_t v=0;v<S.visited.size();++v) EXPECT_EQ(Detect(S.visitedCells[v], S.visited[v], WurtziteZ).Order(), 12u);
}

TEST_P(CellWyckoffMapBoth, ATriclinicStartFindsTheRutileMinimumThroughAllSixStrains)
{
    // Start OFF the group (a sheared, stretched cell, displaced atoms), mapped in P1: 18 positions and 6 strains, 3 translations held.
    const auto& x=RutileDirectMinimum();
    const rmat3d_t Astar=RutileCell(x[0], x[1]);
    const auto Rstar=Cartesian(Astar, RutileFrac(x[2]));
    const PeriodicPairSurface ref(Astar, Rstar, RutileZ, CellIonPair, kCutoff);
    const double Estar=ref.Energy(Astar, Rstar);
    const rmat3d_t F(1.02, 0.01, -0.008,  0.01, 0.985, 0.012,  -0.008, 0.012, 1.015);
    std::vector<rvec3_t> R;
    for (size_t a=0;a<Rstar.size();++a) R.push_back(F*Rstar[a]+0.05*rvec3_t(std::sin(1.0+a), std::cos(2.0+a), std::sin(3.0+2*a)));
    PeriodicPairSurface S(F*Astar, R, RutileZ, CellIonPair, kCutoff);
    const CellWyckoffMap map(S.Reference(), P1());
    ASSERT_EQ(map.Dimension(), 24u);

    OptimizerOptions o=CellTight(); o.maxSteps=300;
    auto out=GeometryOptimizer(S, map, With(o)).Run();
    ASSERT_TRUE(out) << out.Error().details << History(out.Error().history);
    EXPECT_NEAR(out->value, Estar, 1e-9) << "the rutile minimum, reached without its group";
    const rmat3d_t A=*out->geometry.cell;
    EXPECT_NEAR(Length(A, 1), x[0], 1e-4);
    EXPECT_NEAR(Length(A, 2), x[0], 1e-4);
    EXPECT_NEAR(Length(A, 3), x[1], 1e-4);
    const rmat3d_t sigma=S.Stress(A, out->geometry.R);
    for (size_t i=1;i<=3;++i) for (size_t j=1;j<=3;++j) EXPECT_LT(std::fabs(sigma(i,j)), 1e-8) << "σ_" << i << j;
    EXPECT_EQ(Detect(A, out->geometry.R, RutileZ).Order(), 16u) << "the group comes back by itself";
}

INSTANTIATE_TEST_SUITE_P(CellWyckoffMap, CellWyckoffMapBoth, ::testing::Values(Algo::BFGS_RFO, Algo::L_BFGS),
                         [](const ::testing::TestParamInfo<Algo>& i){return Name(i.param);});

INSTANTIATE_TEST_SUITE_P(WyckoffMap, WyckoffMapBoth, ::testing::Values(Algo::BFGS_RFO, Algo::L_BFGS),
                         [](const ::testing::TestParamInfo<Algo>& i){return Name(i.param);});
