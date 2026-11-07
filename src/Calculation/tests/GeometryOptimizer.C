// File: Calculation/tests/GeometryOptimizer.C  The optimizer on ANALYTIC MODEL SURFACES only (user ruling 2026-11-07,
// doc/GeometryOptimizationPlan.md §12.3: no real optimization before analytic forces).  Each model is a PotentialEnergySurface with
// exact forces, so what is tested is the optimizer and the coordinate map alone: step counts on a quadratic, the trust radius on
// Rosenbrock, conditioning on a soft mode, a symmetric step that stays on its group under asymmetric force noise, rigid-mode
// projection, and the refusals.
#include "gtest/gtest.h"
#include <cmath>
#include <functional>
#include <optional>
#include <stdexcept>
#include <vector>
import qchem.GeometryOptimizer;
import qchem.Blaze;

using namespace qchem;

namespace
{
//! A model surface over Cartesian positions: E(R) and F(R) = −∇E from closures; counts its evaluations, records every geometry, and can
//! fail at a chosen evaluation.
class ModelSurface : public virtual PotentialEnergySurface
{
public:
    using E_t = std::function<double(const std::vector<rvec3_t>&)>;
    using F_t = std::function<std::vector<rvec3_t>(const std::vector<rvec3_t>&)>;
    ModelSurface(std::vector<rvec3_t> R0, std::vector<int> Z, E_t e, F_t f) : itsR0(std::move(R0)), itsZ(std::move(Z)), itsE(std::move(e)), itsF(std::move(f)) {}
    PESGeometry Reference() const override {return PESGeometry{std::nullopt, itsR0, itsZ};}
    Outcome<PESPoint, PESFailure> ConvergeAt(const std::optional<rmat3d_t>&, const std::vector<rvec3_t>& R, const PESRequest& want) override
    {
        ++evaluations; visited.push_back(R);
        if (failAt && evaluations==*failAt) return Outcome<PESPoint, PESFailure>::Fail({PESFailure::Why::SCF, "model SCF failure", std::nullopt});
        PESPoint p; p.energy=itsE(R); p.iterations=1;
        if (want.forces) p.forces=itsF(R);
        return Outcome<PESPoint, PESFailure>::Ok(p);
    }
    size_t evaluations=0;
    std::optional<size_t> failAt;
    std::vector<std::vector<rvec3_t>> visited;
private:
    std::vector<rvec3_t> itsR0;
    std::vector<int>     itsZ;
    E_t itsE; F_t itsF;
};

double C(const rvec3_t& v, size_t k) {return k==0 ? v.x : k==1 ? v.y : v.z;}
double& C(rvec3_t& v, size_t k) {return k==0 ? v.x : k==1 ? v.y : v.z;}
std::vector<double> Flat(const std::vector<rvec3_t>& R) {std::vector<double> x; for (const auto& r : R) for (size_t k=0;k<3;++k) x.push_back(C(r,k)); return x;}
std::vector<rvec3_t> Unflat(const std::vector<double>& x) {std::vector<rvec3_t> R; for (size_t a=0;a<x.size()/3;++a) R.emplace_back(x[3*a],x[3*a+1],x[3*a+2]); return R;}

//! E = ½ (x − x0)ᵀ A (x − x0) over the flattened positions.
ModelSurface Bowl(const std::vector<rvec3_t>& R0, const std::vector<double>& x0, const rmat_t& A)
{
    auto dx=[x0](const std::vector<rvec3_t>& R){auto x=Flat(R); for (size_t i=0;i<x.size();++i) x[i]-=x0[i]; return x;};
    auto e=[dx,A](const std::vector<rvec3_t>& R){auto d=dx(R); double E=0; for (size_t i=0;i<d.size();++i) for (size_t j=0;j<d.size();++j) E+=0.5*d[i]*A(i,j)*d[j]; return E;};
    auto f=[dx,A](const std::vector<rvec3_t>& R){auto d=dx(R); std::vector<double> F(d.size(),0.0); for (size_t i=0;i<d.size();++i) for (size_t j=0;j<d.size();++j) F[i]-=A(i,j)*d[j]; return Unflat(F);};
    return ModelSurface(R0, std::vector<int>(R0.size(), 1), e, f);
}

//! A bent triatomic X–O–X in the xz plane: two bond springs, an angle spring and an X–X repulsion -- internal coordinates only, so the
//! surface is translation and rotation invariant and C2v-symmetric.  \a noise adds a FIXED asymmetric vector to the forces only (the
//! model of a grid's force noise: not the gradient of E).
ModelSurface Triatomic(const std::vector<rvec3_t>& R0, double noise)
{
    const double kr=0.5, r0=1.8, kth=0.15, th0=1.82, kxx=0.02;
    auto e=[=](const std::vector<rvec3_t>& R)
    {
        const rvec3_t a=R[1]-R[0], b=R[2]-R[0];
        const double ra=norm(a), rb=norm(b), th=std::acos((a*b)/(ra*rb)), rxx=norm(R[2]-R[1]);
        return kr*(ra-r0)*(ra-r0)+kr*(rb-r0)*(rb-r0)+kth*(th-th0)*(th-th0)+kxx/rxx;
    };
    auto f=[=](const std::vector<rvec3_t>& R)
    {
        // central differences of e at 1e-6 are exact to ~1e-10 here; the model's point is the optimizer, not the force
        std::vector<rvec3_t> F(3, rvec3_t(0,0,0));
        for (size_t a=0;a<3;++a) for (size_t k=0;k<3;++k)
        {
            auto p=R, m=R; C(p[a],k)+=1e-6; C(m[a],k)-=1e-6;
            C(F[a],k)=-(e(p)-e(m))/2e-6;
        }
        const double pattern[9]={0.3,-0.7,0.2, 0.9,0.1,-0.4, -0.5,0.6,0.8};   // asymmetric: no operation of C2v maps it to itself
        for (size_t a=0;a<3;++a) for (size_t k=0;k<3;++k) C(F[a],k)+=noise*pattern[3*a+k];
        return F;
    };
    return ModelSurface(R0, {8,1,1}, e, f);
}

std::vector<rvec3_t> TriatomicStart() {return {rvec3_t(0,0,0), rvec3_t(1.5,0,1.1), rvec3_t(-1.5,0,1.1)};}

//! C2v with the C2 axis along z and the molecule in the xz plane; atom 0 on the axis, atoms 1 and 2 exchanged by C2 and σ(yz).
std::vector<SymmetryOperation> C2v()
{
    auto diag=[](double x, double y, double z){return rmat3d_t(x,0,0, 0,y,0, 0,0,z);};
    return {{diag( 1, 1,1), {0,1,2}}, {diag(-1,-1,1), {0,2,1}}, {diag(1,-1,1), {0,1,2}}, {diag(-1,1,1), {0,2,1}}};
}

OptimizerOptions Tight()
{
    OptimizerOptions o; o.maxGradient=1e-7; o.rmsGradient=1e-7; o.maxStep=1e-6; o.energyChange=1e-14; o.initialCurvature=1.0;
    return o;
}
}

TEST(GeometryOptimizer, AQuadraticBowlConvergesInAboutNSteps)
{
    // n = 6 (two atoms, all Cartesians), condition number 10, off-diagonal coupling: BFGS learns the Hessian in about n steps.  The start
    // is INSIDE the trust radius (0.2 bohr away), so the count measures BFGS, not the radius: from 1.5 bohr away the radius dominates
    // (14 accepted steps, measured 2026-11-07), which is Rosenbrock's test, not this one.
    const size_t n=6;
    rmat_t A(n, n, 0.0);
    for (size_t i=0;i<n;++i) { A(i,i)=1.0+9.0*double(i)/double(n-1); if (i+1<n) { A(i,i+1)=0.3; A(i+1,i)=0.3; } }
    const std::vector<double> x0={0.1,-0.2,0.3,1.0,0.5,-0.4};
    const std::vector<rvec3_t> R0={rvec3_t(0.18,-0.12,0.25), rvec3_t(1.05,0.43,-0.32)};
    ModelSurface S=Bowl(R0, x0, A);
    const CartesianMap map(S.Reference(), CartesianMap::External::None);
    auto R=GeometryOptimizer(S, map, Tight()).Run();
    ASSERT_TRUE(bool(R)) << R.Error().details;
    for (size_t i=0;i<n;++i) EXPECT_NEAR(R->q[i], x0[i], 1e-6);
    EXPECT_LE(R->Steps(), 2*n) << "BFGS on a quadratic: about n steps (measured " << R->Steps() << "; 11 on 2026-11-07, superlinear to |g| 1e-7 without a line search)";
    EXPECT_EQ(S.evaluations, R->history.size());
}

TEST(GeometryOptimizer, RosenbrockIsWalkedByTheTrustRadius)
{
    // E = (1 − x)² + 100 (y − x²)² + z², from the classic start (−1.2, 1): the curved valley needs the trust radius to shrink and grow.
    auto e=[](const std::vector<rvec3_t>& R){const double x=R[0].x, y=R[0].y, z=R[0].z; return (1-x)*(1-x)+100*(y-x*x)*(y-x*x)+z*z;};
    auto f=[](const std::vector<rvec3_t>& R){const double x=R[0].x, y=R[0].y, z=R[0].z;
                                              return std::vector<rvec3_t>{rvec3_t(2*(1-x)+400*x*(y-x*x), -200*(y-x*x), -2*z)};};
    ModelSurface S({rvec3_t(-1.2,1.0,0.5)}, {1}, e, f);
    const CartesianMap map(S.Reference(), CartesianMap::External::None);
    OptimizerOptions o=Tight(); o.maxSteps=200;
    auto R=GeometryOptimizer(S, map, o).Run();
    ASSERT_TRUE(bool(R)) << R.Error().details;
    EXPECT_NEAR(R->q[0], 1.0, 1e-5);
    EXPECT_NEAR(R->q[1], 1.0, 1e-5);
    EXPECT_NEAR(R->q[2], 0.0, 1e-6);
    size_t limited=0, rejected=0;
    for (const auto& s : R->history) { limited+=s.trustLimited; rejected+=!s.accepted; }
    EXPECT_GT(limited, 0u) << "the radius capped at least one step";
    EXPECT_LE(R->Steps(), 50u) << "35 steps, 4 rejected on 2026-11-07; measured " << R->Steps() << " steps, " << rejected << " rejected";
    for (size_t k=1;k<R->history.size();++k)
        if (R->history[k].accepted) EXPECT_LE(R->history[k].stepLength, R->history[k].trust*(1+1e-12)) << "step " << k << " left its radius";
}

TEST(GeometryOptimizer, AVerySoftModeIsReachedAndTheGradientTestIsHonest)
{
    // Curvatures 1, 0.5 and 1e-4 Ha/bohr² from H₀ = 1: the soft mode moves only once BFGS has learned it.  At convergence the soft
    // coordinate's error is bounded by maxGradient / 1e-4 -- the conditioning lesson: a gradient test is loose along a soft mode.
    rmat_t A(3, 3, 0.0); A(0,0)=1.0; A(1,1)=0.5; A(2,2)=1e-4;
    ModelSurface S=Bowl({rvec3_t(1,1,1)}, {0,0,0}, A);
    const CartesianMap map(S.Reference(), CartesianMap::External::None);
    OptimizerOptions o=Tight(); o.maxGradient=1e-8; o.rmsGradient=1e-8;
    auto R=GeometryOptimizer(S, map, o).Run();
    ASSERT_TRUE(bool(R)) << R.Error().details;
    EXPECT_LE(std::fabs(R->q[2]), o.maxGradient/1e-4*(1+1e-9)) << "the soft coordinate within the gradient test's own bound";
    EXPECT_NEAR(R->q[0], 0.0, 1e-8);
    EXPECT_LE(R->Steps(), 30u) << "23 on 2026-11-07; measured " << R->Steps();
}

TEST(GeometryOptimizer, ASymmetrizedStepStaysOnTheGroupUnderAsymmetricForceNoise)
{
    // Force noise 1e-5 (the model of a grid's), the threshold 10x above it: the calibration margin's own rule.
    const double noise=1e-5;
    OptimizerOptions o; o.maxGradient=1e-4; o.rmsGradient=1e-4;
    {
        ModelSurface S=Triatomic(TriatomicStart(), noise);
        const CartesianMap map(S.Reference(), CartesianMap::External::Rigid, C2v());
        auto R=GeometryOptimizer(S, map, o).Run();
        ASSERT_TRUE(bool(R)) << R.Error().details;
        for (const auto& Rv : S.visited)
        {
            double asym=0.0;
            for (const auto& op : C2v()) for (size_t a=0;a<3;++a) asym=std::max(asym, norm(op.rotation*Rv[a]-Rv[op.perm[a]]));
            EXPECT_LT(asym, 1e-12) << "every visited geometry is on C2v";
        }
        EXPECT_LT(map.Asymmetry(R->q), 1e-12);
    }
    {   // the same noise with no symmetrization: the geometry the optimizer stops at is OFF the group -- the drift the map prevents
        ModelSurface S=Triatomic(TriatomicStart(), noise);
        const CartesianMap plain(S.Reference(), CartesianMap::External::Rigid);
        const CartesianMap probe(S.Reference(), CartesianMap::External::Rigid, C2v());
        auto R=GeometryOptimizer(S, plain, o).Run();
        ASSERT_TRUE(bool(R)) << R.Error().details;
        EXPECT_GT(probe.Asymmetry(R->q), 1e-8) << "without the map the noise breaks the symmetry";
    }
}

TEST(GeometryOptimizer, AThresholdBelowTheForceNoiseFailsHonestly)
{
    // Below the noise the optimizer chases it: the noisy gradient's zero is not E's minimum, steps raise E, the radius collapses.  The run
    // must END IN A FAILURE, never report convergence it did not reach -- which is why the threshold is derived with a margin (§8).
    ModelSurface S=Triatomic(TriatomicStart(), 1e-5);
    const CartesianMap map(S.Reference(), CartesianMap::External::Rigid, C2v());
    OptimizerOptions o; o.maxGradient=1e-7; o.rmsGradient=1e-7; o.maxSteps=100;
    auto R=GeometryOptimizer(S, map, o).Run();
    ASSERT_FALSE(bool(R)) << "converged below the noise it carries";
    EXPECT_TRUE(R.Error().why==OptimizeFailure::Why::TrustCollapsed || R.Error().why==OptimizeFailure::Why::MaxSteps) << R.Error().details;
}

TEST(GeometryOptimizer, TheMolecularMapNeverMovesTheMoleculeAsAWhole)
{
    ModelSurface S=Triatomic({rvec3_t(0.1,0.2,0.0), rvec3_t(1.4,0.3,1.0), rvec3_t(-1.6,0.1,1.3)}, 0.0);
    const CartesianMap map(S.Reference(), CartesianMap::External::Rigid);
    OptimizerOptions o; o.maxGradient=1e-7; o.rmsGradient=1e-7;
    auto R=GeometryOptimizer(S, map, o).Run();
    ASSERT_TRUE(bool(R)) << R.Error().details;
    for (size_t k=1;k<S.visited.size();++k)
    {
        const auto& A=S.visited[k-1]; const auto& B=S.visited[k];
        rvec3_t t(0,0,0); for (size_t a=0;a<3;++a) t+=B[a]-A[a];
        EXPECT_LT(norm(t), 1e-10) << "evaluation " << k << ": no net translation";
    }
    // the converged shape: both bonds equal and the angle settled (the X–X repulsion opens it past th0)
    const auto& g=R->geometry.R;
    EXPECT_NEAR(norm(g[1]-g[0]), norm(g[2]-g[0]), 1e-5);
}

TEST(GeometryOptimizer, RefusalsAndFailuresAreOutcomes)
{
    rmat_t A(3, 3, 0.0); A(0,0)=A(1,1)=A(2,2)=1.0;
    {   // a point that fails is the optimizer's failure, naming the step, carrying the point's own reason
        ModelSurface S=Bowl({rvec3_t(1,1,1)}, {0,0,0}, A);
        S.failAt=3;
        const CartesianMap map(S.Reference(), CartesianMap::External::None);
        auto R=GeometryOptimizer(S, map, Tight()).Run();
        ASSERT_FALSE(bool(R));
        EXPECT_EQ(R.Error().why, OptimizeFailure::Why::PES);
        ASSERT_TRUE(R.Error().pes.has_value());
        EXPECT_EQ(R.Error().pes->why, PESFailure::Why::SCF);
        EXPECT_NE(R.Error().details.find("step 2"), std::string::npos) << R.Error().details;
    }
    {   // the step cap
        rmat_t Ar(3, 3, 0.0); Ar(0,0)=1.0; Ar(1,1)=1e-3; Ar(2,2)=1e-6;
        ModelSurface S=Bowl({rvec3_t(1,1,1)}, {0,0,0}, Ar);
        const CartesianMap map(S.Reference(), CartesianMap::External::None);
        OptimizerOptions o=Tight(); o.maxSteps=2;
        auto R=GeometryOptimizer(S, map, o).Run();
        ASSERT_FALSE(bool(R));
        EXPECT_EQ(R.Error().why, OptimizeFailure::Why::MaxSteps);
    }
    {   // the map refuses what it cannot represent
        EXPECT_THROW(CartesianMap(PESGeometry{std::nullopt, {}, {}}), std::invalid_argument);
        EXPECT_THROW(CartesianMap(PESGeometry{rmat3d_t(1,0,0, 0,1,0, 0,0,1), {rvec3_t(0,0,0)}, {1}}), std::invalid_argument);
        auto off=TriatomicStart(); off[1].x+=0.1;
        EXPECT_THROW(CartesianMap(PESGeometry{std::nullopt, off, {8,1,1}}, CartesianMap::External::Rigid, C2v()), std::invalid_argument)
            << "a reference off the stated group";
    }
    {   // a start already at the minimum: no step
        ModelSurface S=Bowl({rvec3_t(0,0,0)}, {0,0,0}, A);
        const CartesianMap map(S.Reference(), CartesianMap::External::None);
        auto R=GeometryOptimizer(S, map, Tight()).Run();
        ASSERT_TRUE(bool(R));
        EXPECT_EQ(R->Steps(), 0u);
        EXPECT_EQ(S.evaluations, 1u);
    }
}

TEST(GeometryOptimizer, TheForceThresholdIsDerivedFromEpsThroughANamedCalibration)
{
    const ForceThresholdCalibration cal;
    EXPECT_DOUBLE_EQ(cal.eta, 15.0);
    EXPECT_DOUBLE_EQ(cal.margin, 3.0);
    const OptimizerOptions o=OptimizerOptions::FromEps(1e-5, cal);
    EXPECT_NEAR(o.maxGradient, 4.5e-4, 1e-15) << "eps 1e-5: the conventional default threshold";
    EXPECT_NEAR(o.rmsGradient, 3.0e-4, 1e-15);
    EXPECT_NEAR(o.maxStep,     1.8e-3, 1e-15);
    EXPECT_DOUBLE_EQ(o.energyChange, 1e-5);
    EXPECT_THROW(OptimizerOptions::FromEps(0.0), std::invalid_argument);
}
