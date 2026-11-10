// File: Optimize/tests/Minimize.C  The Minimize driver and its steppers on ANALYTIC model costs (moved down from
// Calculation/tests/GeometryOptimizer.C, V-OPTIMIZE-DIP §8 step 1; CLAUDE.md: unit over integration).  Each model IS a cost of q
// (doc/Records/VOptimizeDIP_Design.md §7: inheritance is right for a class whose natural parameters are q), so what is tested is the
// algorithm and the driver alone: step counts on a quadratic, the trust radius on Rosenbrock, conditioning on a soft mode, the noise
// rules of §4, the capability faces, and the refusals.
#include "gtest/gtest.h"
#include <cmath>
#include <functional>
#include <memory>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>
import qchem.Optimize;
import qchem.Blaze;

using namespace qchem;
using namespace qchem::Optimize;

namespace
{
//! A model cost from closures: counts its evaluations, records every q, and can fail at a chosen evaluation.
class Model : public virtual Cost
{
public:
    using V_t = std::function<double(const rvec_t&)>;
    using G_t = std::function<rvec_t(const rvec_t&)>;
    Model(size_t n, V_t v, G_t g) : itsN(n), itsV(std::move(v)), itsG(std::move(g)) {}
    size_t Dimension() const override {return itsN;}
    Outcome<CostValue, CostFailure> Evaluate(const rvec_t& q, const Want& want) override
    {
        ++evaluations; visited.push_back(q);
        if (failAt && evaluations==*failAt) return Outcome<CostValue, CostFailure>::Fail({"model failure"});
        CostValue c; c.value=itsV(q);
        if (want.gradient) c.gradient=itsG(q);
        return Outcome<CostValue, CostFailure>::Ok(c);
    }
    size_t evaluations=0;
    std::optional<size_t> failAt;
    std::vector<rvec_t> visited;
private:
    size_t itsN;
    V_t itsV; G_t itsG;
};

//! E = ½ (q − q0)ᵀ A (q − q0) + cᵀ q-noise: \a gradientNoise adds a FIXED vector to the gradient only (the model of a grid's force
//! noise: not the gradient of E).
Model Bowl(const rmat_t& A, const rvec_t& q0, const rvec_t& gradientNoise={})
{
    auto v=[=](const rvec_t& q){const rvec_t d=q-q0; return 0.5*blazem::dot(d, rvec_t(A*d));};
    auto g=[=](const rvec_t& q){rvec_t r=A*rvec_t(q-q0); if (gradientNoise.size()) r+=gradientNoise; return r;};
    return Model(q0.size(), v, g);
}

rmat_t Diag(std::vector<double> d) {rmat_t A(d.size(), d.size(), 0.0); for (size_t i=0;i<d.size();++i) A(i,i)=d[i]; return A;}
rvec_t Vec(std::vector<double> x)  {rvec_t v(x.size()); for (size_t i=0;i<x.size();++i) v[i]=x[i]; return v;}

//! The tight rule of the moved tests: |g| 1e-7, step 1e-6, value change 1e-14.
GradientStepValue Tight(size_t maxIterations=200, double g=1e-7) {return GradientStepValue(g, g, 1e-6, 1e-14, maxIterations);}

BFGS_RFO::Options UnitCurvature() {BFGS_RFO::Options o; o.initialCurvature=1.0; return o;}

//! THE LSP PAIR (V-OPTIMIZE-DIP §8 step 3): every model-cost test below runs, unchanged, under each stepper.
enum class Algo {BFGS_RFO, L_BFGS};
std::string Name(Algo a) {return a==Algo::BFGS_RFO ? "BFGS_RFO" : "L_BFGS";}
//! gtest prints the parameter with this (ctest names the case by it, not by its bytes).
void PrintTo(Algo a, std::ostream* os) {*os << Name(a);}
std::unique_ptr<Stepper> Make(Algo a, double initialCurvature=0.5)
{
    if (a==Algo::BFGS_RFO) { BFGS_RFO::Options o; o.initialCurvature=initialCurvature; return std::make_unique<BFGS_RFO>(o); }
    L_BFGS::Options o; o.initialCurvature=initialCurvature; return std::make_unique<L_BFGS>(o);
}

//! The 6-d coupled bowl, condition number 10.
rmat_t Coupled6()
{
    const size_t n=6;
    rmat_t A(n, n, 0.0);
    for (size_t i=0;i<n;++i) { A(i,i)=1.0+9.0*double(i)/double(n-1); if (i+1<n) { A(i,i+1)=0.3; A(i+1,i)=0.3; } }
    return A;
}
}

class BothSteppers : public ::testing::TestWithParam<Algo>
{
protected:
    std::unique_ptr<Stepper> Make(double initialCurvature=0.5) const {return ::Make(GetParam(), initialCurvature);}
    //! A step-count bound per algorithm (the counts differ; the claim, convergence on the same cost under the same rule, does not).
    size_t Bound(size_t bfgs, size_t lbfgs) const {return GetParam()==Algo::BFGS_RFO ? bfgs : lbfgs;}
};

TEST_P(BothSteppers, AQuadraticBowlConvergesInAboutNSteps)
{
    // n = 6, condition number 10, off-diagonal coupling: BFGS learns the Hessian in about n steps.  The start is INSIDE the trust radius
    // (0.2 away), so the count measures BFGS, not the radius.
    const rvec_t q0=Vec({0.1,-0.2,0.3,1.0,0.5,-0.4});
    Model S=Bowl(Coupled6(), q0);
    auto stepper=Make(1.0);
    auto R=Minimize(S, *stepper, Tight(), Vec({0.18,-0.12,0.25,1.05,0.43,-0.32}));
    ASSERT_TRUE(bool(R)) << R.Error().details;
    for (size_t i=0;i<6;++i) EXPECT_NEAR(R->q[i], q0[i], 1e-6);
    EXPECT_LE(R->Steps(), Bound(12, 14)) << "about n steps to |g| 1e-7 without a line search (2026-11-09: BFGS-RFO 11, L-BFGS 12); measured " << R->Steps();
    EXPECT_EQ(S.evaluations, R->history.size());
}

TEST_P(BothSteppers, RosenbrockIsWalkedByTheTrustRadius)
{
    // E = (1 − x)² + 100 (y − x²)² + z², from the classic start (−1.2, 1): the curved valley needs the radius to shrink and grow.
    auto v=[](const rvec_t& q){const double x=q[0], y=q[1], z=q[2]; return (1-x)*(1-x)+100*(y-x*x)*(y-x*x)+z*z;};
    auto g=[](const rvec_t& q){const double x=q[0], y=q[1], z=q[2]; return Vec({-2*(1-x)-400*x*(y-x*x), 200*(y-x*x), 2*z});};
    Model S(3, v, g);
    auto stepper=Make(1.0);
    auto R=Minimize(S, *stepper, Tight(), Vec({-1.2,1.0,0.5}));
    ASSERT_TRUE(bool(R)) << R.Error().details;
    EXPECT_NEAR(R->q[0], 1.0, 1e-5);
    EXPECT_NEAR(R->q[1], 1.0, 1e-5);
    EXPECT_NEAR(R->q[2], 0.0, 1e-6);
    size_t limited=0, rejected=0;
    for (const auto& s : R->history) { limited+=s.limited; rejected+=!s.accepted; }
    EXPECT_GT(limited, 0u) << "the radius capped at least one step";
    EXPECT_LE(R->Steps(), Bound(50, 60)) << "2026-11-09: BFGS-RFO 35 steps, 4 rejected; L-BFGS 48, 6 rejected; measured " << R->Steps() << " steps, " << rejected << " rejected";
    for (size_t k=1;k<R->history.size();++k)
        EXPECT_LE(R->history[k].stepLength, R->history[k].bound*(1+1e-12)) << "step " << k << " left its radius";
}

TEST_P(BothSteppers, AVerySoftModeIsReachedAndTheGradientTestIsHonest)
{
    // Curvatures 1, 0.5 and 1e-4 from H₀ = 1: the soft mode moves only once BFGS has learned it.  At convergence the soft coordinate's
    // error is bounded by maxGradient / 1e-4 -- the conditioning lesson: a gradient test is loose along a soft mode.
    Model S=Bowl(Diag({1.0, 0.5, 1e-4}), Vec({0,0,0}));
    auto stepper=Make(1.0);
    const GradientStepValue rule=Tight(200, 1e-8);
    auto R=Minimize(S, *stepper, rule, Vec({1,1,1}));
    ASSERT_TRUE(bool(R)) << R.Error().details;
    EXPECT_LE(std::fabs(R->q[2]), rule.maxGradient/1e-4*(1+1e-9)) << "the soft coordinate within the gradient test's own bound";
    EXPECT_NEAR(R->q[0], 0.0, 1e-8);
    EXPECT_LE(R->Steps(), Bound(30, 30)) << "2026-11-09: BFGS-RFO 23, L-BFGS 19; measured " << R->Steps();
}

TEST_P(BothSteppers, AThresholdBelowTheGradientNoiseBuysNoPrecision)
{
    // A fixed 1e-5 offset c on the gradient only: its zero, x* = −A⁻¹c, is not the value's minimum (0).  Above the noise (threshold 1e-4,
    // the calibration margin's own rule) the run converges.  Below it (1e-7) the PRECISION the threshold asks for is not available, and
    // no stepper may pretend otherwise: it either ENDS IN A FAILURE (BFGS-RFO, 2026-11-09: GradientDisagrees at step 8, the trial that met
    // the rule read above the incumbent), or it lands on x* (L-BFGS, 2026-11-09: converged in 8 steps AT x*, 1.2e-5 from the minimum,
    // 40x the 3e-7 the threshold would imply).  Which one is the path's luck, not a property of the algorithm, so the claim is the pair.
    // That is why the geometry client derives its threshold with a margin above the grid's force noise (ForceThresholdCalibration).
    const rmat_t A=Diag({1.0, 0.6, 0.3});
    const rvec_t noise=Vec({0.3e-5,-0.7e-5,0.2e-5});
    const rvec_t xStar=Vec({-0.3e-5/1.0, 0.7e-5/0.6, -0.2e-5/0.3});
    {
        Model S=Bowl(A, Vec({0,0,0}), noise);
        auto stepper=Make();
        auto R=Minimize(S, *stepper, GradientStepValue(1e-4, 1e-4, 1.8e-3, 1e-6, 100), Vec({0.5,-0.4,0.3}));
        ASSERT_TRUE(bool(R)) << R.Error().details;
    }
    {
        Model S=Bowl(A, Vec({0,0,0}), noise);
        auto stepper=Make();
        auto R=Minimize(S, *stepper, GradientStepValue(1e-7, 1e-7, 1e-6, 1e-14, 100), Vec({0.5,-0.4,0.3}));
        if (R)
        {
            EXPECT_LT(blazem::norm(rvec_t(R->q-xStar)), 1e-6) << "a 'converged' answer is the NOISY gradient's zero";
            EXPECT_GT(blazem::norm(R->q), 1e-5) << "... which is not the minimum: the threshold bought no precision";
        }
        else
        {
            const auto why=R.Error().why;
            EXPECT_TRUE(why==Failure::Why::GradientDisagrees || why==Failure::Why::TrustCollapsed || why==Failure::Why::MaxSteps) << R.Error().details;
        }
    }
}

namespace
{
//! The soft-mode bowl (curvatures 1, 0.5, 1e-4) whose VALUE reads high by up to \a amp, at random (a hash of q and the evaluation
//! count): the model of an SCF stopped at its ΔE gate, which is VARIATIONAL, so it reads ABOVE its converged value.  The gradient is
//! exact.  This class does NOT state its noise.
class NoisyValueBowl : public virtual Cost
{
public:
    explicit NoisyValueBowl(double amp) : itsAmp(amp) {}
    size_t Dimension() const override {return 3;}
    Outcome<CostValue, CostFailure> Evaluate(const rvec_t& q, const Want&) override
    {
        ++evaluations;
        double h=0.0; for (size_t i=0;i<3;++i) h+=std::sin(1e4*q[i]*(i+1.37));
        const double u=0.5*(1.0+std::sin(7.1*h+double(evaluations)));
        CostValue c;
        c.value=0.5*(q[0]*q[0]+0.5*q[1]*q[1]+1e-4*q[2]*q[2])+itsAmp*u;
        c.gradient=Vec({q[0], 0.5*q[1], 1e-4*q[2]});
        return Outcome<CostValue, CostFailure>::Ok(c);
    }
    size_t evaluations=0;
protected:
    double itsAmp;
};
//! The same cost, stating its bound.
class DeclaredNoisyValueBowl : public NoisyValueBowl, public virtual ValueNoise
{
public:
    using NoisyValueBowl::NoisyValueBowl;
    double Noise(const rvec_t&, double) const override {return itsAmp;}
};
}

TEST_P(BothSteppers, ARiseWithinTheValueNoiseIsNotARise)
{
    // §4 rules 1-3 (lane H's rutile finding, 2026-11-09): the value reads high by up to 1e-6, the gradient is exact.  Once the true
    // change per step falls below the noise, the accepted points are the ones that happened to read LOW, so most trials "rise"; the old
    // rule (reject any rise, cut the radius by 4) then walks the radius to trustMin at a point the gradient says is nearly converged.
    // Measured 2026-11-09: declared, CONVERGED in 23 steps; the same cost without the face, TrustCollapsed at step 30 (at noise 1e-8:
    // GradientDisagrees at step 23, the gradient met the rule at a trial that read above the incumbent).
    const GradientStepValue rule(1e-7, 1e-7, 1e-6, 0.0, 200);
    {
        DeclaredNoisyValueBowl S(1e-6);
        auto stepper=Make(1.0);
        auto R=Minimize(S, *stepper, rule, Vec({1,1,1}));
        ASSERT_TRUE(bool(R)) << R.Error().details;
        EXPECT_NEAR(R->q[0], 0.0, 1e-7);
        EXPECT_LE(std::fabs(R->q[2]), 1e-7/1e-4*(1+1e-9));
    }
    {
        NoisyValueBowl S(1e-6);
        auto stepper=Make(1.0);
        auto R=Minimize(S, *stepper, rule, Vec({1,1,1}));
        ASSERT_FALSE(bool(R)) << "the undeclared twin converged: the test no longer shows what the face is for";
        EXPECT_EQ(R.Error().why, Failure::Why::TrustCollapsed) << R.Error().details;
    }
}

namespace
{
//! A bowl that holds coordinate 2 fixed (its step space is coordinates 0 and 1), and optionally offers its exact Hessian.
class HeldBowl : public virtual Cost, public virtual Projected
{
public:
    explicit HeldBowl(rmat_t A) : itsA(std::move(A)) {}
    size_t Dimension() const override {return itsA.rows();}
    Outcome<CostValue, CostFailure> Evaluate(const rvec_t& q, const Want&) override
    {
        ++evaluations; visited.push_back(q);
        return Outcome<CostValue, CostFailure>::Ok(CostValue{0.5*blazem::dot(q, rvec_t(itsA*q)), itsA*q});
    }
    rvec_t Project(const rvec_t&, const rvec_t& v) const override {rvec_t p=v; p[2]=0.0; return p;}
    size_t evaluations=0;
    std::vector<rvec_t> visited;
protected:
    rmat_t itsA;
};
class ExactHessianBowl : public virtual Cost, public virtual ModelHessian
{
public:
    explicit ExactHessianBowl(rmat_t A) : itsA(std::move(A)) {}
    size_t Dimension() const override {return itsA.rows();}
    Outcome<CostValue, CostFailure> Evaluate(const rvec_t& q, const Want&) override
    {
        ++evaluations;
        return Outcome<CostValue, CostFailure>::Ok(CostValue{0.5*blazem::dot(q, rvec_t(itsA*q)), itsA*q});
    }
    rvec_t ApplyH0(const rvec_t&, const rvec_t& v) const override {return itsA*v;}
    size_t evaluations=0;
private:
    rmat_t itsA;
};
}

TEST_P(BothSteppers, TheCostsStepSpaceHoldsWhatItFixes)
{
    // Coordinate 2 is held: no evaluated point moves it, and the rule judges the PROJECTED gradient (whose held component, 0.7, is
    // never zero) -- so the run converges.
    rmat_t A=Diag({1.0, 2.0, 3.0}); A(0,2)=A(2,0)=0.2;
    HeldBowl S(A);
    auto stepper=Make(1.0);
    auto R=Minimize(S, *stepper, Tight(), Vec({0.5,-0.4,0.7}));
    ASSERT_TRUE(bool(R)) << R.Error().details;
    for (const auto& q : S.visited) EXPECT_EQ(q[2], 0.7);
    EXPECT_GT(R->Steps(), 0u);
}

TEST(Minimize, AModelHessianIsUsedWhenTheCostOffersOne)
{
    // The exact Hessian as the model: RFO on a quadratic with the exact H steps nearly to Newton's from the start.
    const rmat_t A=Coupled6();
    const rvec_t start=Vec({0.18,-0.12,0.25,1.05,0.43,-0.32});
    ExactHessianBowl exact(A);
    BFGS_RFO a(UnitCurvature());
    auto R=Minimize(exact, a, Tight(), start);
    ASSERT_TRUE(bool(R)) << R.Error().details;
    Model plain=Bowl(A, Vec({0,0,0,0,0,0}));
    BFGS_RFO b(UnitCurvature());
    auto P=Minimize(plain, b, Tight(), start);
    ASSERT_TRUE(bool(P)) << P.Error().details;
    EXPECT_LT(R->Steps(), P->Steps()) << "with the exact Hessian " << R->Steps() << " steps, from c·1 " << P->Steps();
    EXPECT_LE(R->Steps(), 7u) << "6 on 2026-11-09 (RFO with the exact H is not the Newton step); measured " << R->Steps();
}

namespace
{
struct TypedFailure { int code=0; std::string details; };
//! A typed model cost (the PESCost shape): its own failure type rides back through Minimize<Fail>.
class TypedBowl : public virtual tCost<TypedFailure>
{
public:
    size_t Dimension() const override {return 3;}
    Outcome<CostValue, TypedFailure> EvaluateTyped(const rvec_t& q, const Want&) override
    {
        if (++evaluations==3) return Outcome<CostValue, TypedFailure>::Fail({42, "typed model failure"});
        return Outcome<CostValue, TypedFailure>::Ok(CostValue{0.5*blazem::dot(q, q), q});
    }
    size_t evaluations=0;
};
//! A stepper that needs a face no model offers.
class NeedsLeastSquares : public virtual Stepper
{
public:
    std::string                Name() const override {return "test-LM";}
    std::optional<std::string> Refuse(const Cost&) const override {return "needs LeastSquares; this cost offers Cost";}
    void     Start  (const Cost&, const rvec_t&, const CostValue&) override {}
    Proposal Propose(const rvec_t&, const CostValue&) override {return {};}
    Verdict  Tell   (const rvec_t&, const CostValue&, double) override {return {};}
};
}

TEST(Minimize, RefusalsAndFailuresAreOutcomes)
{
    {   // an evaluation that fails is the minimizer's failure, naming the step; a typed cost gets its own failure back by value
        TypedBowl S;
        BFGS_RFO bfgs(UnitCurvature());
        auto R=Minimize(S, bfgs, Tight(), Vec({1,1,1}));
        ASSERT_FALSE(bool(R));
        EXPECT_EQ(R.Error().why, Failure::Why::Cost);
        ASSERT_TRUE(R.Error().cost.has_value());
        EXPECT_EQ(R.Error().cost->code, 42);
        EXPECT_NE(R.Error().details.find("step 2"), std::string::npos) << R.Error().details;
    }
    {   // the rule's budget
        Model S=Bowl(Diag({1.0, 1e-3, 1e-6}), Vec({0,0,0}));
        BFGS_RFO bfgs(UnitCurvature());
        auto R=Minimize(S, bfgs, Tight(2), Vec({1,1,1}));
        ASSERT_FALSE(bool(R));
        EXPECT_EQ(R.Error().why, Failure::Why::MaxSteps);
        EXPECT_EQ(S.evaluations, 3u);
    }
    {   // a stepper that needs a missing face refuses BEFORE any evaluation
        Model S=Bowl(Diag({1.0, 1.0, 1.0}), Vec({0,0,0}));
        NeedsLeastSquares lm;
        auto R=Minimize(S, lm, Tight(), Vec({1,1,1}));
        ASSERT_FALSE(bool(R));
        EXPECT_EQ(R.Error().why, Failure::Why::MissingCapability);
        EXPECT_NE(R.Error().details.find("LeastSquares"), std::string::npos) << R.Error().details;
        EXPECT_EQ(S.evaluations, 0u);
    }
    {   // a start already at the minimum: no step
        Model S=Bowl(Diag({1.0, 1.0, 1.0}), Vec({0,0,0}));
        BFGS_RFO bfgs;
        auto R=Minimize(S, bfgs, Tight(), Vec({0,0,0}));
        ASSERT_TRUE(bool(R));
        EXPECT_EQ(R->Steps(), 0u);
        EXPECT_EQ(S.evaluations, 1u);
    }
    {   // inconsistent inputs are broken invariants: they throw
        Model S=Bowl(Diag({1.0, 1.0, 1.0}), Vec({0,0,0}));
        BFGS_RFO bfgs;
        EXPECT_THROW((void)Minimize(S, bfgs, Tight(), Vec({1,1})), std::invalid_argument);
        BFGS_RFO::Options o; o.trust.min=0.0;
        EXPECT_THROW(BFGS_RFO{o}, std::invalid_argument);
        EXPECT_THROW(GradientStepValue(0.0, 1e-4, 1e-3, 1e-6, 10), std::invalid_argument);
    }
}

TEST(Minimize, TheConventionalRuleJudgesGradientThenStepOrValue)
{
    const GradientStepValue rule(1e-4, 1e-4, 1e-3, 1e-6, 5);
    const rvec_t small=Vec({1e-5, -1e-5}), big=Vec({1e-3, 0.0}), tiny=Vec({1e-4, 0.0}), long_=Vec({1e-2, 0.0});
    using D=StoppingRule::Decision;
    EXPECT_EQ(rule.Judge({0, 1.0, &small, nullptr, 0.0}),   D::Converged) << "at the start the gradient alone decides";
    EXPECT_EQ(rule.Judge({1, 1.0, &small, &tiny, 1e-3}),    D::Converged) << "a short step";
    EXPECT_EQ(rule.Judge({1, 1.0, &small, &long_, 1e-7}),   D::Converged) << "a long step with a small value change";
    EXPECT_EQ(rule.Judge({1, 1.0, &small, &long_, 1e-3}),   D::Continue)  << "a long step and a large change";
    EXPECT_EQ(rule.Judge({1, 1.0, &big, &tiny, 0.0}),       D::Continue)  << "the gradient fails";
    EXPECT_EQ(rule.Judge({5, 1.0, &big, &tiny, 0.0}),       D::Budget);
}

INSTANTIATE_TEST_SUITE_P(Minimize, BothSteppers, ::testing::Values(Algo::BFGS_RFO, Algo::L_BFGS),
                         [](const ::testing::TestParamInfo<Algo>& i){return Name(i.param);});
