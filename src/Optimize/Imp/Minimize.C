// File: Optimize/Imp/Minimize.C  The Minimize driver and the conventional stopping rule.
module;
#include <cmath>
#include <cstdio>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.Optimize.Minimize;
import qchem.Blaze;
import qchem.Reporting;

namespace qchem::Optimize
{

namespace
{
std::string Sci(double x) {char b[32]; std::snprintf(b, sizeof b, "%.3g", x); return b;}
double MaxAbs(const rvec_t& v) {double m=0.0; for (double x : v) m=std::max(m, std::fabs(x)); return m;}
double Rms(const rvec_t& v)    {return v.size() ? blazem::norm(v)/std::sqrt(double(v.size())) : 0.0;}
}

size_t Minimum::Steps() const
{
    size_t n=0;
    for (size_t k=1;k<history.size();++k) if (history[k].accepted) ++n;
    return n;
}

double RoundoffNoise(double value) {return kRoundoffUnits*std::numeric_limits<double>::epsilon()*std::fabs(value);}

//------------------------------------------------------------------------------------------------------------------------------------
// GradientStepValue

GradientStepValue::GradientStepValue(double g, double rms, double s, double dv, size_t n)
    : maxGradient(g), rmsGradient(rms), maxStep(s), valueChange(dv), maxIterations(n)
{
    if (!(g>0.0 && rms>0.0 && s>0.0 && dv>=0.0))
        throw std::invalid_argument("GradientStepValue: the gradient and step thresholds must be > 0 and the value change >= 0");
}

StoppingRule::Decision GradientStepValue::Judge(const IterationState& s) const
{
    if (!s.gradient) throw std::logic_error("GradientStepValue::Judge: the state carries no gradient");
    const bool gradient = MaxAbs(*s.gradient)<=maxGradient && Rms(*s.gradient)<=rmsGradient;
    const bool moved    = !s.step || MaxAbs(*s.step)<=maxStep || std::fabs(s.valueChange)<=valueChange;
    if (gradient && moved)       return Decision::Converged;
    if (s.index>=maxIterations)  return Decision::Budget;
    return Decision::Continue;
}

//------------------------------------------------------------------------------------------------------------------------------------
// Minimize

Outcome<Minimum, Failure> Minimize(Cost& cost, Stepper& stepper, const StoppingRule& rule, rvec_t q)
{
    using O=Outcome<Minimum, Failure>;
    using Decision=StoppingRule::Decision;
    namespace rpt=qchem::report;
    const size_t n=cost.Dimension();
    if (q.size()!=n) throw std::invalid_argument("Minimize: the start has "+std::to_string(q.size())+" coordinates for the cost's dimension "+std::to_string(n));
    std::vector<Iteration> history;
    auto fail=[&](Failure::Why why, std::string details) {return O::Fail(Failure{why, std::move(details), history});};
    if (auto why=stepper.Refuse(cost)) return fail(Failure::Why::MissingCapability, stepper.Name()+": "+*why);

    const auto* projected=dynamic_cast<const Projected*>(&cost);
    const auto* noisy    =dynamic_cast<const ValueNoise*>(&cost);
    auto noise=[&](const rvec_t& x, double v) {return noisy ? noisy->Noise(x, v) : RoundoffNoise(v);};
    // One evaluation, the gradient projected onto the step space.  The CALLER holds the report row open around it (the cost's facts
    // land in it) and around the stepper's verdict.
    auto evaluate=[&](const rvec_t& x, Iteration row, std::string& why)->std::optional<CostValue>
    {
        rpt::Set("step", (long)row.index);
        auto e=cost.Evaluate(x, Want{});
        if (!e) { why=e.Error().details; row.accepted=false; rpt::Set("failed", why); history.push_back(row); return std::nullopt; }
        CostValue v=e.TakeValue();
        if (v.gradient.size()!=n) throw std::logic_error("Minimize: the cost returned a gradient of "+std::to_string(v.gradient.size())+" for dimension "+std::to_string(n));
        if (projected) v.gradient=projected->Project(x, v.gradient);
        row.value=v.value; row.maxGradient=MaxAbs(v.gradient); row.rmsGradient=Rms(v.gradient);
        rpt::Set("value", v.value); rpt::Set("maxGradient", row.maxGradient); rpt::Set("rmsGradient", row.rmsGradient);
        rpt::Set("stepLength", row.stepLength); rpt::Set("bound", row.bound); rpt::Set("limited", row.limited);
        history.push_back(row);
        return v;
    };
    auto setAccepted=[&](bool a) {history.back().accepted=a; rpt::Set("accepted", a);};
    auto done=[&](const rvec_t& x, double v) {return O::Ok(Minimum{x, v, history});};

    std::string why;
    std::optional<CostValue> here;
    { rpt::Row r("optimizer"); here=evaluate(q, Iteration{}, why); if (here) setAccepted(true); }
    if (!here) return fail(Failure::Why::Cost, "the starting point: "+why);
    switch (rule.Judge(IterationState{0, here->value, &here->gradient, nullptr, 0.0}))
    {
        case Decision::Converged: return done(q, here->value);
        case Decision::Budget:    return fail(Failure::Why::MaxSteps, "the rule's budget allows no step");
        case Decision::Continue:  break;
    }
    stepper.Start(cost, q, *here);
    for (size_t index=1;;++index)
    {
        const Proposal p=stepper.Propose(q, *here);
        if (p.step.size()==0) return done(q, here->value);   // nothing may move
        if (p.step.size()!=n) throw std::logic_error("Minimize: "+stepper.Name()+" proposed a step of "+std::to_string(p.step.size())+" for dimension "+std::to_string(n));
        Iteration row; row.index=index; row.stepLength=blazem::norm(p.step); row.bound=p.bound; row.limited=p.limited;
        const rvec_t q1=q+p.step;
        rpt::Row r("optimizer");
        auto there=evaluate(q1, row, why);
        if (!there) return fail(Failure::Why::Cost, "step "+std::to_string(index)+": "+why);
        const double change=there->value-here->value;
        const Verdict v=stepper.Tell(p.step, *there, noise(q, here->value)+noise(q1, there->value));
        const bool accepted = v.kind==Verdict::Kind::Accept;
        setAccepted(accepted);
        // RULE 3 (§4): every evaluated point is judged, accepted or not.
        const Decision d=rule.Judge(IterationState{index, there->value, &there->gradient, &p.step, change});
        if (accepted)
        {
            q=q1; here=std::move(there);
            if (d==Decision::Converged) return done(q, here->value);
        }
        else if (d==Decision::Converged)
            return fail(Failure::Why::GradientDisagrees, "step "+std::to_string(index)+": the gradient meets the rule (max "
                        +Sci(history.back().maxGradient)+") at a trial "+Sci(change)+" above the incumbent, beyond the "
                        "value noise: the gradient's zero is not the value's minimum");
        if (v.kind==Verdict::Kind::Collapsed) return fail(Failure::Why::TrustCollapsed, "step "+std::to_string(index)+": "+v.details);
        if (d==Decision::Budget) return fail(Failure::Why::MaxSteps, "not converged after "+std::to_string(index)+" steps");
    }
}

} // namespace qchem::Optimize
