// File: Optimize/Minimize.C  The minimizer: an ask/tell STEPPER per algorithm, a STOPPING RULE from the client, and the \c Minimize
// driver between them (V-OPTIMIZE-DIP, doc/Records/VOptimizeDIP_Design.md §2-§4; user ruling 2026-11-09, decisions 3 and 4).
//
// WHY A STEPPER.  Each algorithm is written ONCE as ask/tell: \c Propose a step from where it stands, be \c Told what the cost said at
// the trial.  \c Minimize is the driver for a client that hands over its loop (geometry, the grid tool); a client that owns its loop
// (the SCF's GDM: MOM, the Ladder, held fills) holds a Stepper and calls it from inside that loop.  So any algorithm serves both.
//
// WHO OWNS WHAT.
//   - the STEPPER owns the algorithm's knobs and state (trust radii, the Hessian, a memory), set in its constructor;
//   - the client's STOPPING RULE owns convergence and the budget (the geometry client derives its thresholds from ε through a named
//     calibration that stays in Calculation);
//   - the DRIVER owns the sequence: evaluate, project the gradient, judge EVERY evaluated point by the rule (§4 rule 3), tell the
//     stepper, report.
//
// REPORTING (CLAUDE.md): the driver writes one \c optimizer row per evaluation at the caller's cursor, and evaluates INSIDE that row,
// so a cost's own facts (a PES point's SCF) land in the row of the step that asked for them.
module;
#include <limits>
#include <optional>
#include <string>
#include <vector>
export module qchem.Optimize.Minimize;
export import qchem.Optimize.Cost;

export namespace qchem::Optimize
{

//! One evaluation in the history (a VALUE: the report and the tests read it).
struct Iteration
{
    size_t index       = 0;
    double value       = 0.0;
    double maxGradient = 0.0;    //!< of the PROJECTED gradient
    double rmsGradient = 0.0;
    double stepLength  = 0.0;    //!< |s| of the step that led here (0 at the start)
    double bound       = 0.0;    //!< the bound the stepper held the step to (a trust radius); 0: none
    bool   accepted    = true;   //!< false: the stepper rejected the trial and stayed where it was
    bool   limited     = false;  //!< the step was cut to \c bound
};

struct Minimum
{
    rvec_t                 q;
    double                 value = 0.0;
    std::vector<Iteration> history;   //!< every evaluation in order, the start first
    size_t Steps() const;             //!< accepted steps
};

struct Failure
{
    //! \c Cost: an evaluation failed.  \c MaxSteps: the rule's budget ran out.  \c TrustCollapsed: the stepper could not find a step that
    //! the values accept, while the gradient did NOT meet the rule.  \c GradientDisagrees: a trial's gradient met the rule but its value
    //! sits above the incumbent's by more than the noise: the gradient's zero is not the value's minimum (gradient noise above the
    //! threshold, or a gradient that is not the value's).  \c MissingCapability: the stepper needs a face the cost does not offer.
    enum class Why {Cost, MaxSteps, TrustCollapsed, GradientDisagrees, MissingCapability};
    Why                    why = Why::Cost;
    std::string            details;
    std::vector<Iteration> history;
};

//! The client's typed twin of \c Failure: the cost's own failure, by value, when \c why is \c Cost.
template <class Fail> struct MinimizeFailure : Failure
{
    std::optional<Fail> cost;
};

//! What a stopping rule sees of one evaluated point.
struct IterationState
{
    size_t        index       = 0;
    double        value       = 0.0;
    const rvec_t* gradient    = nullptr;   //!< PROJECTED
    const rvec_t* step        = nullptr;   //!< the step that led here; null at the start
    double        valueChange = 0.0;       //!< value − the incumbent's (0 at the start)
};

//! \brief The CLIENT's stopping rule.  \c Converged: stop with this point.  \c Budget: stop, unconverged.
class StoppingRule
{
public:
    enum class Decision {Continue, Converged, Budget};
    virtual ~StoppingRule() = default;
    virtual Decision Judge(const IterationState&) const = 0;
};

//! \brief The four conventional tests.  Converged when max|g| <= maxGradient and rms(g) <= rmsGradient, and (at the start, or)
//! max|step| <= maxStep or |Δvalue| <= valueChange.  \c Budget once \c maxIterations steps have been evaluated unconverged.
struct GradientStepValue : public virtual StoppingRule
{
    GradientStepValue(double maxGradient, double rmsGradient, double maxStep, double valueChange, size_t maxIterations);
    Decision Judge(const IterationState&) const override;
    double maxGradient, rmsGradient, maxStep, valueChange;
    size_t maxIterations;
};

//! A proposed step.  An EMPTY \c step means nothing may move (the step space is empty): the driver stops, converged.
struct Proposal
{
    rvec_t step;
    double bound   = 0.0;     //!< the bound it was held to (0: none)
    bool   limited = false;   //!< it was cut to that bound
};

struct Verdict
{
    enum class Kind {Accept, Reject, Collapsed};
    Kind        kind = Kind::Accept;
    std::string details;      //!< why, for \c Collapsed
};

//! \brief One algorithm, ask/tell.  One instance per run: it holds the algorithm's state.
class Stepper
{
public:
    virtual ~Stepper() = default;
    virtual std::string Name() const = 0;
    //! Before any evaluation: a reason when the cost lacks a face this algorithm needs, else none.
    virtual std::optional<std::string> Refuse(const Cost&) const = 0;
    //! The run starts at \a q with \a here (gradient projected); the stepper reads the cost's capability faces here.
    virtual void Start(const Cost&, const rvec_t& q, const CostValue& here) = 0;
    //! The next trial step from the incumbent \a q, \a here.
    virtual Proposal Propose(const rvec_t& q, const CostValue& here) = 0;
    //! The cost's answer at the incumbent + \a step (gradient projected there).  \a noise bounds the value difference that is NOT a
    //! change (δ(incumbent) + δ(trial)).  Accept moves the incumbent to the trial; Reject keeps it.
    virtual Verdict Tell(const rvec_t& step, const CostValue& there, double noise) = 0;
};

//! The default value noise of a cost with no \c ValueNoise face: roundoff in the value itself, \c kRoundoffUnits x ε_mach x |value|.
//! EMPIRICAL: the rutile model's value rose by ~5 ε relative from one evaluation to the next at |g| 2e-9 (lane H, 2026-11-09);
//! 16 units bounds a short sum's accumulated rounding.
inline constexpr double kRoundoffUnits = 16.0;
double RoundoffNoise(double value);

//! \brief THE DRIVER: minimise \a cost from \a start with \a stepper until \a rule stops it.  Non-template, compiled once.
Outcome<Minimum, Failure> Minimize(Cost& cost, Stepper& stepper, const StoppingRule& rule, rvec_t start);

//! The typed client's driver: as above, with the cost's own failure handed back by value.
template <class Fail> Outcome<Minimum, MinimizeFailure<Fail>> Minimize(tCost<Fail>& cost, Stepper& stepper, const StoppingRule& rule, rvec_t start)
{
    using O=Outcome<Minimum, MinimizeFailure<Fail>>;
    auto r=Minimize(static_cast<Cost&>(cost), stepper, rule, std::move(start));
    if (r) return O::Ok(r.TakeValue());
    MinimizeFailure<Fail> f;
    static_cast<Failure&>(f)=r.Error();
    if (f.why==Failure::Why::Cost) f.cost=cost.LastFailure();
    return O::Fail(std::move(f));
}

} // namespace qchem::Optimize
