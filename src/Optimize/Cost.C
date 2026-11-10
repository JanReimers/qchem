// File: Optimize/Cost.C  The COST face every minimizer drives, and the CAPABILITY faces beside it (V-OPTIMIZE-DIP,
// doc/Records/VOptimizeDIP_Design.md §1; design approved by the user 2026-11-09).
//
// A cost is a function of its PARAMETERS q: give it q, get back a value and, when asked, a gradient.  That is all every algorithm
// needs.  What only SOME costs can answer sits beside it as a separate data-free face (ISP), reached by a cross-cast between pure faces:
//   Projected     the step space at q (rigid modes, a polar origin, a symmetry's totally symmetric part): an orthogonal projector
//   ModelHessian  a model curvature to start a quasi-Newton history from (Lindh for a PES, later)
//   ValueNoise    a bound on the value's evaluation error (an SCF's ΔE gate): a rise inside it is not a rise (§4)
// Later, each with its first client: LeastSquares (the grid tool), Retraction (GDM), Curvature (Hv, the Newton family).
//
// ⚠ THE FACES ARE EMERGENT (user, 2026-11-09): a best guess, refactored freely until two optimizers are LSP-substitutable on one cost.
// Not tier-1 hubs; no hub-interfaces.snapshot entry.
//
// THE ADAPTER (user ruling 2026-11-09, decision 1).  A calculation is NOT a cost: a cost is the calculation PLUS a parameterisation
// (a PES through a CoordinateMap, a Hamiltonian through exp(κ)).  So a client writes a small adapter class that holds both and
// implements these faces; neither the calculation nor the minimizer knows the other.
module;
#include <optional>
#include <string>
export module qchem.Optimize.Cost;
export import qchem.Types;      // rvec_t
export import qchem.Outcome;

export namespace qchem::Optimize
{

struct Want      { bool gradient = true; };                    //!< a value-only ask sets gradient = false
struct CostValue { double value = 0.0; rvec_t gradient; };      //!< \c gradient is empty iff it was not asked for
//! What the ALGORITHM sees of a failed evaluation.  The client's own typed failure rides beside it (\c tCost, \c Minimize<Fail>).
struct CostFailure { std::string details; };

//! \brief The one face every algorithm needs.  Data-free and abstract.
class Cost
{
public:
    virtual ~Cost() = default;
    virtual size_t Dimension() const = 0;
    //! Non-const: a PES advances the density it carries from point to point.
    virtual Outcome<CostValue, CostFailure> Evaluate(const rvec_t& q, const Want& want) = 0;
};

//! \brief The step space at q: \c Project is an ORTHOGONAL projector (P = Pᵀ = P²) applied as an operator, never formed by the cost.
//! Absent: every direction is free.  The minimizer applies it to every gradient it judges and every step it takes.
class Projected
{
public:
    virtual ~Projected() = default;
    virtual rvec_t Project(const rvec_t& q, const rvec_t& v) const = 0;
};

//! \brief A model Hessian at q, applied to a vector (symmetric, Ha/q² for a PES).  Absent: the algorithm's own c·1.
class ModelHessian
{
public:
    virtual ~ModelHessian() = default;
    virtual rvec_t ApplyH0(const rvec_t& q, const rvec_t& v) const = 0;
};

//! \brief A bound δ(q) on the evaluation error of the value at q (the value's own units).  Two evaluations differ by up to δ₁ + δ₂ with
//! no change in the true value: an SCF stops at its ΔE gate, a grid has its delivered error.  Absent: exact to roundoff
//! (\c RoundoffNoise).
class ValueNoise
{
public:
    virtual ~ValueNoise() = default;
    virtual double Noise(const rvec_t& q, double value) const = 0;
};

//! \brief THE TYPED ADAPTER BASE (not a face: it holds the last failure).  A client whose evaluation fails with its own type \a Fail
//! (\c PESFailure) implements \c EvaluateTyped; the algorithms, compiled once and non-template, see only \c CostFailure, and
//! \c Minimize<Fail> hands the typed failure back BY VALUE (no cast, no std::any).  \a Fail needs a \c details string.
template <class Fail> class tCost : public virtual Cost
{
public:
    virtual Outcome<CostValue, Fail> EvaluateTyped(const rvec_t& q, const Want& want) = 0;
    Outcome<CostValue, CostFailure> Evaluate(const rvec_t& q, const Want& want) final
    {
        auto r=EvaluateTyped(q, want);
        if (r) { itsFailure.reset(); return Outcome<CostValue, CostFailure>::Ok(r.TakeValue()); }
        itsFailure=r.Error();
        return Outcome<CostValue, CostFailure>::Fail(CostFailure{r.Error().details});
    }
    //! The failure of the LAST evaluation, if it failed.
    const std::optional<Fail>& LastFailure() const {return itsFailure;}
private:
    std::optional<Fail> itsFailure;
};

} // namespace qchem::Optimize
