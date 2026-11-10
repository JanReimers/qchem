// File: Optimize/BFGS_RFO.C  BFGS + a rational-function (RFO) step + a trust radius, as an ask/tell Stepper.  MOVED from
// qchem::GeometryOptimizer (lane B, 2026-11-07) with the three noise rules of doc/Records/VOptimizeDIP_Design.md §4 (lane H, 2026-11-09).
//
// THE METHOD.  BFGS on the Hessian in q, initialised from the cost's \c ModelHessian when it has one; an RFO step in the cost's step
// space (\c Projected), which is a descent direction for any Hessian; a trust radius that grows on a well-predicted step and shrinks on
// a poor one.  No line search: an evaluation can be a whole SCF, and the radius is the cheaper safeguard.  The step space is found by
// diagonalising the projector, O(n³) per step: fine for a geometry (n up to ~10³), not for ~10⁵ parameters (L-BFGS).
//
// THE NOISE RULES (§4; the rutile model collapsed its radius at a converged point because each retry rose by roundoff): rules 1 and 2
// are the shared \c TrustRegion's (qchem.Optimize.TrustRegion), rule 3 the driver's.
module;
#include <optional>
#include <string>
export module qchem.Optimize.BFGS_RFO;
export import qchem.Optimize.TrustRegion;

export namespace qchem::Optimize
{

class BFGS_RFO : public virtual Stepper
{
public:
    struct Options
    {
        TrustOptions trust;
        double       initialCurvature = 0.5;   //!< H₀ = this x 1 when the cost offers no ModelHessian
    };
    //! THROWS on inconsistent options (the trust radii: \c TrustRegion; initialCurvature > 0).
    explicit BFGS_RFO(const Options& options);
    BFGS_RFO();   //!< the default options (a separate constructor: a default argument cannot use the nested struct's initialisers here)
    std::string                Name() const override {return "BFGS-RFO";}
    std::optional<std::string> Refuse(const Cost&) const override {return std::nullopt;}   // needs nothing beyond Cost
    void     Start  (const Cost&, const rvec_t& q, const CostValue& here) override;
    Proposal Propose(const rvec_t& q, const CostValue& here) override;
    Verdict  Tell   (const rvec_t& step, const CostValue& there, double noise) override;
private:
    Options             itsOpts;
    const Projected*    itsProjected = nullptr;   //!< the cost's step space, if it has one (non-owning: the run's cost outlives the run)
    rmat_t              itsH;
    TrustRegion         itsTrust;
    rvec_t              itsQ, itsG;               //!< the incumbent and its projected gradient
    double              itsValue = 0.0;
    double              itsPredicted = 0.0;       //!< the model's change for the last proposal
    bool                itsLimited = false;
};

} // namespace qchem::Optimize
