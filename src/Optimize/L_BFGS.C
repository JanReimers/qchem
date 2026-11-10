// File: Optimize/L_BFGS.C  Limited-memory BFGS + a trust radius, as an ask/tell Stepper (V-OPTIMIZE-DIP §8 step 3: the second
// optimizer that must SUBSTITUTE for BFGS-RFO on the same cost, the LSP proof that hardens the emergent faces).
//
// THE METHOD.  The two-loop recursion (Nocedal & Wright, Alg. 7.4) applies the inverse-Hessian estimate of the last \c memory
// (s, y) pairs to the gradient, H₀ = γ·1 with γ = sᵀy / yᵀy from the newest pair (1 / initialCurvature before the first); the step is
// d = −H g, held to the trust radius.  No line search, like BFGS-RFO: the shared \c TrustRegion and its noise rules globalise it.
// O(memory · n) per step and NO dense matrix: the cost's step space enters only through \c Project applied as an operator, which is
// what the grid tool's ~1e5 parameters need.  A \c ModelHessian is not used (its inverse would be needed); the algorithm works without.
module;
#include <deque>
#include <optional>
#include <string>
export module qchem.Optimize.L_BFGS;
export import qchem.Optimize.TrustRegion;

export namespace qchem::Optimize
{

class L_BFGS : public virtual Stepper
{
public:
    struct Options
    {
        TrustOptions trust;
        double       initialCurvature = 0.5;   //!< H₀ = 1 / this until the first curvature pair
        size_t       memory           = 8;     //!< (s, y) pairs kept
    };
    //! THROWS on inconsistent options (the trust radii: \c TrustRegion; initialCurvature > 0; memory > 0).
    explicit L_BFGS(const Options& options);
    L_BFGS();   //!< the default options
    std::string                Name() const override {return "L-BFGS";}
    std::optional<std::string> Refuse(const Cost&) const override {return std::nullopt;}   // needs nothing beyond Cost
    void     Start  (const Cost&, const rvec_t& q, const CostValue& here) override;
    Proposal Propose(const rvec_t& q, const CostValue& here) override;
    Verdict  Tell   (const rvec_t& step, const CostValue& there, double noise) override;
private:
    rvec_t InverseHessianTimes(const rvec_t& g) const;   //!< the two-loop recursion
    Options             itsOpts;
    const Projected*    itsProjected = nullptr;   //!< the cost's step space, if it has one (non-owning: the run's cost outlives the run)
    TrustRegion         itsTrust;
    std::deque<rvec_t>  itsS, itsY;               //!< the curvature pairs, newest last
    rvec_t              itsQ, itsG;               //!< the incumbent and its projected gradient
    double              itsValue = 0.0;
    double              itsPredicted = 0.0;
    bool                itsLimited = false;
};

} // namespace qchem::Optimize
