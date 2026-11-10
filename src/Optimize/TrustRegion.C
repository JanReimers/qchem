// File: Optimize/TrustRegion.C  The trust radius every no-line-search stepper shares, with the noise rules of
// doc/Records/VOptimizeDIP_Design.md §4 in ONE place (so BFGS-RFO and L-BFGS judge a rise identically: the LSP pair of §8 step 3).
//
//   RULE 1: only a rise beyond the noise (δ(here) + δ(there)) is a rise.  It rejects the step and cuts the radius x 1/4; a radius cut
//           below its minimum COLLAPSES the run.
//   RULE 2: the ratio actual / predicted means something only when the prediction exceeds the noise.  Above it: grow x2 on a
//           well-predicted step the radius held back (ratio > 3/4), shrink x1/2 on a poor one (ratio < 1/4).  Below it: HOLD.
// The constants 1/4, 3/4, x2, x1/2, x1/4 are the textbook trust-region schedule (Nocedal & Wright, Alg. 4.1), carried over unchanged
// from GeometryOptimizer (lane B, 2026-11-07).
module;
#include <string>
export module qchem.Optimize.TrustRegion;
export import qchem.Optimize.Minimize;   // Verdict

export namespace qchem::Optimize
{

struct TrustOptions
{
    double initial = 0.3;
    double max     = 1.0;
    double min     = 1.0e-6;   //!< a rejection that cuts the radius below this ends the run: Verdict Collapsed
};

class TrustRegion
{
public:
    //! THROWS unless 0 < min and 0 < initial <= max.
    explicit TrustRegion(const TrustOptions& options);
    void   Reset() {itsRadius=itsOpts.initial;}
    double Radius() const {return itsRadius;}
    //! Judge a trial (rules 1 and 2) and update the radius.  \a change = value(there) − value(here); \a predicted = the model's change;
    //! \a limited = the step was cut to the radius.
    Verdict Judge(double change, double predicted, double noise, bool limited);
private:
    TrustOptions itsOpts;
    double       itsRadius;
};

} // namespace qchem::Optimize
