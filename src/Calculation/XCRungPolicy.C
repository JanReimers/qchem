// File: Calculation/XCRungPolicy.C  How the uniform XC feed sees the completion rung's pairs: B (pairs onto L0) or C (the rung on a multiple
// of L0), a COST choice between two exact routes, made at setup.
//
// THE CHOICE (user ruling 2026-11-26, via the controller: "I think B and C can each be valid in certain circumstances.  I am leaning
// towards supporting both options as an advance user knob." / "Auto, with B as the fallback"; lane A).  The uniform XC raster is
// \f$\rho[D]\f$ at L0's points, exactly.  The ladder's completion rung (present only below the rung gate) is a raster of its own, and
// band-truncating it down to L0 was a 9e-4 Ha XC error (the O2 box at density x 32).  Two exact routes reach L0's points
// (\c BasisSet::Gaussian::XCRungFeed, the realization; GPW_Evaluator's XC feed):
//   B  PairsOnL0     Hartree keeps its rung; XC collocates the rung's pairs a second time, directly onto L0's raster.
//   C  RungMultiple  the rung's raster is rounded up to an integer multiple m of L0's; XC reads L0's points off Hartree's own collocation.
// Both give rho[D] at L0's points to the collocation floor; they differ in COST, and C also gives Hartree a finer rung (measured: <= 6e-10 Ha).
// So the user may state one (the deck's `solid.xcRungFeed`) and a stated value WINS; \c Auto (default) takes the cheaper on the cost model
// below, and B when the model cannot be formed.  Every decision files a `policy` row with both predictions.
//
// RIGOROUS vs EMPIRICAL kept apart (CLAUDE.md): the CHOICE is structural (\c XCRungCheaperRoute: two raster shapes, no constant), the
// absolute seconds are the named calibration \c kXCRungScatterSecondsPerPoint, which decides nothing.  They meet in ONE place,
// \c DecideXCRung, and the [policy] line prints x, m, both predicted extras, kappa and the choice.
module;
#include <cmath>
#include <format>
#include <iostream>
#include <string>
export module qchem.XCRungPolicy;
export import qchem.BasisSet.Gaussian.Lattice.RasterLadder;   // RungShape (the ladder's answer), XCRungFeed (the realization)
import qchem.Policy;            // policy::Report
import qchem.Reporting;

export namespace qchem::policy
{

//! \brief What the USER states for the XC feed's view of the completion rung (the solid options' \c xcRungFeed).
enum class XCRungRoute { Auto, PairsOnL0, RungMultiple };
inline const char* NameOf(XCRungRoute r)
{
    switch (r)
    {
        case XCRungRoute::PairsOnL0:    return "pairsOnL0";
        case XCRungRoute::RungMultiple: return "rungMultiple";
        default:                        return "auto";
    }
}
inline const char* NameOf(BasisSet::Gaussian::XCRungFeed f)
{return f==BasisSet::Gaussian::XCRungFeed::RungMultiple ? "rungMultiple (C)" : "pairsOnL0 (B)";}

//! \brief THE COST MODEL, rigorous in form (no constant).  The collocation of a Gaussian pair scatters onto every raster point inside its
//! screened support; that support is a fixed region of space (its radius is set by the pair's exponents and the collocation floor), so the
//! number of points it holds -- and the scatter's cost, which is a fixed work per point -- is proportional to the raster's point DENSITY:
//! \f[ c(N)\;=\;c(N^{own})\,\prod_{i=1}^{3}\frac{N_i}{N^{own}_i}. \f]
//! Let \f$R\f$ be the rung pairs' scatter on the rung's OWN raster \f$N^{own}\f$ (Hartree pays it in either route), \f$N^0\f$ L0's raster,
//! \f$x_i=N^{own}_i/N^0_i\f$ and \f$m_i=\lceil x_i\rceil\f$ (raised to the FFT menu).  Per iteration, on top of Hartree's own collocation:
//! \f[ \Delta_B = R\prod_i\frac{N^0_i}{N^{own}_i} = \frac{R}{\prod_i x_i}, \qquad
//!     \Delta_C = R\Big(\prod_i\frac{m_iN^0_i}{N^{own}_i}-1\Big) = R\,\frac{\prod_i m_i-\prod_i x_i}{\prod_i x_i}. \f]
//! (B scatters the rung pairs once more, on L0; C scatters them once, on \f$mN^0\f$ instead of \f$N^{own}\f$; the gathers mirror these.)
//! \f$R\f$ CANCELS from the choice, which is THE RULE:
//! \f[ \Delta_C<\Delta_B\iff\prod_i m_i-\prod_i x_i<1 \f]
//! (a product of per-axis LINEAR ratios: the cost follows points per volume).  Its ISOTROPIC CASE, the rasters the ladder builds
//! (\f$x_i=x\f$, \f$m_i=m\f$), is \f$m^3-x^3<1\f$: C whenever \f$x\le1\f$ (then \f$m=1\f$ and \f$\Delta_C=R(1-x^3)/x^3<R/x^3\f$); for
//! \f$k-1<x\le k\f$ (\f$m=k\f$) only when \f$x>(k^3-1)^{1/3}\f$ -- the 1.913 threshold is \f$k=2\f$, \f$7^{1/3}\f$ (2.962 at \f$k=3\f$):
//! C wins only just below an integer multiple, B elsewhere.  The ladder rounds \f$m\f$ up to the FFT menu, which only raises C's side.
//! Not modelled: Hartree's FFT of the rung level (\f$O(N\log N)\f$; C's raised raster pays more).  MEASURED 2026-11-26 on Si Γ x16/16 (rung
//! 24 -> 30): 0.23 ms of C's 2.6 ms extra per density, ~10%.  Realness is not a factor: one scatter kernel serves every block (the phase
//! enters as c = Re[D phase]), and a k mesh collocates ONCE through the k-fold (P(R) summed over the blocks) -- Si k411 (Γ, a TRIM, complex
//! ±1/4) measured C:B 2.5 against Γ's 2.3.
//! MEASURED (the calibration's runs, below): NaF default (L0 25, rung 36, C at 50) predicts B 0.27 x (25/36)^3 = 0.090 against 0.10 measured,
//! C 0.27 x ((50/36)^3 - 1) = 0.45 against 0.48 (fractions of the base levels' scatter).
struct XCRungCost
{
    double x=0.0;          //!< \f$\prod_i x_i\f$ (the rung's own raster over L0's, in points)
    double m=0.0;          //!< \f$\prod_i m_i\f$ (C's multiple, in points)
    double relB=0.0;       //!< \f$\Delta_B/R\f$
    double relC=0.0;       //!< \f$\Delta_C/R\f$
};
inline XCRungCost RungCost(const BasisSet::Gaussian::RungShape& s)
{
    XCRungCost c;
    c.x=(double(s.Nown.x)/s.N0.x)*(double(s.Nown.y)/s.N0.y)*(double(s.Nown.z)/s.N0.z);
    c.m=double(s.m.x)*s.m.y*s.m.z;
    c.relB=1.0/c.x;
    c.relC=(c.m-c.x)/c.x;
    return c;
}
//! The structural choice: C iff \f$\prod m-\prod x<1\f$ (\c XCRungCost; the derivation there).  Ties go to B (the fallback).
inline BasisSet::Gaussian::XCRungFeed XCRungCheaperRoute(const XCRungCost& c)
{return c.m-c.x<1.0 ? BasisSet::Gaussian::XCRungFeed::RungMultiple : BasisSet::Gaussian::XCRungFeed::PairsOnL0;}

//! \brief THE MEASURED CALIBRATION kappa: the rung pairs' scatter+gather seconds per point of the rung's OWN raster, \f$R\approx\kappa\prod_i
//! N^{own}_i\f$.  It feeds ONLY the absolute seconds of the [policy] line; the choice never reads it.  Fitted 2026-11-26 (lane A,
//! ~/Code/qchem6-runs/anchors/lane-a-pin10/cost_*.log): per SCF call, the collocate + gather timed five ways back to back inside the SCF,
//! R = (all pairs) - (rung pairs dropped), paired median over the calls, on a shared box (load ~10 on 16 threads):
//!   case (L0 / rung own)      calls   base c+g/call   R/base   R (s)     Nown^3     kappa (s/pt)
//!   MnO 20/8.3 (90 / 81)        64       30.8 s        0.018   0.55     531,441     1.0e-6   (R inside the +-1% noise)
//!   NaF default (25 / 36)       74       0.040 s       0.27    0.0107    46,656     2.3e-7
//!   NaF (36 / 36)               13       0.122 s       0.12    0.0146    46,656     3.1e-7
//!   NaF (54 / 36)               14       0.427 s       0.067   0.0286    46,656     6.1e-7
//!   NaF x32/8 (72 / 50)         66       0.654 s       0.024   0.0157   125,000     1.3e-7
//! A factor ~8 scatter (the rung pairs' supports are not the cell, and the MnO row is noise-limited); the geometric mean is taken.  An
//! ORDER-OF-MAGNITUDE estimate, labelled so on the line -- refit it if the line's absolute seconds ever decide anything.
inline constexpr double kXCRungScatterSecondsPerPoint=3.5e-7;

//! The verdict before anything is reported.
struct XCRungDecision
{
    BasisSet::Gaussian::XCRungFeed feed=BasisSet::Gaussian::XCRungFeed::PairsOnL0;
    BasisSet::Gaussian::XCRungFeed autoPick=BasisSet::Gaussian::XCRungFeed::PairsOnL0;
    bool        hasRung=false;   //!< no rung: either route is the plain feed, nothing to choose
    bool        warn=false;      //!< a stated route predicted to cost more than twice the other's extra
    XCRungCost  cost;
    double      R=0.0;           //!< kappa x the rung's own points (s per iteration, estimate)
    std::string why;
};

//! \brief THE ONE PLACE the structural rule and the calibration combine.  Pure (no report): \c XCRungPolicy reports what it returns.
inline XCRungDecision DecideXCRung(XCRungRoute stated, const BasisSet::Gaussian::RungShape& s)
{
    using F=BasisSet::Gaussian::XCRungFeed;
    XCRungDecision d;
    d.hasRung=s.HasRung();
    if (!d.hasRung)   // no rung: the feed is L0 plus the coarse levels, whichever route is named
    {
        d.feed = stated==XCRungRoute::RungMultiple ? F::RungMultiple : F::PairsOnL0;
        d.why  = "the ladder has no completion rung (L0 is the finest level): every pair is collocated on its level and nothing is rounded";
        return d;
    }
    d.cost=RungCost(s);
    d.R=kXCRungScatterSecondsPerPoint*double(s.Nown.x)*s.Nown.y*s.Nown.z;
    d.autoPick=XCRungCheaperRoute(d.cost);
    const std::string rule=std::format("prod m - prod x = {:.3f} - {:.3f} = {:.3f} {} 1", d.cost.m, d.cost.x, d.cost.m-d.cost.x,
                                       d.cost.m-d.cost.x<1.0 ? "<" : ">=");
    switch (stated)
    {
        case XCRungRoute::Auto:
            d.feed=d.autoPick;
            d.why="the cheaper extra: "+rule+(d.autoPick==F::RungMultiple ? " (C: the rung on m x L0)" : " (B: the rung pairs onto L0)");
            break;
        case XCRungRoute::PairsOnL0:
        case XCRungRoute::RungMultiple:
        {
            d.feed = stated==XCRungRoute::RungMultiple ? F::RungMultiple : F::PairsOnL0;
            const double mine  = d.feed==F::RungMultiple ? d.cost.relC : d.cost.relB;
            const double other = d.feed==F::RungMultiple ? d.cost.relB : d.cost.relC;
            d.warn = mine>2.0*other;
            d.why  = std::string("stated ")+NameOf(stated)+(d.warn ? std::format(": predicted {:.2f} R against {:.2f} R for the other (auto would pick {}; {})",
                                                                                mine, other, NameOf(d.autoPick), rule)
                                                                    : std::string());
            break;
        }
    }
    return d;
}

//! \brief The policy door the solid setup calls once its ladder is built: decides, prints the [policy] line, files the row, and returns the
//! realization the ladder and the XC feed must be built for.
inline BasisSet::Gaussian::XCRungFeed XCRungPolicy(XCRungRoute stated, const BasisSet::Gaussian::RungShape& s)
{
    const XCRungDecision d=DecideXCRung(stated, s);
    report::json detail={{"stated", NameOf(stated)}, {"chosen", NameOf(d.feed)}, {"hasRung", d.hasRung},
                         {"N0", {s.N0.x, s.N0.y, s.N0.z}}};
    if (d.hasRung)
    {
        detail["Nown"]={s.Nown.x, s.Nown.y, s.Nown.z};
        detail["m"]={s.m.x, s.m.y, s.m.z};
        detail["x"]=d.cost.x; detail["prodM"]=d.cost.m;
        detail["autoPick"]=NameOf(d.autoPick);
        detail["extraB_overR"]=d.cost.relB; detail["extraC_overR"]=d.cost.relC;
        detail["R_seconds_estimate"]=d.R; detail["extraB_seconds_estimate"]=d.R*d.cost.relB; detail["extraC_seconds_estimate"]=d.R*d.cost.relC;
        detail["kappa"]=kXCRungScatterSecondsPerPoint;
        std::cout<<"[policy] XC rung feed ("<<NameOf(stated)<<"): "<<NameOf(d.feed)<<" -- "<<d.why
                 <<std::format("  [L0 {}^3, rung own {}^3 (x {:.3f} in points), m {} (prod {:.0f}); extra per iteration B {:.2f} R, C {:.2f} R, "
                               "R ~ kappa x {}^3 = {:.2g} s (kappa {:.2g} s/pt, an order-of-magnitude calibration)]",
                               s.N0.x, s.Nown.x, d.cost.x, s.m.x, d.cost.m, d.cost.relB, d.cost.relC, s.Nown.x, d.R, kXCRungScatterSecondsPerPoint)
                 <<std::endl;
        if (d.warn) std::cout<<"[policy] XC rung feed WARNING: "<<d.why<<std::endl;
    }
    Report("XCRungFeed", std::string("the uniform XC feed sees the completion rung as ")+NameOf(d.feed),
           d.warn ? Verdict::Warn : Verdict::Choose, d.why,
           d.hasRung ? std::optional<Cost>(Cost{.wallSeconds=d.R*(d.feed==BasisSet::Gaussian::XCRungFeed::RungMultiple ? d.cost.relC : d.cost.relB)})
                     : std::nullopt,
           detail);
    return d.feed;
}

//! \brief The run's XC is NOT the plane-wave fit on the raster (a Becke mesh, a delta fit): there is no raw XC feed, so the plain ladder
//! (B: no work) is kept.  A STATED route is reported as having nothing to act on; Auto files nothing (the policy did not fire).
inline BasisSet::Gaussian::XCRungFeed XCRungNotOnRaster(XCRungRoute stated)
{
    if (stated!=XCRungRoute::Auto)
        Report("XCRungFeed", std::string("stated ")+NameOf(stated)+" has nothing to act on", Verdict::Choose,
               "the run's XC is not the plane-wave fit on the GPW raster (a Becke mesh or a delta fit), so there is no raw XC feed for the "
               "completion rung to reach; the ladder is built plain", std::nullopt, report::json{{"stated", NameOf(stated)}, {"chosen", "none"}});
    return BasisSet::Gaussian::XCRungFeed::PairsOnL0;
}

} // namespace qchem::policy
