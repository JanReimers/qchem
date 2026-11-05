// File: LASolver/OverlapConditioning.C  How well conditioned one overlap matrix S is: the ONE analysis every facade's basis
// report rests on (doc/OpenWork.md "The basis conditioning report", doc/Pins.md pin 7).
//
// Pure linear algebra: it knows no basis type and edits nothing.  Report only; a trim is a separate, stated decision.
module;
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>
export module qchem.LASolver.OverlapConditioning;
export import qchem.Types;
import qchem.Blaze;
import qchem.LASolver;

export namespace qchem
{

//! \brief The pivot floor a basis report counts near-dependent AOs at: the AOs \c Ortho::CholeskyPivoted WOULD drop at this
//! \c orthoTol.  A WARNING level and an empirical calibration, not a safety line and not a derived bound.
//!
//! Measured (2026-11-04): Li bcc q3 k12 with a diffuse s 0.04 had a min kept pivot of 1.26e-7 and collapsed to -1418 Ha per
//! atom; beta-MnO2 at lambda_min 6.75e-7 ran green.  1e-7 would miss the Li run; 1e-6 flags both, so no single floor separates
//! healthy from broken, and a flagged run is a run to look at, not a run that is wrong.
inline constexpr double ConditioningFloor=1e-6;

//! What one overlap matrix says about its basis.
struct OverlapConditioning
{
    double lambdaMin=0.0;   //!< Smallest eigenvalue of S (negative when roundoff has made S indefinite).
    double lambdaMax=0.0;   //!< Largest eigenvalue of S.
    double cond     =0.0;   //!< \f$\lambda_{\max}/\min_i|\lambda_i|\f$; \c numeric_limits::max() for an exactly singular S.
    size_t wouldDrop=0;     //!< AOs a pivoted Cholesky at \a floor would drop.
};

//! Analyse the Hermitian overlap \a S at pivot floor \a floor.  An empty S answers all zeros.
template <class T> OverlapConditioning AnalyzeOverlap(const hmat_t<T>& S, double floor=ConditioningFloor)
{
    if (S.rows()==0) return {};
    rvec_t d; mat_t<T> U;
    blazem::eigen(S, d, U);   // ascending
    double minAbs=std::fabs(d[0]);
    for (double v : d) minAbs=std::min(minAbs, std::fabs(v));
    const double lmax=d[d.size()-1];
    return {.lambdaMin=d[0], .lambdaMax=lmax, .cond=minAbs>0.0 ? lmax/minAbs : std::numeric_limits<double>::max(),
            .wouldDrop=PivotedCholeskyDrops<T>(S, floor).size()};
}

} // namespace qchem
