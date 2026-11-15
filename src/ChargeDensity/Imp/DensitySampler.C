// File: ChargeDensity/Imp/DensitySampler.C  The factory -- the ONE place that names both concrete strategies.
//
// ★ THIS UNIT IS WHY THE CONCRETES CAN BE INTERNAL.  It belongs to the PUBLIC module
// `qchem.ChargeDensity.DensitySampler` (so `MakeDensitySampler` is declared and defined in one module, as
// module linkage requires) while being free to import the Internal one -- which is the project's standing
// "abstract interface + factory public, concrete impl internal" shape.  A client above gets the face and
// the factory and never sees a strategy type.
module;
#include <cassert>
#include <memory>
#include <stdexcept>
module qchem.ChargeDensity.DensitySampler;
import qchem.ChargeDensity.Internal.DensitySampler;   // the two strategies this chooses between
import qchem.BasisSet.G_FieldEvaluator;               // G_RasterTransform -- the capability that decides
import qchem.ScalarFunction;                          // RasterPointCapture: the raster hands its points to one bulk call
import qchem.Types;

namespace qchem::ChargeDensity
{

// The grid-charge toggle lives WITH its declaration (module linkage): declared in the public interface
// unit, defined here.  It used to sit in the Pair strategy's impl, which is where it is USED -- but that
// unit now belongs to the Internal module, and a definition cannot cross modules.
// The grid-charge diagnostic's process-wide toggle.  It lives with the PAIR route because that is the only
// route that collocates rho onto a raster and can therefore report what the raster lost.
// ⚠ Process-wide MUTABLE state in a library -- flagged in doc/CleanupCandidates.md R1.0e as something the
// run report should own instead (theRunPolicy() already carries every other run-scoped switch).
bool& ReportGridCharge() { static bool on = false; return on; }

// CAPABILITY DECIDES (doc/OpenWork.md): a delta basis carries points and nothing else -> singles; a
// raster-backed one carries the FFT transforms and keys the 3-centre tensor -> pair.  One decision, taken
// once, and latched for the run by the simple fact that the Hamiltonian builds this object once.
std::shared_ptr<const DensitySampler>
MakeDensitySampler(const fitbasis_t& fb, BasisSet::FitQuadrature quad)
{
    assert(fb);
    // CAPABILITY decides.  A basis that carries the {r}<->{G} transforms is raster-backed, so its
    // collocation pair (Overlap3C's applyRaw/applyRawAdjoint) exists and the PAIR route is available --
    // and preferred, being the production GPW path.  Anything else can only be contracted through a Phi
    // table: SINGLES.  Note this asks what the basis CAN do, never what it IS.
    // A ball that holds every orbital pair (a plane-wave orbital basis at >= 4 Ecut) makes P rho = rho[D] exactly: the PROJECTED
    // strategy is the exact route there, and the pair route would have nothing to collocate (a plane-wave block has no raw pair).
    // Any other raster ball (GPW) can only reach rho[D] by collocation: PAIR.  (V-XC-ONE-TERM step 2, 2026-11-09: the ball used
    // to be the pair route's per-call fallback.)
    if (dynamic_cast<const BasisSet::G_RasterTransform*>(fb.get()))
    {
        auto* ball=dynamic_cast<const BasisSet::G_OrbitalPairBall*>(fb.get());   // abstract -> abstract
        std::shared_ptr<const DensitySampler> raster;
        if (ball && ball->HoldsEveryOrbitalPair()) raster=std::make_shared<const ProjectedDensitySampler>(fb);
        else                                       raster=std::make_shared<const PairDensitySampler>(fb);
        // A POINT SET handed in with a ball (a Becke mesh, I3): XC is evaluated there, on top of the raster strategy's exact pair.
        if (quad.GetMesh()) return std::make_shared<const FitAtPointsDensitySampler>(std::make_shared<const BallFitSpace>(fb, quad), raster, quad);
        return raster;
    }
    return std::make_shared<const SinglesDensitySampler>(fb, std::move(quad));
}

// The raster's points, in raster order: the one bulk call G_RasterTransform::Sample makes, captured.
namespace
{
class RasterPointCapture : public virtual ScalarFunction<double>
{
public:
    double  operator()(const rvec3_t&) const override {throw std::logic_error("RasterPointCapture: bulk only");}
    rvec3_t Gradient  (const rvec3_t&) const override {throw std::logic_error("RasterPointCapture: no gradient");}
    rvec_t  operator()(const rvec3vec_t& rs) const override {itsPoints=rs; return rvec_t(rs.size(), 0.0);}
    mutable rvec3vec_t itsPoints;
};
}

std::shared_ptr<const DensitySampler>
MakeAuxiliaryFitSampler(const fitbasis_t& ball, const fitbasis_t& delta, BasisSet::FitQuadrature points,
                        std::shared_ptr<const BasisSet::Orbital_1E_IBS<double>> a1, bool projected)
{
    if (!points.GetMesh()) throw std::logic_error("MakeAuxiliaryFitSampler: the periodic A1 fit needs the XC point set");
    if (!projected)
    {
        auto fit=std::make_shared<const GaussianFitSpace>(std::move(a1), rvec3vec_t{}, 0.0, points);
        return std::make_shared<const FittedPotentialDensitySampler>(fit, MakeDensitySampler(delta, points));
    }
    const auto* raster=dynamic_cast<const BasisSet::G_RasterTransform*>(ball.get());   // abstract -> abstract
    if (!raster) throw std::logic_error("MakeAuxiliaryFitSampler: the projected cell needs the plane-wave fit's raster");
    RasterPointCapture capture;
    (void)raster->Sample(capture);
    const double weight=raster->Integral(rvec_t(raster->RasterSize(), 1.0))/double(raster->RasterSize());   // Omega/N
    auto fit=std::make_shared<const GaussianFitSpace>(std::move(a1), std::move(capture.itsPoints), weight, points);
    return std::make_shared<const FitAtPointsDensitySampler>(fit, std::make_shared<const PairDensitySampler>(ball), points);
}

std::shared_ptr<const tDensitySampler<double>>
MakeDensitySampler(const std::shared_ptr<const BasisSet::rFIT_SF_ABS>& fb, BasisSet::FitQuadrature quad)
{
    assert(fb);
    // A non-orthogonal (A1) fit is a PROJECTION P onto its span, and the projected strategy evaluates P rho on the XC mesh handed in
    // (V-XC-ONE-TERM step 3).  The caller reaches it only for xc.rho "fit": the "exact" A1 cell is the fitted-v term, not a sampler.
    if (!fb->isOrtho()) return std::make_shared<const OverlapProjectedDensitySampler>(fb, std::move(quad));
    return std::make_shared<const tSinglesDensitySampler<double>>(fb, std::move(quad));
}

} //namespace
