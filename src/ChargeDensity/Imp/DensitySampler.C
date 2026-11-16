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
#include <string>
module qchem.ChargeDensity.DensitySampler;
import qchem.ChargeDensity.Internal.DensitySampler;   // the two strategies this chooses between
import qchem.BasisSet.G_FieldEvaluator;               // G_RasterTransform -- the capability that decides
import qchem.BasisSet.IrrepBasisSet;    // Evaluatable_IBS<double>: the A1 fit functions at points
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
        // A POINT SET handed in with a ball (a Becke mesh, I3): XC is evaluated there, on top of the raster strategy's exact pair -- ONLY where
        // the ball holds every pair, so its synthesis IS rho[D].  Over GPW orbitals it would be the raster's rho[D] truncated to the ball: no
        // path may project rho before XC (pin 28, the rho axis retired 2026-11-15).
        if (quad.GetMesh())
        {
            if (!(ball && ball->HoldsEveryOrbitalPair()))
                throw std::invalid_argument("MakeDensitySampler: the plane-wave ball at a point set over Gaussian (GPW) orbitals synthesizes the raster's "
                                            "rho[D] truncated to the ball, not rho[D] (doc/Pins.md pin 28) -- the delta fit on those points is that cell");
            return std::make_shared<const FitAtPointsDensitySampler>(std::make_shared<const BallFitSpace>(fb, quad), raster, quad);
        }
        return raster;
    }
    return std::make_shared<const SinglesDensitySampler>(fb, std::move(quad));
}

std::shared_ptr<const DensitySampler>
MakeAuxiliaryFitSampler(const fitbasis_t& delta, BasisSet::FitQuadrature points,
                        std::shared_ptr<const BasisSet::Orbital_1E_IBS<double>> a1)
{
    if (!points.GetMesh()) throw std::logic_error("MakeAuxiliaryFitSampler: the periodic A1 fit needs the XC point set");
    if (!delta) throw std::logic_error("MakeAuxiliaryFitSampler: the periodic A1 fit samples rho[D] through the delta basis on the XC points");
    auto fit=std::make_shared<const GaussianFitSpace>(std::move(a1), points);
    auto inner=MakeDensitySampler(delta, points);
    return std::make_shared<const FittedPotentialDensitySampler>(fit, inner);
}

std::shared_ptr<const tDensitySampler<double>>
MakeDensitySampler(const std::shared_ptr<const BasisSet::rFIT_SF_ABS>& fb, BasisSet::FitQuadrature quad)
{
    assert(fb);
    // A non-orthogonal (A1) fit's projection P rho is not rho[D]: it never reaches XC (pin 28; its projected strategy is deleted with the rho
    // axis, 2026-11-15).  The A1 cell is the fitted-potential sampler, which samples rho[D] through the delta basis.
    if (!fb->isOrtho())
        throw std::invalid_argument("MakeDensitySampler: the non-orthogonal fit "+fb->Name()+" would hand XC its projection P rho, not rho[D] "
                                    "(doc/Pins.md pin 28) -- its cell is MakeFittedPotentialSampler");
    return std::make_shared<const tSinglesDensitySampler<double>>(fb, std::move(quad));
}

std::shared_ptr<const tDensitySampler<double>>
MakeFittedPotentialSampler(const std::shared_ptr<const BasisSet::rFIT_SF_ABS>& delta, BasisSet::FitQuadrature points,
                           const std::shared_ptr<const BasisSet::rFIT_SF_ABS>& fit)
{
    if (!fit) throw std::logic_error("MakeFittedPotentialSampler: no fit basis");
    // Abstract -> abstract: the analytic metric and the values at points, the two faces the fit space asks of the set.
    const auto* metric=dynamic_cast<const BasisSet::FIT_SF_NonOrtho*>(fit.get());
    if (!metric) throw std::logic_error("MakeFittedPotentialSampler: the fit basis "+fit->Name()+" has no overlap metric (FIT_SF_NonOrtho)");
    if (!dynamic_cast<const BasisSet::Evaluatable_IBS<double>*>(fit.get()))
        throw std::logic_error("MakeFittedPotentialSampler: the fit basis "+fit->Name()+" cannot evaluate its functions at points (Evaluatable_IBS)");
    auto eval=[fit](const rvec3vec_t& R) {return rmat_t((*dynamic_cast<const BasisSet::Evaluatable_IBS<double>*>(fit.get()))(R));};
    auto space=std::make_shared<const GaussianFitSpace>(metric->Overlap(), eval, points, "the A1 Gaussian set");
    return std::make_shared<const tFittedPotentialDensitySampler<double>>(space, MakeDensitySampler(delta, std::move(points)));
}

std::shared_ptr<const tDensitySampler<double>>
MakeMultipoleFittedPotentialSampler(const std::shared_ptr<const BasisSet::rFIT_SF_ABS>& delta, BasisSet::FitQuadrature points,
                                    std::shared_ptr<const qcMesh::SiteMultipoles> sites)
{
    auto space=std::make_shared<const MultipoleFitSpace>(std::move(sites), points);
    return std::make_shared<const tFittedPotentialDensitySampler<double>>(space, MakeDensitySampler(delta, std::move(points)));
}

} //namespace
