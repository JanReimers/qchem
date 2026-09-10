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
module qchem.ChargeDensity.DensitySampler;
import qchem.ChargeDensity.Internal.DensitySampler;   // the two strategies this chooses between
import qchem.BasisSet.G_FieldEvaluator;               // G_RasterTransform -- the capability that decides

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
    if (dynamic_cast<const BasisSet::G_RasterTransform*>(fb.get()))
        return std::make_shared<const PairDensitySampler>(fb);
    return std::make_shared<const SinglesDensitySampler>(fb, std::move(quad));
}

} //namespace
