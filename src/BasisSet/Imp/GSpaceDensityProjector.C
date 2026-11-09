// File: BasisSet/Imp/GSpaceDensityProjector.C  The ρ̃ projector's factory: where the concrete pairing lives.
module;
#include <memory>
#include <stdexcept>
#include <string>
module qchem.BasisSet.GSpaceDensityProjector;
import qchem.BasisSet.Internal.GSpaceDensityProjector;   // the concrete class
import qchem.BasisSet.G_FieldEvaluator;                  // G_RasterTransform (the face the factory gathers)

namespace qchem::BasisSet
{

std::unique_ptr<PairProjector<dcmplx>> MakeGSpaceDensityProjector(std::shared_ptr<const cFIT_SF_ABS> sf, std::string starAverageBucket)
{
    if (!sf) throw std::logic_error("MakeGSpaceDensityProjector: no fit basis");
    auto* rt=dynamic_cast<const G_RasterTransform*>(sf.get());   // the factory gathers the raster face of the same object
    if (!rt) throw std::logic_error("MakeGSpaceDensityProjector: the fit basis '"+sf->Name()+"' has no raster (G_RasterTransform)");
    return std::make_unique<GSpaceDensityProjector>(std::move(sf), *rt, std::move(starAverageBucket));
}

} //namespace
