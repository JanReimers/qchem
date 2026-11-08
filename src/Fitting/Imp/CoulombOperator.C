// File: Fitting/Imp/CoulombOperator.C  The factories of the grid density route.
module;
#include <memory>
#include <stdexcept>
#include <string>
module qchem.Fitting.CoulombOperator;
import qchem.Fitting.Internal.DirectSumCoulomb;
import qchem.Fitting.Internal.GridDensityFitter;
import qchem.Fitting.Internal.GSpaceDensityFitter;
import qchem.BasisSet.G_FieldEvaluator;   // the plane-wave fit basis faces the factory gathers
import qchem.Fitting.Internal.RadialPoissonCoulomb;

namespace qchem::Fitting
{

std::unique_ptr<CoulombOperator> MakeDirectSumCoulomb(const qcMesh::Mesh& mesh)
{
    return std::make_unique<DirectSumCoulomb>(mesh);
}

std::unique_ptr<CoulombOperator> MakeRadialPoissonCoulomb(std::shared_ptr<const qcMesh::SiteMultipoles> sites)
{
    return std::make_unique<RadialPoissonCoulomb>(std::move(sites));
}

std::unique_ptr<CoulombOperator> MakeMultipoleCoulomb(std::shared_ptr<const qcMesh::SiteMultipoles> sites, rvec_t metric)
{
    return std::make_unique<MultipoleCoulomb>(std::move(sites), std::move(metric));
}

std::unique_ptr<GridDensityFit>
MakeGridDensityFitter(std::shared_ptr<const BasisSet::rFIT_SF_ABS> deltaBasis, std::unique_ptr<CoulombOperator> op)
{
    return std::make_unique<GridDensityFitter>(std::move(deltaBasis), std::move(op));
}

std::unique_ptr<FirstOrderDensityFit<dcmplx>> MakeGSpaceDensityFitter(std::shared_ptr<const BasisSet::cFIT_CD_ABS> pw)
{
    // The factory gathers the plane-wave fit basis's faces (the concrete association lives here, not in the class).
    if (!pw) throw std::logic_error("MakeGSpaceDensityFitter: no fit basis");
    auto* sf=dynamic_cast<const BasisSet::cFIT_SF_ABS*>(pw.get());
    auto* pk=dynamic_cast<const BasisSet::G_PoissonKernel*>(pw.get());
    auto* rt=dynamic_cast<const BasisSet::G_RasterTransform*>(pw.get());
    if (!sf || !pk || !rt)
        throw std::logic_error("MakeGSpaceDensityFitter: the fit basis '"+pw->Name()+"' lacks a G-space face (overlap, Poisson kernel, raster)");
    return std::make_unique<GSpaceDensityFitter>(std::move(pw), *sf, *pk, *rt);
}

} //namespace
