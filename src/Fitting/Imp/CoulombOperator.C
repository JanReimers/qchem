// File: Fitting/Imp/CoulombOperator.C  The factories of the grid density route.
module;
#include <memory>
#include <stdexcept>
#include <string>
module qchem.Fitting.CoulombOperator;
import qchem.Fitting.Internal.GridDensityFitter;
import qchem.Fitting.Internal.GSpaceDensityFitter;
import qchem.Fitting.Internal.FunctionFitterImp;   // IntegralConstrainedFF (the Gaussian Dunlap fit)
import qchem.BasisSet.G_FieldEvaluator;
import qchem.BasisSet.CoulombDeltaFit_IBS;   // PointCoulombIntegrals   // the plane-wave fit basis faces the factory gathers

namespace qchem::Fitting
{

std::unique_ptr<GridDensityFit>
MakeGridDensityFitter(std::shared_ptr<const BasisSet::rFIT_SF_ABS> deltaBasis)
{
    // The factory gathers the basis's Coulomb-integral face (the concrete association lives here, not in the class).
    if (!deltaBasis) throw std::logic_error("MakeGridDensityFitter: no fit basis");
    auto* coulomb=dynamic_cast<const BasisSet::PointCoulombIntegrals*>(deltaBasis.get());
    if (!coulomb)
        throw std::logic_error("MakeGridDensityFitter: the fit basis '"+deltaBasis->Name()+"' carries no Coulomb integrals "
                               "(wrap a δ basis in BasisSet::CoulombDeltaFit_IBS with its operator)");
    return std::make_unique<GridDensityFitter>(std::move(deltaBasis), *coulomb);
}

std::unique_ptr<DunlapDensityFit<double>> MakeDunlapDensityFitter(std::shared_ptr<const BasisSet::rFIT_CD_ABS> fbs)
{
    if (!fbs) throw std::logic_error("MakeDunlapDensityFitter: no fit basis");
    auto nonOrtho=std::dynamic_pointer_cast<const BasisSet::FIT_CD_NonOrtho>(fbs);   // its Coulomb metric-solve face
    if (fbs->isOrtho() || !nonOrtho)
        throw std::logic_error("MakeDunlapDensityFitter: the fit basis '"+fbs->Name()+"' has no Coulomb metric to solve in");
    return std::make_unique<IntegralConstrainedFF<double>>(nonOrtho);
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
