// File: Fitting/Imp/CoulombOperator.C  The factories of the grid density route.
module;
#include <memory>
module qchem.Fitting.CoulombOperator;
import qchem.Fitting.Internal.DirectSumCoulomb;
import qchem.Fitting.Internal.GridDensityFitter;
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

} //namespace
