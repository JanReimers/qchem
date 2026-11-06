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

std::unique_ptr<CoulombOperator> MakeRadialPoissonCoulomb(const qcMesh::Mesh& mesh, const rvec3vec_t& centres, const rvec_t& radii, int lMax)
{
    return std::make_unique<RadialPoissonCoulomb>(mesh, centres, radii, lMax);
}

std::unique_ptr<FunctionFitter_Density_NonOrtho<double>>
MakeGridDensityFitter(std::shared_ptr<const BasisSet::rFIT_SF_ABS> deltaBasis, std::unique_ptr<CoulombOperator> op)
{
    return std::make_unique<GridDensityFitter>(std::move(deltaBasis), std::move(op));
}

} //namespace
