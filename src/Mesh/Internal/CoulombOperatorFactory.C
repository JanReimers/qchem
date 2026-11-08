// File: Mesh/Internal/CoulombOperatorFactory.C  The factories of the point-set Coulomb operators (moved from qcFitting with them).
module;
#include <memory>
module qchem.Mesh.CoulombOperator;
import qchem.Mesh.Internal.DirectSumCoulomb;
import qchem.Mesh.Internal.RadialPoissonCoulomb;

namespace qchem::qcMesh
{

std::unique_ptr<CoulombOperator> MakeDirectSumCoulomb(const Mesh& mesh)
{
    return std::make_unique<DirectSumCoulomb>(mesh);
}

std::unique_ptr<CoulombOperator> MakeRadialPoissonCoulomb(std::shared_ptr<const SiteMultipoles> sites)
{
    return std::make_unique<RadialPoissonCoulomb>(std::move(sites));
}

std::unique_ptr<CoulombOperator> MakeMultipoleCoulomb(std::shared_ptr<const SiteMultipoles> sites, rvec_t metric)
{
    return std::make_unique<MultipoleCoulomb>(std::move(sites), std::move(metric));
}

} //namespace
