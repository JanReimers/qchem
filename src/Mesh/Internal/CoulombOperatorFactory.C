// File: Mesh/Internal/CoulombOperatorFactory.C  The factories of the point-set Coulomb operators (moved from qcFitting with them).
module;
#include <memory>
#include <stdexcept>
module qchem.Mesh.CoulombOperator;
import qchem.Mesh.Internal.DirectSumCoulomb;
import qchem.Mesh.Internal.RadialPoissonCoulomb;

namespace qchem::qcMesh
{

rvec3vec_t CoulombOperator::MakeSelfRepulsionPullback(const rvec_t&) const
{
    throw std::logic_error("CoulombOperator: this operator gives no self-repulsion pullback (ask CanPullbackSelfRepulsion first)");
}

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
