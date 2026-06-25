// File: Internal/AngularFactory.C  MakeAngular implementation (impl unit of qchem.Mesh1.Angular).
module;
#include <memory>
#include <stdexcept>
module qchem.Mesh1.Angular;
import qchem.Mesh1.Angular.Internal;

namespace qcMesh1
{

std::unique_ptr<AngularMesh> MakeAngular(const MeshParams& p)
{
    switch (p.angular)
    {
    case AngularKind::Gauss:
        return std::make_unique<GaussAngularMesh>(p.nAngular);
    case AngularKind::GaussLegendre:
        return std::make_unique<GaussLegendreAngularMesh>(p.nAngular);
    case AngularKind::EulerMaclaren:
        return std::make_unique<EulerMaclarenAngularMesh>(p.nAngular, p.em_m);
    }
    throw std::runtime_error("MakeAngular: unknown AngularKind");
}

} //namespace qcMesh1
