// File: Internal/AngularFactory.C  MakeAngular -- dispatch to the per-scheme builder functions.
module;
#include <stdexcept>
module qchem.Mesh.Angular;

namespace qchem::qcMesh
{

AngularMesh MakeAngular(const MeshParams& p)
{
    switch (p.angular)
    {
    case AngularKind::Gauss:         return GaussAngular(p.nAngular);
    case AngularKind::GaussLegendre: return GaussLegendreAngular(p.nAngular);
    case AngularKind::EulerMaclaren: return EulerMaclarenAngular(p.nAngular, p.em_m);
    }
    throw std::runtime_error("MakeAngular: unknown AngularKind");
}

} //namespace qchem::qcMesh
