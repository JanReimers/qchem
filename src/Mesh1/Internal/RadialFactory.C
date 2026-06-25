// File: Internal/RadialFactory.C  MakeRadial implementation.  This is an implementation unit of
// the qchem.Mesh1.Radial interface module, so it may reach into the Internal concretes.
module;
#include <memory>
#include <stdexcept>
module qchem.Mesh1.Radial;
import qchem.Mesh1.Radial.Internal;

namespace qcMesh1
{

std::unique_ptr<RadialMesh> MakeRadial(const MeshParams& p)
{
    switch (p.radial)
    {
    case RadialKind::MHL:
        return std::make_unique<MHLRadialMesh>(p.nRadial, p.mhl_m, p.mhl_alpha);
    case RadialKind::Log:
        return std::make_unique<LogRadialMesh>(p.logStart, p.logStop, p.nRadial);
    case RadialKind::Linear:
        return std::make_unique<LinearRadialMesh>(0.0, p.logStop, p.nRadial);
    }
    throw std::runtime_error("MakeRadial: unknown RadialKind");
}

} //namespace qcMesh1
