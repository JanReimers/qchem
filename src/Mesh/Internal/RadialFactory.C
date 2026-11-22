// File: Internal/RadialFactory.C  MakeRadial -- dispatch to the per-scheme builder functions.
module;
#include <stdexcept>
module qchem.Mesh.Radial;

namespace qchem::qcMesh
{

RadialMesh MakeRadial(const MeshParams& p)
{
    switch (p.radial)
    {
    case RadialKind::MHL:    return MHLRadial(p.nRadial, p.mhl_m, p.mhl_alpha);
    case RadialKind::Log:    return LogRadial(p.logStart, p.logStop, p.nRadial);
    case RadialKind::Linear: return LinearRadial(0.0, p.logStop, p.nRadial);
    case RadialKind::DE2:    return DE2Radial(p.nRadial, p.de2_alpha, p.de2_rMin, p.de2_rMax);
    }
    throw std::runtime_error("MakeRadial: unknown RadialKind");
}

} //namespace qchem::qcMesh
