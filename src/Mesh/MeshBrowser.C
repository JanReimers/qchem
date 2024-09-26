// File: MeshBrowser.C  Radial mesh browser



#include "Mesh/MeshBrowser.H"
#include "Mesh/Mesh.H"
#include "Mesh/MeshImplementation.H"

MeshBrowser::MeshBrowser(const Mesh& m)
    : Rb(dynamic_cast<const MeshImplementation&>(m).itsPoints .begin())
    , Wb(dynamic_cast<const MeshImplementation&>(m).itsWeights.begin())
    , RbEnd(dynamic_cast<const MeshImplementation&>(m).itsPoints .end())
    , WbEnd(dynamic_cast<const MeshImplementation&>(m).itsWeights.end())
{};
