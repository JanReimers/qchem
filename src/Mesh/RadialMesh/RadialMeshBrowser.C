// File: RadialMeshBrowser.C  Radial mesh browser



#include "Mesh/RadialMesh/RadialMeshBrowser.H"
#include "Mesh/RadialMesh/RadialMesh.H"
#include "Mesh/RadialMesh/RadialMeshImplementation.H"

RadialMeshBrowser::RadialMeshBrowser(const RadialMesh& m)
    : Rb(dynamic_cast<const RadialMeshImplementation*>(&m)->itsPoints.begin() )
    , Wb(dynamic_cast<const RadialMeshImplementation*>(&m)->itsWeights.begin())
    , RbEnd(dynamic_cast<const RadialMeshImplementation*>(&m)->itsPoints.end() )
    , WbEnd(dynamic_cast<const RadialMeshImplementation*>(&m)->itsWeights.end())
{};
