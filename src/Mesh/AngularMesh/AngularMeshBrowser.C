// File: AngularMeshBrowser.C  Radial angular mesh browser



#include "Mesh/AngularMesh/AngularMeshBrowser.H"
#include "Mesh/AngularMesh/AngularMesh.H"
#include "Mesh/AngularMesh/AngularMeshImplementation.H"

AngularMeshBrowser::AngularMeshBrowser(const AngularMesh& m)
    : Db(dynamic_cast<const AngularMeshImplementation*>(&m)->itsDirections.begin())
    , Wb(dynamic_cast<const AngularMeshImplementation*>(&m)->itsWeights.begin())
    , DbEnd(dynamic_cast<const AngularMeshImplementation*>(&m)->itsDirections.end())
    , WbEnd(dynamic_cast<const AngularMeshImplementation*>(&m)->itsWeights.end())
{};
