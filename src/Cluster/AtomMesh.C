// File: AtomMesh.C  mesh implementation


#include "oml/vector3d.h"
#include "AtomMesh.H"
import RadialMesh;


//
//  The full mesh is just a direct product of radial and ungular meshes.
//
AtomMesh::AtomMesh(RadialMesh* rm, Mesh* am, const RVec3& R)
{
    for (auto arw: *am)
        for (auto rrw: *rm)
            push_back(r(rrw) * normalize(r(arw)),w(rrw) *  w(arw));
    ShiftOrigin(R);
    delete am;
    delete rm;
}


Mesh* AtomMesh::Clone() const
{
    return new AtomMesh(*this);
}


