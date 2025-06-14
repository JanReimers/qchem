// File: AtomMesh.C  mesh implementation



#include "Imp/Cluster/AtomMesh.H"
#include <Mesh/RadialMesh.H>
#include <typeinfo>

//
//  The full mesh is just a direct product of radial and ungular meshes.
//
AtomMesh::AtomMesh(const RadialMesh& rm, const Mesh& am, const RVec3& R)
{
    for (auto arw: am)
        for (auto rrw: rm)
            push_back(r(rrw) * normalize(r(arw)),w(rrw) *  w(arw));
    ShiftOrigin(R);
}


Mesh* AtomMesh::Clone() const
{
    return new AtomMesh(*this);
}


