// File: AtomMesh.C  mesh implementation



#include "Mesh/AtomMesh.H"
#include "Mesh/RadialMesh/RadialMesh.H"
#include "Mesh/AngularMesh/GaussAngularMesh.H"
#include <typeinfo>

//
//  The full mesh is just a direct product of radial and ungular meshes.
//
AtomMesh::AtomMesh(const RadialMesh& rm, const Mesh& am)
{
    for (auto arw: am)
        for (auto rrw: rm)
            push_back(r(rrw) * normalize(r(arw)),w(rrw) *  w(arw));
        
}


Mesh* AtomMesh::Clone() const
{
    return new AtomMesh(*this);
}


