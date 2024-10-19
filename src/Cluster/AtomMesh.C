// File: AtomMesh.C  mesh implementation



#include "Imp/Cluster/AtomMesh.H"
#include <RadialMesh.H>
#include "Imp/Mesh/GaussAngularMesh.H"
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


