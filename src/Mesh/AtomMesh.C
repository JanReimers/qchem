// File: AtomMesh.C  mesh implementation



#include "Mesh/AtomMesh.H"
#include "Mesh/AngularMesh/AngularMesh.H"
#include "Mesh/RadialMesh/RadialMesh.H"
#include "Mesh/AngularMesh/AngularMeshImplementation.H"
#include "Mesh/AngularMesh/GaussAngularMesh.H"
#include <typeinfo>

//
//  The full mesh is just a direct product of radial and ungular meshes.
//
AtomMesh::AtomMesh(const RadialMesh& rm, const AngularMesh& am)
    : MeshImplementation(rm.size() * am.size())
{
    index_t n=rm.size() * am.size();
    Vector<RVec3>  Points (n);
    Vector<double> Weights(n);

    Vector<RVec3> ::iterator ip(Points.begin());
    Vector<double>::iterator iw(Weights.begin());

    for (auto arw: am)
    {
        for (auto rrw: rm)
        {
            *ip= r(rrw) * normalize(r(arw)); //Use normalized version of the direction d.
            *iw= w(rrw) *  w(arw);
            iw++;
            ip++;
           // itsRWs.push_back(std::make_tuple(rb.R() * normalize(ab.D()),rb.W() *  ab.W()));
        }
    }
    //for (auto i:R.indices()) itsRWs.push_back(std::make_tuple(R(i),W(i)));
    Initialize(Points,Weights);
}


Mesh* AtomMesh::Clone() const
{
    return new AtomMesh(*this);
}


