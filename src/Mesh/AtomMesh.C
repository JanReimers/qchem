// File: AtomMesh.C  mesh implementation



#include "Mesh/AtomMesh.H"
#include "Mesh/AngularMesh/AngularMesh.H"
#include "Mesh/AngularMesh/AngularMeshBrowser.H"
#include "Mesh/RadialMesh/RadialMesh.H"
#include "Mesh/RadialMesh/RadialMeshBrowser.H"
#include "Mesh/AngularMesh/AngularMeshImplementation.H"
#include "Mesh/AngularMesh/GaussAngularMesh.H"
#include <typeinfo>

//
//  The full mesh is just a direct product of radial and ungular meshes.
//
AtomMesh::AtomMesh(const RadialMesh& r, const AngularMesh& a)
    : MeshImplementation(r.GetNumPoints() * a.GetNumDirections())
{
    index_t n=r.GetNumPoints() * a.GetNumDirections();
    Vector<RVec3>  Points (n);
    Vector<double> Weights(n);

    Vector<RVec3> ::iterator ip(Points.begin());
    Vector<double>::iterator iw(Weights.begin());

    for (AngularMeshBrowser ab(a); ab; ab++)
    {
        for (RadialMeshBrowser  rb(r); rb; rb++,ip++,iw++)
        {
            *ip= rb.R() * normalize(ab.D()); //Use normalized version of the direction d.
            *iw= rb.W() *  ab.W();
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


