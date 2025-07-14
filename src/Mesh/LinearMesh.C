// File: LinearMesh.C  Linear mesh implementation.
module;
#include "oml/vector.h"
#include <cmath>
export module Mesh.LinearMesh;
export import Mesh;

export class LinearMesh : public Mesh
{
public:
    LinearMesh(double start, double stop, const RVec3& direction, index_t numPoints);
    Mesh*  Clone() const;
};


LinearMesh::LinearMesh(double start, double stop, const RVec3& direction, index_t NumPoints)
{
    RVec3 nd=normalize(direction); //Make sure its normailized.

    Vector<double> R(NumPoints);
    Vector<double> W(NumPoints);
    FillLinear(R,start,stop);
    Fill(W,1.0/NumPoints);
    
    for (auto i:R.indices())
        push_back(R(i)*nd,W(i));

//    Vector<RVec3> Rv(NumPoints);
//    Vector<double>::const_iterator b(R.begin());
//    for(Vector<RVec3> ::iterator i(Rv.begin()); i!=Rv.end(); i++,b++) *i=*b*nd;
//
//    Initialize(Rv,W);
}

Mesh* LinearMesh::Clone() const
{
    return new LinearMesh(*this);
}

