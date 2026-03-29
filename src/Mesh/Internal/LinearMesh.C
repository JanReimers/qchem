// File: LinearMesh.C  Linear mesh implementation.
module;
#include <cmath>
// #include <blaze/math/expressions/DVecGenExpr.h>
#include <blaze/Math.h>
module qchem.Mesh.Internal.Types;

LinearMesh::LinearMesh(double start, double stop, const rvec3_t& direction, int NumPoints)
{
    rvec3_t nd=normalize(direction); //Make sure its normailized.
    rvec_t R=blaze::linspace(NumPoints,start,stop);
    
    for (auto r:R)
        push_back(r*nd,1.0/NumPoints);

}

Mesh* LinearMesh::Clone() const
{
    return new LinearMesh(*this);
}

