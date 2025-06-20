// File: LinearMesh.H  Linear mesh, mostly used for plotting.
#ifndef _LinearMesh_H_
#define _LinearMesh_H_

#include <Mesh/Mesh.H>

class LinearMesh : public Mesh
{
public:
    LinearMesh(double start, double stop, const RVec3& direction, index_t numPoints);
    Mesh*  Clone() const;
};

#endif //_LinearRadialMesh_H_
