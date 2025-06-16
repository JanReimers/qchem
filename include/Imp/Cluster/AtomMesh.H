// File: AtomMesh.H   Mesh implemented as an vector of MeshPoint's.
#ifndef _AtomMesh_H_
#define _AtomMesh_H_

#include <Mesh/Mesh.H>

class AtomMesh : public Mesh
{
public:
    AtomMesh(                              ) {};
    AtomMesh(RadialMesh*, Mesh*, const RVec3& R);

    virtual Mesh*  Clone(        ) const;
};

#endif // _AtomMesh_H_
