// File: AtomMesh.H   Mesh implemented as an vector of MeshPoint's.
#ifndef _AtomMesh_H_
#define _AtomMesh_H_

#include <Mesh.H>

class AtomMesh : public Mesh
{
public:
    AtomMesh(                                     ) {};
    AtomMesh(const RadialMesh&, const Mesh&);

    virtual Mesh*  Clone(        ) const;
};

#endif // _AtomMesh_H_
