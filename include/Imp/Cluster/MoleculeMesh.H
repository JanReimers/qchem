// File: MoleculeMesh.H   Create a molecular mesh using Becke's fuzzy polyhedorns.
#ifndef _MoleculeMesh_H_
#define _MoleculeMesh_H_


#include <Mesh/fwd.H>

import Mesh;

class Cluster;

class MoleculeMesh : public Mesh
    
{
public:
    MoleculeMesh(                     ) {};
    MoleculeMesh(const Cluster&, const MeshParams& mp);

    virtual Mesh*  Clone() const;
};

#endif // _MoleculeMesh_H_
