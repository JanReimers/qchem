// File: MoleculeMesh.H   Create a molecular mesh using Becke's fuzzy polyhedorns.
#ifndef _MoleculeMesh_H_
#define _MoleculeMesh_H_



#include "Mesh/MeshImplementation.H"

class Cluster;

class MoleculeMesh
    : public virtual Mesh
    , public MeshImplementation
{
public:
    MoleculeMesh(                     ) {};
    MoleculeMesh(const Cluster&, int m);

    virtual Mesh*  Clone() const;
private:
    void LoadFuzzyPoints  (const Atom&, const Cluster&, int);

};

#endif // _MoleculeMesh_H_
