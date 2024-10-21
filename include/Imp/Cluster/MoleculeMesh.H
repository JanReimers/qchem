// File: MoleculeMesh.H   Create a molecular mesh using Becke's fuzzy polyhedorns.
#ifndef _MoleculeMesh_H_
#define _MoleculeMesh_H_



#include <Mesh.H>

class Cluster;

class MoleculeMesh : public Mesh
    
{
public:
    MoleculeMesh(                     ) {};
    MoleculeMesh(const Cluster&, int m,const MeshParams& mp);

    virtual Mesh*  Clone() const;
private:
    void LoadFuzzyPoints  (const Atom&, const Cluster&, int m, const MeshParams&);

};

#endif // _MoleculeMesh_H_
