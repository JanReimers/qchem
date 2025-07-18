// File: MoleculeMesh.C  mesh implementation

module;
export module qchem.Cluster.MoleculeMesh;
import qchem.Cluster;
import Mesh;

export class MoleculeMesh : public Mesh    
{
public:
    MoleculeMesh(                     ) {};
    MoleculeMesh(const Cluster&, const MeshParams& mp);

    virtual Mesh*  Clone() const;
};
