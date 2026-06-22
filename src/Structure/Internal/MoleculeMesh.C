// File: MoleculeMesh.C  mesh implementation

module;
export module qchem.Structure.MoleculeMesh;
import qchem.Structure;
import qchem.Mesh;

export class MoleculeMesh : public Mesh    
{
public:
    MoleculeMesh(                     ) {};
    MoleculeMesh(const Structure&, const MeshParams& mp);

    virtual Mesh*  Clone() const;
};
