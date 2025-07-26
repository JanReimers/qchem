// File: AtomMesh.C  Atom centered mesh interface.
export module qchem.Cluster.AtomMesh;
import qchem.Mesh;
import qchem.RadialMesh;

export class AtomMesh : public Mesh
{
public:
    AtomMesh(                              ) {};
    AtomMesh(RadialMesh*, Mesh*, const RVec3& R);

    virtual Mesh*  Clone(        ) const;
};

