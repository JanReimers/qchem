// File: AtomMesh.C  Atom centered mesh interface.
export module qchem.Cluster.AtomMesh;
import qchem.Mesh;
import qchem.RadialMesh;

export class AtomMesh : public Mesh
{
public:
    AtomMesh(                              ) {};
    AtomMesh(RadialMesh*, Mesh*, const rvec3_t& R);

    virtual Mesh*  Clone(        ) const;
};

