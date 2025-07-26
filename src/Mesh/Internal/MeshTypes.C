// File: MeshTypes.C  Define some concrete classes for Mesh.
export module qchem.Mesh.Internal.Types;
export import qchem.Mesh;

export class EulerMaclarenAngularMesh : public  Mesh
{
public:
    EulerMaclarenAngularMesh(int L, int m);
};

export class GaussAngularMesh : public  Mesh
{
public:
    GaussAngularMesh(int numDir=1);
};

export class GaussLegendreAngularMesh :  public  Mesh
{
public:
    GaussLegendreAngularMesh(int L, int m);
};


export class LinearMesh : public Mesh
{
public:
    LinearMesh(double start, double stop, const RVec3& direction, int numPoints);
    Mesh*  Clone() const;
};
