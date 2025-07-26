// File: RadialMeshTypes.C  Define some concrete radial mesh classes.
export module qchem.Mesh.Internal.RadialTypes;
export import qchem.RadialMesh;

//
// C.W. Murray,N.C. Handy and G.J. Laming, Mol. Phys.78 (1993) 997. 
//
export class MHLRadialMesh : public  RadialMesh
{
public:
    MHLRadialMesh(int NumPoints, int m, double alpha);
};

export class LogRadialMesh : public   RadialMesh
{
public:
    LogRadialMesh(double start, double stop, int numPoints);
};
