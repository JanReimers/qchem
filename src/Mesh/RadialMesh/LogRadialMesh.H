// File: LogRadialMesh.H  Radial Logarithmic Mesh.
#ifndef _LogRadialMesh_H_
#define _LogRadialMesh_H_



#include "Mesh/RadialMesh/RadialMeshImplementation.H"

class LogRadialMesh
    : public  virtual RadialMesh
    , public          RadialMeshImplementation
{
public:
    LogRadialMesh(double start, double stop, index_t numPoints);
};

#endif //_LogRadialMesh_H_
