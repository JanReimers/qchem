// File: LogRadialMesh.H  Radial Logarithmic Mesh.
#ifndef _LogRadialMesh_H_
#define _LogRadialMesh_H_

#include <Mesh/RadialMesh.H>

class LogRadialMesh : public   RadialMesh
{
public:
    LogRadialMesh(double start, double stop, index_t numPoints);
};

#endif //_LogRadialMesh_H_
