// File: MHLRadialMesh.H  Murry, Handy Laming radial mesh.
#ifndef _MHLRadialMesh_H_
#define _MHLRadialMesh_H_



#include "Mesh/RadialMesh/RadialMeshImplementation.H"

class MHLRadialMesh
    : public  virtual RadialMesh
    , public          RadialMeshImplementation
{
public:
    MHLRadialMesh(index_t NumPoints, int m, double alpha);
};

#endif //_MHLRadialMesh_H_
