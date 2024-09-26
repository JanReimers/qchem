// File: GaussLegendreAngularMesh.H  Angular mesh.
#ifndef _GaussLegendreAngularMesh_H_
#define _GaussLegendreAngularMesh_H_



#include "Mesh/AngularMesh/AngularMeshImplementation.H"

class GaussLegendreAngularMesh :
    public  virtual AngularMesh,
    public          AngularMeshImplementation
{
public:
    GaussLegendreAngularMesh(int L, int m);
};

#endif //_GaussLegendreAngularMesh_H_
