// File: GaussLegendreAngularMesh.H  Angular mesh.
#ifndef _GaussLegendreAngularMesh_H_
#define _GaussLegendreAngularMesh_H_



#include <Mesh.H>

class GaussLegendreAngularMesh :  public  Mesh
{
public:
    GaussLegendreAngularMesh(int L, int m);
};

#endif //_GaussLegendreAngularMesh_H_
