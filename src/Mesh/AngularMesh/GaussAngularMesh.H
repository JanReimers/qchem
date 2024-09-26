// File: GaussAngularMesh.H  Angular mesh.
#ifndef _GaussAngularMesh_H_
#define _GaussAngularMesh_H_



#include "Mesh/AngularMesh/AngularMeshImplementation.H"

class GaussAngularMesh :
    public  virtual AngularMesh,
    public          AngularMeshImplementation
{
public:
    GaussAngularMesh(index_t numDir=1);
};

#endif //_GaussAngularMesh_H_
