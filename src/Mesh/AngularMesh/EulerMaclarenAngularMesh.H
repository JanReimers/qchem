// File: EulerMaclarenAngularMesh.H  Angular mesh.
#ifndef _EulerMaclarenAngularMesh_H_
#define _EulerMaclarenAngularMesh_H_



#include "Mesh/AngularMesh/AngularMeshImplementation.H"

class EulerMaclarenAngularMesh :
    public  virtual AngularMesh,
    public          AngularMeshImplementation
{
public:
    EulerMaclarenAngularMesh(int L, int m);
};

#endif //_EulerMaclarenAngularMesh_H_
