// File: GaussAngularMesh.H  Angular mesh.
#ifndef _GaussAngularMesh_H_
#define _GaussAngularMesh_H_



#include <Mesh/Mesh.H>

class GaussAngularMesh : public  Mesh
{
public:
    GaussAngularMesh(index_t numDir=1);
};

#endif //_GaussAngularMesh_H_
