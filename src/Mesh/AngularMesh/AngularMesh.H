//  file AngularMesh.H  Abstract interface for all angular integration meshes.
#ifndef _AngularMesh_H_
#define _AngularMesh_H_



#include "oml/imp/index_t.h"

template <class T> class Vector3D;
typedef Vector3D<double> RVec3;

class AngularMesh
{
public:
    virtual ~AngularMesh() {};

    virtual index_t GetNumDirections() const=0;
};

#endif //_AngularMesh_H_

