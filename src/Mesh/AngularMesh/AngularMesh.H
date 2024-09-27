//  file AngularMesh.H  Abstract interface for all angular integration meshes.
#ifndef _AngularMesh_H_
#define _AngularMesh_H_

#include "types.H"

class AngularMesh
{
public:
    virtual ~AngularMesh() {};

    virtual size_t GetNumDirections() const=0;
};

#endif //_AngularMesh_H_

