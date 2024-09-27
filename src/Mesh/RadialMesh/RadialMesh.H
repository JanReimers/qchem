// File: RadialMesh.H  Abstract interface for a radial mesh.
#ifndef _RadialMesh_H_
#define _RadialMesh_H_

#include "types.H"

class RadialMesh
{
public:
    virtual ~RadialMesh() {};
    virtual index_t GetNumPoints() const=0;
};

#endif //_RadialMesh_H_

