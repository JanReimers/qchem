// File: RadialMesh.H  Abstract interface for a radial mesh.
#ifndef _RadialMesh_H_
#define _RadialMesh_H_



#include "oml/imp/index_t.h"

class RadialMesh
{
public:
    virtual ~RadialMesh() {};
    virtual index_t GetNumPoints() const=0;
};

#endif //_RadialMesh_H_

