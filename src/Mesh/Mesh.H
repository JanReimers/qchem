// File: Mesh.H  Abstract interface for all integration Meshes.
#ifndef _Mesh_H_
#define _Mesh_H_



#include "Misc/Persistent/Persistent.H"
#include "oml/imp/index_t.h"

template <class T> class Vector3D;
typedef Vector3D<double> RVec3;

class Mesh
    : public virtual Persistent
{
public:
    virtual ~Mesh () {};

    virtual void    ShiftOrigin(const RVec3&)=0;
    virtual index_t GetNumPoints     () const=0;
    virtual Mesh*   Clone            () const=0;
    static  Mesh*   Factory  (std::istream&)        ;
};

#endif //_Mesh_H_

