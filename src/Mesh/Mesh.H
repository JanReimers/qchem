// File: Mesh.H  Abstract interface for all integration Meshes.
#ifndef _Mesh_H_
#define _Mesh_H_

#include "types.H"
#include "Misc/Persistent/Persistent.H"

class Mesh
    : public virtual Persistent
{
public:
    virtual ~Mesh () {};

    virtual void    ShiftOrigin(const RVec3&)=0;
    virtual size_t  GetNumPoints     () const=0;
    virtual Mesh*   Clone            () const=0;
    static  Mesh*   Factory  (std::istream&)        ;
};

#endif //_Mesh_H_

