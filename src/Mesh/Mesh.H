// File: Mesh.H  Abstract interface for all integration Meshes.
#ifndef _Mesh_H_
#define _Mesh_H_

#include "types.H"
#include "Misc/Persistent/Persistent.H"
#include "oml/vector3d.h"
#include <tuple>
#include <vector>

class Mesh
    : public virtual Persistent
{
public:
    virtual ~Mesh () {};

    virtual void    ShiftOrigin(const RVec3&)=0;
    virtual size_t  GetNumPoints     () const=0;
    virtual Mesh*   Clone            () const=0;
    static  Mesh*   Factory  (std::istream&);
    
    typedef std::tuple<RVec3,double> rw_t;
    typedef std::vector<rw_t>         vec_t;
    typedef vec_t::const_iterator const_iterator;
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
    virtual index_t        size () const=0;
};

inline const RVec3 & r(const Mesh::rw_t& rw) {return std::get<0>(rw);}
inline const double& w(const Mesh::rw_t& rw) {return std::get<1>(rw);}


#endif //_Mesh_H_

