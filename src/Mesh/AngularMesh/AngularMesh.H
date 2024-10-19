//  file AngularMesh.H  Abstract interface for all angular integration meshes.
#ifndef _AngularMesh_H_
#define _AngularMesh_H_

#include "types.H"
#include "oml/vector3d.h"
#include <tuple>
#include <vector>

class AngularMesh
{
public:
    virtual ~AngularMesh() {};

    virtual size_t GetNumDirections() const=0;

    typedef std::tuple<RVec3,double> rw_t;
    typedef std::vector<rw_t>         vec_t;
    typedef vec_t::const_iterator const_iterator;
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
    virtual index_t        size () const=0;
    
};

//inline const RVec3&  r(const AngularMesh::rw_t& rw) {return std::get<0>(rw);}
//inline const double& w(const AngularMesh::rw_t& rw) {return std::get<1>(rw);}

#endif //_AngularMesh_H_

