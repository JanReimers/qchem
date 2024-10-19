// File: RadialMesh.H  Abstract interface for a radial mesh.
#ifndef _RadialMesh_H_
#define _RadialMesh_H_

#include "types.H"
#include <tuple>
#include <vector>

class RadialMesh
{
public:
    virtual ~RadialMesh() {};
    
    typedef std::tuple<double,double> rw_t;
    typedef std::vector<rw_t>         vec_t;
    typedef vec_t::const_iterator const_iterator;
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
    virtual index_t        size () const=0;
    
};

inline const double& r(const RadialMesh::rw_t& rw) {return std::get<0>(rw);}
inline const double& w(const RadialMesh::rw_t& rw) {return std::get<1>(rw);}

#endif //_RadialMesh_H_

