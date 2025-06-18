// File: RadialMesh.H  Abstract interface for a radial mesh.
#ifndef _RadialMesh_H_
#define _RadialMesh_H_

#include "omlfwd.H"
#include <tuple>
#include <vector>


class RadialMesh
{
public:
    virtual ~RadialMesh() {};
    
    typedef std::tuple<double,double> rw_t;
    typedef std::vector<rw_t>         vec_t;
    typedef vec_t::const_iterator const_iterator;
    virtual const_iterator begin() const {return itsRWs.begin();}
    virtual const_iterator end  () const {return itsRWs.end  ();}
    virtual index_t        size () const {return itsRWs.size ();}
protected:
    void push_back(const double& r, const double& w) {itsRWs.push_back(std::make_tuple(r,w));}

private:
    vec_t itsRWs;    
};

inline const double& r(const RadialMesh::rw_t& rw) {return std::get<0>(rw);}
inline const double& w(const RadialMesh::rw_t& rw) {return std::get<1>(rw);}

#endif //_RadialMesh_H_

