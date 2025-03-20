// File: Mesh.H  Abstract interface for all integration Meshes.
#ifndef _Mesh_H_
#define _Mesh_H_

#include "types.H"
#include "oml/vector3d.h"
#include <tuple>
#include <vector>

class Mesh 
{ 
public:
    virtual ~Mesh () {};

    virtual void    ShiftOrigin(const RVec3&);
    virtual Mesh*   Clone      () const;
    
    typedef std::tuple<RVec3,double> rw_t;
    typedef std::vector<rw_t>         vec_t;
    typedef vec_t::const_iterator const_iterator;
    virtual const_iterator begin() const {return itsRWs.begin();}
    virtual const_iterator end  () const {return itsRWs.end  ();}
    virtual index_t        size () const {return itsRWs.size ();}

protected:
    void push_back(const RVec3& r, const double& w) {itsRWs.push_back(std::make_tuple(r,w));}
private:
    vec_t itsRWs;
};

inline const RVec3 & r(const Mesh::rw_t& rw) {return std::get<0>(rw);}
inline const double& w(const Mesh::rw_t& rw) {return std::get<1>(rw);}


#endif //_Mesh_H_

