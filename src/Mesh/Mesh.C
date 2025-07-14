//  File: Mesh.C  mesh implementation
module;
#include "oml/vector3d.h"
#include <tuple>
#include <vector>

export module Mesh;

export class Mesh 
{ 
public:
    typedef Vector3D<double> RVec3; 
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

export inline const Mesh::RVec3 & r(const Mesh::rw_t& rw) {return std::get<0>(rw);}
export inline const double& w(const Mesh::rw_t& rw) {return std::get<1>(rw);}


void Mesh::ShiftOrigin(const RVec3& r)
{
    for (auto& rw:itsRWs) std::get<0>(rw)+=r;
}

Mesh*   Mesh::Clone      () const
{
    return new Mesh(*this);
}
