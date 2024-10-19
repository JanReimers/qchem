// File: MeshImplementation.H   Mesh implemented as an vector of MeshPoint's.
#ifndef _MeshImplementation_H_
#define _MeshImplementation_H_

#include "Mesh/Mesh.H"
#include "oml/vector3d.h"
#include "oml/vector.h"

template <class T> class MeshIntegrator;

class MeshImplementation : public virtual Mesh
{
public:
    MeshImplementation(                 ) {};
   
    virtual void    ShiftOrigin (const RVec3&);
    
    virtual const_iterator begin() const {return itsRWs.begin();}
    virtual const_iterator end  () const {return itsRWs.end  ();}
    virtual index_t        size () const {return itsRWs.size ();}

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

protected:
    void push_back(const RVec3& r, const double& w) {itsRWs.push_back(std::make_tuple(r,w));}
private:
    vec_t itsRWs;
};

#endif // _MeshImplementation_H_
