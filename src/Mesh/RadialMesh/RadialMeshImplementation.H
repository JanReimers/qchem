// File: RadialMeshImplementation.H  Radial Mesh Implementation, radial distances and weights.
#ifndef _RadialMeshImplementation_H_
#define _RadialMeshImplementation_H_



#include "Mesh/RadialMesh/RadialMesh.H"
#include "oml/vector.h"

class RadialMeshImplementation
    : public virtual RadialMesh
{
public:
    RadialMeshImplementation();

    virtual const_iterator begin() const {return itsRWs.begin();}
    virtual const_iterator end  () const {return itsRWs.end  ();}
    virtual index_t        size () const {return itsRWs.size ();}
    
protected:
    void push_back(const double& r, const double& w) {itsRWs.push_back(std::make_tuple(r,w));}

private:
    vec_t itsRWs;
};

#endif // _RadialMeshImplementation_H_
