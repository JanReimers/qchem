// File: AngularMeshImp.H  Angular Mesh implemented as an vector of directions and weights.
#ifndef _AngularMeshImp_H_
#define _AngularMeshImp_H_

#include "Mesh/AngularMesh/AngularMesh.H"
#include "oml/vector3d.h"
#include "oml/vector.h"

class AngularMeshImplementation : public virtual AngularMesh
{
public:
    AngularMeshImplementation();

    virtual const_iterator begin() const {return itsRWs.begin();}
    virtual const_iterator end  () const {return itsRWs.end  ();}
    virtual index_t        size () const {return itsRWs.size ();}

protected:
    void Initialize(const Vector<RVec3>& D, const Vector<double>& W);

private:
    vec_t itsRWs;
};

#endif // _AngularMeshImp_H_
