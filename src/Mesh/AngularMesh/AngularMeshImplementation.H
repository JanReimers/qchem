// File: AngularMeshImp.H  Angular Mesh implemented as an vector of directions and weights.
#ifndef _AngularMeshImp_H_
#define _AngularMeshImp_H_

#include "Mesh/AngularMesh/AngularMesh.H"
#include "oml/vector3d.h"
#include "oml/vector.h"

class AngularMeshImplementation : public virtual AngularMesh
{
public:
    AngularMeshImplementation(index_t NumDirections);

    virtual size_t GetNumDirections () const
    {
        return itsDirections.size();
    }

protected:
    void Initialize(const Vector<RVec3>& D, const Vector<double>& W);

private:
    friend  class AngularMeshBrowser;

    Vector<RVec3>  itsDirections;
    Vector<double> itsWeights;
};

#endif // _AngularMeshImp_H_
