// File: RadialMeshImplementation.H  Radial Mesh Implementation, radial distances and weights.
#ifndef _RadialMeshImplementation_H_
#define _RadialMeshImplementation_H_



#include "Mesh/RadialMesh/RadialMesh.H"
#include "oml/vector.h"

class RadialMeshImplementation
    : public virtual RadialMesh
{
public:
    RadialMeshImplementation(index_t NumPoints);

    virtual index_t  GetNumPoints () const
    {
        return itsPoints.size();
    }

protected:
    void Initialize(const Vector<double>& R,const Vector<double>& W);

private:
    friend  class RadialMeshBrowser;

    Vector<double> itsPoints;
    Vector<double> itsWeights;
};

#endif // _RadialMeshImplementation_H_
