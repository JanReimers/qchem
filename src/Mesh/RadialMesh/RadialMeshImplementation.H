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

    virtual const_iterator begin() const {return itsRWs.begin();}
    virtual const_iterator end  () const {return itsRWs.end  ();}
    virtual index_t        size () const {return itsRWs.size ();}
    
protected:
    void Initialize(const Vector<double>& R,const Vector<double>& W);

private:

    Vector<double> itsPoints;
    Vector<double> itsWeights;
    
    vec_t itsRWs;
};

#endif // _RadialMeshImplementation_H_
