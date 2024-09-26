// File: RadialMeshImplementation.C  mesh implementation



#include "Mesh/RadialMesh/RadialMeshImplementation.H"
#include <cassert>

RadialMeshImplementation::RadialMeshImplementation(index_t NumPoints)
    : itsPoints (NumPoints)
    , itsWeights(NumPoints)
{};

void RadialMeshImplementation::Initialize(const Vector<double>& R, const Vector<double>& W)
{
    assert(R.size()==W        .size());
    assert(R.size()==itsPoints.size());
    itsPoints =R;
    itsWeights=W;
}

