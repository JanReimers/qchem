// File: RadialMeshImplementation.C  mesh implementation



#include "Mesh/RadialMesh/RadialMeshImplementation.H"
#include <cassert>

RadialMeshImplementation::RadialMeshImplementation()
    
{};

void RadialMeshImplementation::Initialize(const Vector<double>& R, const Vector<double>& W)
{
    assert(R.size()==W.size());
    for (auto i:R.indices()) itsRWs.push_back(std::make_tuple(R(i),W(i)));
}

