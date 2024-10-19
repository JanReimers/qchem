// File: AngularMeshImp.C  mesh implementation



#include "Mesh/AngularMesh/AngularMeshImplementation.H"
#include "Mesh/Mesh.H"
#include "oml/io3d.h"
#include <iostream>
#include <iomanip>
#include <cassert>

AngularMeshImplementation::AngularMeshImplementation()
    {};

void AngularMeshImplementation::Initialize(const Vector<RVec3>& D, const Vector<double>& W)
{
    assert(W.size()==D     .size());
    for (auto i:D.indices()) itsRWs.push_back(std::make_tuple(D(i),W(i)));

    for (auto rw: itsRWs)
        if (norm(r(rw)) < 0.9999 || norm(r(rw))>1.0001)
            std::cerr << "AngularMeshImplementation::Initialize Direction " << r(rw) << " not as unit vector, mag=" << norm(r(rw)) << std::endl;
}


