// File: AngularMeshImp.C  mesh implementation



#include "Mesh/AngularMesh/AngularMeshImplementation.H"
#include "Mesh/AngularMesh/AngularMeshBrowser.H"
#include "oml/io3d.h"
#include <iostream>
#include <iomanip>
#include <cassert>

AngularMeshImplementation::AngularMeshImplementation(index_t NumDirections)
    : itsDirections(NumDirections)
    , itsWeights   (NumDirections)
{};

void AngularMeshImplementation::Initialize(const Vector<RVec3>& D, const Vector<double>& W)
{
    assert(D.size()==itsDirections.size());
    assert(W.size()==D            .size());
    itsDirections=D;
    itsWeights   =W;
    Vector<RVec3>::const_iterator b(D.begin());
    for (; b!=D.end(); b++)
        if (!*b < 0.9999 || !*b>1.0001)
            std::cerr << "AngularMeshImplementation::Initialize Direction " << *b << " not as unit vector, mag=" << !*b << std::endl;
}


