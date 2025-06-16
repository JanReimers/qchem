//  File: Mesh.C  mesh implementation

#include <Mesh/Mesh.H>

void Mesh::ShiftOrigin(const RVec3& r)
{
    for (auto& rw:itsRWs) std::get<0>(rw)+=r;
}

Mesh*   Mesh::Clone      () const
{
    return new Mesh(*this);
}
