// File: MeshImplementation.C  mesh implementation



#include "Mesh/MeshImplementation.H"
//#include "oml/vector_io.h"
#include "oml/io3d.h"
#include <iostream>
#include <iomanip>

//
//  The full mesh is just a direct product of radial and ungular meshes.
//
MeshImplementation::MeshImplementation(index_t numPoints)
    : itsPoints (numPoints)
    , itsWeights(numPoints)
{};

void MeshImplementation::Initialize(const Vector<RVec3>& R,const Vector<double>& W)
{
    itsPoints=R;
    itsWeights=W;
}

std::ostream& MeshImplementation::Write(std::ostream& os) const
{
    if (!StreamableObject::Pretty())
    {
        UniqueID::Write(os);
        os << itsPoints << itsWeights;
    }
    else
    {
        std::ios::fmtflags f=os.setf(std::ios::left,std::ios::adjustfield);
//    os << setw(15) << Type();
        os.flags(f);
        os << std::setw(8) << itsPoints.size();
    }
    return os;
}

std::istream& MeshImplementation::Read(std::istream& is)
{
    UniqueID::Read(is);
    is >> itsPoints >> itsWeights;
    return is;
}

void MeshImplementation::ShiftOrigin(const RVec3& r)
{
    Vector<RVec3>::iterator i(itsPoints.begin());
    for(; i!=itsPoints.end(); i++) (*i)+=r;
    NewID();
}

