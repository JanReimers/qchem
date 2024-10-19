// File: MeshImplementation.C  mesh implementation



#include "Mesh/MeshImplementation.H"
#include "oml/io3d.h"
#include <iostream>
#include <iomanip>

//
//  The full mesh is just a direct product of radial and ungular meshes.
//
MeshImplementation::MeshImplementation(index_t numPoints)
{};

void MeshImplementation::Initialize(const Vector<RVec3>& R,const Vector<double>& W)
{
    for (auto i:R.indices()) itsRWs.push_back(std::make_tuple(R(i),W(i)));
}

std::ostream& MeshImplementation::Write(std::ostream& os) const
{
    if (!StreamableObject::Pretty())
    {
        UniqueID::Write(os);
       // os << itsRWs;
    }
    else
    {
        std::ios::fmtflags f=os.setf(std::ios::left,std::ios::adjustfield);
//    os << setw(15) << Type();
        os.flags(f);
        os << std::setw(8) << size();
    }
    return os;
}

std::istream& MeshImplementation::Read(std::istream& is)
{
    UniqueID::Read(is);
    //is >> itsRWs;
    return is;
}

void MeshImplementation::ShiftOrigin(const RVec3& r)
{
    for (auto& rw:itsRWs) std::get<0>(rw)+=r;
    NewID();
}

