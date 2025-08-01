module;
#include <cassert>
#include <nlohmann/json.hpp>
module qchem.Atom;
import qchem.Cluster.AtomMesh;
import qchem.Mesh.Factory;
import oml.Vector3D;
using json = nlohmann::json;


Atom::Atom()
    : itsZ(0)
    , itsCharge(0)
    , itsR( )
{};

Atom::Atom(int Z, double charge)
    : itsZ(Z)
    , itsCharge(charge)
    , itsR(0,0,0)
{
    assert(itsZ>0);
    assert(itsZ<150); //Maybe there is an island of stability at Z=140!!!!
};

Atom::Atom(int Z, double charge, const RVec3& R)
    : itsZ(Z)
    , itsCharge(charge)
    , itsR(R)
{
    assert(itsZ>0);
    assert(itsZ<150); //Maybe there is an island of stability at Z=140!!!!
};

Mesh* Atom::CreateMesh(const MeshParams& mp) const
{
    json js={{"N",mp.Nradial},{"m",mp.MHL_m},{"alpha",mp.MHL_alpha}};
    RadialMesh* rm=MeshF::Factory(qchem::MHL,js);
    js={{"Nangle",mp.Nangle}};
    Mesh* am=MeshF::Factory(qchem::Gauss,js);  
    return new AtomMesh(rm,am,itsR); 
}


double Atom::GetNumElectrons() const
{
    return itsZ-itsCharge;
}

std::ostream& Atom::Write  (std::ostream& os) const
{
    os.setf(std::ios::fixed,std::ios::floatfield);
    os << std::setw(4) << itsZ << "    "
    << std::setw(5) << std::setprecision(2) << itsR << "     ";
    // os << std::endl; let the caller decide
    return os;
}


