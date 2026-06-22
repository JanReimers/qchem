module;
#include <cassert>
#include <nlohmann/json.hpp>
module qchem.Structure;
import qchem.Structure.AtomMesh;
import qchem.Mesh.Factory;
import qchem.Vector3D;
using json = nlohmann::json;


Atom::Atom(int Z)
    : Atom(Z,0.0,{0,0,0})
{};
Atom::Atom(int Z, double charge)
    : Atom(Z,charge,{0,0,0})
{};
Atom::Atom(int Z, const rvec3_t& R)
    : Atom(Z,0.0,R)
{};

Atom::Atom(int Z, double charge, const rvec3_t& R)
    : itsZ(Z)
    , itsCharge(charge)
    , itsR(R)
{
    assert(itsZ>0);
    assert(itsZ<150); //Maybe there is an island of stability at Z=140!!!!
    dummy.push_back(this);
};

Atom::Atom(const Atom& a) : Atom(a.itsZ,a.itsCharge,a.itsR) {}

Mesh* Atom::CreateMesh(const MeshParams& mp) const
{
    json js={{"N",mp.Nradial},{"m",mp.MHL_m},{"alpha",mp.MHL_alpha}};
    RadialMesh* rm=MeshF::Factory(qchem::MHL,js);
    js={{"Nangle",mp.Nangle}};
    Mesh* am=MeshF::Factory(qchem::Gauss,js);  
    return new AtomMesh(rm,am,itsR); 
}


std::string Atom::ID() const
{
    std::ostringstream os;
    os << "Z=" << itsZ << " R=" << itsR;
    return os.str();
}

std::ostream& Atom::Write  (std::ostream& os) const
{
    os.setf(std::ios::fixed,std::ios::floatfield);
    os << std::setw(4) << itsZ << "    "
    << std::setw(5) << std::setprecision(2) << itsR << "     ";
    // os << std::endl; let the caller decide
    return os;
}


