module;
#include <cassert>
#include <sstream>
#include <iomanip>
module qchem.Structure;
import qchem.Vector3D;

namespace qchem {

// A single atom's integration mesh is its atom-centred radial x angular grid.  The shared MakeMolecularMesh
// (Imp/Molecule.C) over a one-atom structure yields exactly that -- the Becke cell weight is identically 1
// with no other centres, so there is no multi-centre overhead.  (A dedicated pure-radial quadrature that
// exploits spherically-symmetric integrands would be more efficient still, but needs the AO's (l,m) exposed
// -- a later specialization.)
qcMesh::Mesh Atom::CreateIntegrationMesh(const qcMesh::MeshParams& mp) const
{
    return MakeMolecularMesh(*this, mp);
}

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
    : Atom(Site{Z,R,charge})
{};

Atom::Atom(const Site& s)
    : itsSite(s)
{
    assert(itsSite.Z>0);
    assert(itsSite.Z<150); //Maybe there is an island of stability at Z=140!!!!
};

std::string Atom::ID() const
{
    return SiteID(itsSite);
}

std::ostream& Atom::Write  (std::ostream& os) const
{
    return WriteSite(os,itsSite);
}

std::string SiteID(const Site& s)
{
    std::ostringstream os;
    os << "Z=" << s.Z << " R=" << s.R;
    return os.str();
}

std::ostream& WriteSite(std::ostream& os, const Site& s)
{
    os.setf(std::ios::fixed,std::ios::floatfield);
    os << std::setw(4) << s.Z << "    "
    << std::setw(5) << std::setprecision(2) << s.R << "     ";
    // os << std::endl; let the caller decide
    return os;
}

} // namespace qchem