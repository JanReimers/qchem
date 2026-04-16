module;
#include <iomanip>
#include <cassert>
#include <memory>

module qchem.Cluster;
import qchem.Cluster.MoleculeMesh;
import qchem.stl_io;
import qchem.Streamable;

Molecule::Molecule(const Molecule& m)
{
    for (auto a:m)
    {
        itsAtoms.push_back(new Atom(a->itsZ,a->GetNetCharge(),a->itsR));
    }
}

Molecule::~Molecule()
{
    for (auto a:itsAtoms) delete a;
}

void Molecule::Insert(Atom* a, double charge)
{
    itsAtoms.push_back(a);
    itsCharge+=charge;
}

size_t Molecule::GetNumAtoms() const
{
    return itsAtoms.size();
}

Mesh*  Molecule::CreateMesh(const MeshParams& mp) const
{
    return new MoleculeMesh(*this,mp);
}

std::ostream& Molecule::Write(std::ostream& os) const
{
    os << "Molecule with " << GetNumAtoms() << " atoms"
    << ", nuclear charge " << GetNuclearCharge() << "(e)"
    << ", net charge "<< GetNetCharge() << "(e)" << std::endl;
    os << "Atom #  Element  Position vector" << std::endl;
    int i=1;
    for (auto& b:*this) os << std::setw(5) << i++ << "   " << *b << std::endl;
    // os << std::endl;

    return os;
}

