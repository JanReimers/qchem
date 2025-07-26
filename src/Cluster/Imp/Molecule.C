module;
#include <iomanip>
#include <cassert>
#include <memory>

module qchem.Molecule;
import qchem.Cluster.MoleculeMesh;
import qchem.stl_io;
import qchem.Streamable;

Molecule::Molecule()
    : itsNumElectrons(0)
    , itsAtoms    ( )
{};

void Molecule::Insert(Atom* a)
{
    itsAtoms.push_back(std::unique_ptr<Atom>(a));
    itsNumElectrons+=a->GetNumElectrons();
}

size_t Molecule::GetNumAtoms() const
{
    return itsAtoms.size();
}

int Molecule::GetNuclearCharge() const
{
    int chg=0;
    for(auto& b:*this) chg+=b->itsZ;
    return chg;
}

double Molecule::GetNetCharge() const
{
    return GetNuclearCharge()-itsNumElectrons;
}

double Molecule::GetNumElectrons() const
{
    return itsNumElectrons;
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

