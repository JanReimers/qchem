module;
#include <iomanip>
#include <cassert>
#include <memory>

module qchem.Structure;
import qchem.stl_io;
import qchem.Streamable;

Molecule::Molecule(const Structure& atoms)
{
    for (auto a:atoms)
    {
        itsAtoms.push_back(new Atom(*a));
    }
}

Molecule::~Molecule()
{
    for (auto a:itsAtoms) delete a;
}

void Molecule::Insert(Atom* a)
{
    itsAtoms.push_back(a);
}

size_t Molecule::GetNumAtoms() const
{
    return itsAtoms.size();
}

std::ostream& Molecule::Write(std::ostream& os) const
{
    os << "Molecule with " << GetNumAtoms() << " atoms"
    << ", nuclear charge " << GetNuclearCharge() << "(e)"
    << ", net charge "<< GetNetCharge() << "(e)" << std::endl;
    os << "Atom #  Element  Position vector" << std::endl;
    int i=1;
    for (auto b:*this) os << std::setw(5) << i++ << "   " << *b << std::endl;
    // os << std::endl;

    return os;
}

