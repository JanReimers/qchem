// File: UnitSymmetryQN  Quantum number for a raw basis set.



#include "BasisSetImplementation/UnitSymmetryQN.H"

UnitSymmetryQN::UnitSymmetryQN()
{};

int UnitSymmetryQN::GetDegeneracy() const
{
    return 1;
}


std::ostream& UnitSymmetryQN::Write(std::ostream& os) const
{
    UniqueID::Write(os);
    return os;
}

std::istream& UnitSymmetryQN::Read (std::istream& is)
{
    UniqueID::Read(is);
    return is;
}

QuantumNumber* UnitSymmetryQN::Clone() const
{
    return new UnitSymmetryQN;
}

