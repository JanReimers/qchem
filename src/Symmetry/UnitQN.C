// File: UnitSymmetryQN  Quantum number for a raw basis set.



#include "Imp/Symmetry/UnitQN.H"
#include <iostream>

UnitQN::UnitQN()
{};

int UnitQN::GetDegeneracy() const
{
    return 1;
}

QNs* UnitQN::AddPrincipleQN(int index) const
{
    return new UnitnQN(index);
}


std::ostream& UnitQN::Write(std::ostream& os) const
{
    return os;
}

std::istream& UnitQN::Read (std::istream& is)
{
    return is;
}

QNs* UnitQN::Clone() const
{
    return new UnitQN;
}

UnitnQN::UnitnQN() : n(0) {};
UnitnQN::UnitnQN(int _n) : n(_n) {};

std::ostream& UnitnQN::Write(std::ostream& os) const
{
    os << n << " ";
    return UnitQN::Write(os);
}

QNs* UnitnQN::Clone() const
{
    return new UnitnQN(*this);
}
