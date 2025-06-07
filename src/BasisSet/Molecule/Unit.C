// File: UnitSymmetryQN  Quantum number for a raw basis set.



#include "Imp/BasisSet/Molecule/Unit.H"
#include <iostream>

UnitQN::UnitQN()
{};


int UnitQN::GetDegeneracy() const
{
    return 1;
}


std::ostream& UnitQN::Write(std::ostream& os) const
{
    return os;
}

std::istream& UnitQN::Read (std::istream& is)
{
    return is;
}


UnitnQN::UnitnQN() : n(0) {};
UnitnQN::UnitnQN(int _n) : n(_n) {};

std::ostream& UnitnQN::Write(std::ostream& os) const
{
    os << n << " ";
    return UnitQN::Write(os);
}

