// File: UnitSymmetryQN  Quantum number for a raw basis set.



#include "Imp/BasisSet/Molecule/Unit.H"
#include <iostream>

UnitQN::UnitQN()
{};

bool UnitQN::MatchType(const Symmetry& b) const
{
    return dynamic_cast<const UnitQN*>(&b)!=0;
}

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

Symmetry* UnitQN::Clone() const
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

Symmetry* UnitnQN::Clone() const
{
    return new UnitnQN(*this);
}
