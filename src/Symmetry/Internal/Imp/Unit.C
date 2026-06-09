module;
#include <iostream>
module qchem.Symmetry.Unit;

UnitQN::UnitQN()
{};


size_t UnitQN::GetDegeneracy() const
{
    return 1;
}


std::ostream& UnitQN::Write(std::ostream& os) const
{
    return os << "No Symmetry";
}

