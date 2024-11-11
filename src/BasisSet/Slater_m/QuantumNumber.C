// File: Slater_m/QuantumNumer.H  Spherical harmonic Ylm symmetry.



#include "Imp/BasisSet/Slater_m/QuantumNumber.H"
#include <iostream>
#include <cassert>

YlmQN::YlmQN(): SphericalSymmetryQN(0),m(0) {};

YlmQN::YlmQN(int _l, int _m) : SphericalSymmetryQN(_l), m(_m) {};

bool YlmQN::Match(const QuantumNumber& qn) const
{
    const YlmQN* yqn = dynamic_cast<const YlmQN*>(&qn);
    assert(yqn);
    return itsL==yqn->itsL && m==yqn->m;
}

int YlmQN::GetDegeneracy() const
{
    return 1;
}

extern std::string SPDFG[];

std::ostream& YlmQN::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
        os << SPDFG[itsL] << "_" << m << " ";
    return os;
}

std::istream& YlmQN::Read (std::istream& is)
{
    return is;
}

QuantumNumber* YlmQN::Clone() const
{
    return new YlmQN(*this);
}

