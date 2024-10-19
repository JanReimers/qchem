// File: BlochQN.N  A Quantum Number for atomic (spherical) symmetry



#include "BasisSetImplementation/PlaneWave/BlochQN.H"
#include "oml/imp/binio.h"
#include "oml/io3d.h"
#include <iostream>
#include <cassert>

BlochQN::BlochQN()
    : itsK(0,0,0)
{};

BlochQN::BlochQN(RVec3 theK)
    : itsK(theK)
{};

bool BlochQN::Match(const QuantumNumber& qn) const
{
    const BlochQN* aqn = dynamic_cast<const BlochQN*>(&qn);
    assert(aqn);
    return itsK==aqn->itsK;
}

int BlochQN::GetDegeneracy() const
{
    return 1;
}

std::ostream& BlochQN::Write(std::ostream& os) const
{
    UniqueID::Write(os);
    if (StreamableObject::Binary())
        os << itsK;
    else
        os << itsK << " ";
    return os;
}

std::istream& BlochQN::Read (std::istream& is)
{
    UniqueID::Read(is);
    return is >> itsK;
}

QuantumNumber* BlochQN::Clone() const
{
    return new BlochQN(*this);
}

