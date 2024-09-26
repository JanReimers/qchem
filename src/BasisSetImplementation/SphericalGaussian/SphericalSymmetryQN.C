// File: SphericalSymmetryQN.N  A Quantum Number for atomic (spherical) symmetry



#include "BasisSetImplementation/SphericalGaussian/SphericalSymmetryQN.H"
#include "oml/imp/binio.h"
#include <iostream>
#include <cassert>

SphericalSymmetryQN::SphericalSymmetryQN()
    : itsL(0)
{};

SphericalSymmetryQN::SphericalSymmetryQN(int theL)
    : itsL(theL)
{};

bool SphericalSymmetryQN::Match(const QuantumNumber& qn) const
{
    const SphericalSymmetryQN* aqn = dynamic_cast<const SphericalSymmetryQN*>(&qn);
    assert(aqn);
    return itsL==aqn->itsL;
}

int SphericalSymmetryQN::GetDegeneracy() const
{
    return 2*itsL+1;
}

std::ostream& SphericalSymmetryQN::Write(std::ostream& os) const
{
    UniqueID::Write(os);
    if (StreamableObject::Binary())
        BinaryWrite(itsL,os);
    else
        os << itsL << " ";
    return os;
}

std::istream& SphericalSymmetryQN::Read (std::istream& is)
{
    UniqueID::Read(is);
    if (StreamableObject::Binary())
        BinaryRead(itsL,is);
    else
    {
        is >> itsL;
    }
    return is;
}

QuantumNumber* SphericalSymmetryQN::Clone() const
{
    return new SphericalSymmetryQN(*this);
}

