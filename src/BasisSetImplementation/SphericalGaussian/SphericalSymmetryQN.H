// File: SphericalSymmetryQN.H  A Quantum Number for atomic (spherical) symmetry
#ifndef _SphericalSymmetryQN_H_
#define _SphericalSymmetryQN_H_



#include "BasisSet/QuantumNumber.H"

//---------------------------------------------------------------------------------
//
//  The relevant quantum number for spherical symmetry is the magnitude of the
//  orbital angular momentum L.
//

class SphericalSymmetryQN
    : public QuantumNumber
{
public:
    SphericalSymmetryQN(        );
    SphericalSymmetryQN(int theL);

    virtual bool Match(const QuantumNumber&) const;
    virtual int GetDegeneracy() const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;
    virtual QuantumNumber* Clone(        ) const;

    int     GetL() const
    {
        return itsL;
    }

private:
    int itsL;
};

#endif
