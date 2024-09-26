// File: BlochQN.H  A Quantum Number for atomic (spherical) symmetry
#ifndef _BlochQN_H_
#define _BlochQN_H_



#include "BasisSet/QuantumNumber.H"
#include "oml/vector3d.h"

template <class T> class Vector3D;
typedef Vector3D<double> RVec3;

//---------------------------------------------------------------------------------
//
//  Translational symmetry, Bloch function wave vector.
//

class BlochQN
    : public QuantumNumber
{
public:
    BlochQN(         );
    BlochQN(RVec3 theK);

    virtual bool Match(const QuantumNumber&) const;
    virtual int GetDegeneracy() const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;
    virtual QuantumNumber* Clone(        ) const;

    RVec3   GetK() const
    {
        return itsK;
    }

private:
    RVec3 itsK;
};

#endif
