// File: Symmetry/YlmQN.H  A Quantum Number for atomic (spherical) symmetry
#ifndef _Symmetry_YlQN_H_
#define _Symmetry_YlQN_H_



#include <QuantumNumber.H>
#include <tuple>

//---------------------------------------------------------------------------------
//
//  The relevant quantum number for spherical symmetry is the magnitude of the
//  orbital angular momentum L.
//

class YlQN
    : public virtual QuantumNumber
{
public:
    YlQN(        );
    YlQN(int theL);

    virtual bool Match(const QuantumNumber&) const;
    virtual int GetDegeneracy() const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;
    virtual QuantumNumber* Clone(        ) const;

    int     GetL() const {return itsL;}
    virtual std::pair<int,int> GetN(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const;

protected:
    int itsL;
    static const int LMax=3;
};

#endif //_Symmetry_YlQN_H_
