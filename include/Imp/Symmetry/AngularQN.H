// File: AngularQN.H A Quantum Number for atomic (spherical) symmetry
#ifndef _Symmetry_AngularQN_H_
#define _Symmetry_AngularQN_H_

#include <QuantumNumber.H>
#include <tuple>
//---------------------------------------------------------------------------------
//
// AngularQN for atoms
//

class AngularQN
    : public virtual QuantumNumber
{
public:
    virtual std::pair<int,int> GetN(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const=0;
    virtual int GetL() const=0;
    virtual AngularQN* Clone() const=0;
};

#endif //_Symmetry_AngularQN_H_