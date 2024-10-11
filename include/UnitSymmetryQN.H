// File: UnitSymmetryQN.H  Abstract interface for a Quantum Number
#ifndef _UnitSymmetryQN_H_
#define _UnitSymmetryQN_H_



#include "QuantumNumber.H"

class UnitSymmetryQN
    : public virtual QuantumNumber
{
public:
    UnitSymmetryQN();
    virtual bool Match(const QuantumNumber&) const
    {
        return true;
    }
    virtual int GetDegeneracy() const;

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;
    virtual QuantumNumber* Clone() const;
};

#endif
