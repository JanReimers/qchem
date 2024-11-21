// File: UnitSymmetryQN.H  Abstract interface for a Quantum Number
#ifndef _UnitSymmetryQN_H_
#define _UnitSymmetryQN_H_



#include <QuantumNumber.H>

class UnitQN
    : public virtual QuantumNumber
{
public:
    UnitQN();
    virtual bool Match(const QuantumNumber&) const
    {
        return true;
    }
    virtual int GetDegeneracy() const;
    virtual QuantumNumber* AddPrincipleQN(int index) const;

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;
    virtual QuantumNumber* Clone() const;
};

class UnitnQN
    : public virtual QuantumNumber
    , UnitQN
{
public:
    UnitnQN() ;
    UnitnQN(int n);

    virtual bool Match(const QuantumNumber&) const
    {
        return true;
    }
    virtual std::ostream& Write(std::ostream&) const;
    virtual QuantumNumber* Clone() const;
private:
    int n; //Principle QN
};

#endif
