// File: UnitSymmetryQN.H  Abstract interface for a Quantum Number
#ifndef _UnitSymmetryQN_H_
#define _UnitSymmetryQN_H_



#include <QuantumNumber.H>

class UnitQN
    : public virtual QNs
{
public:
    UnitQN();
    virtual bool Match(const QNs&) const
    {
        return true;
    }
    virtual int GetDegeneracy() const;
    virtual QNs* AddPrincipleQN(int index) const;

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;
    virtual QNs* Clone() const;
};

class UnitnQN
    : public virtual QNs
    , UnitQN
{
public:
    UnitnQN() ;
    UnitnQN(int n);

    virtual bool Match(const QNs&) const
    {
        return true;
    }
    virtual std::ostream& Write(std::ostream&) const;
    virtual QNs* Clone() const;
private:
    int n; //Principle QN
};

#endif
