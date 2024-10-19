// File: qnum.h  Abstract interface for a Quantum Number
#ifndef _qnum_h_
#define _qnum_h_



#include "Imp/Misc/Persistent/Persistent.H"

class QuantumNumber
    : public virtual Persistent
{
public:
    virtual bool Match(const QuantumNumber&) const=0;
    friend  bool operator==(const QuantumNumber& q1, const QuantumNumber& q2)
    {
        return q1.Match(q2);
    }
    virtual int GetDegeneracy() const=0;

    virtual QuantumNumber* Clone() const=0;
    static  QuantumNumber* Factory(std::istream&);
};

#endif
