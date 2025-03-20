// File: qnum.h  Abstract interface for a Quantum Number
#ifndef _qnum_h_
#define _qnum_h_

#include "Imp/Misc/pmstream.h"

class QuantumNumber
    : public virtual PMStreamableObject
{
public:
    virtual ~QuantumNumber() {};
    virtual bool Match(const QuantumNumber&) const=0;
    friend  bool operator==(const QuantumNumber& q1, const QuantumNumber& q2)
    {
        return q1.Match(q2);
    }
    //! Does not include spin degeneracy which is handles separately
    virtual int GetDegeneracy() const=0;
    //! Combine with a principle QN from eigen levels with this Irrep.
    virtual QuantumNumber* AddPrincipleQN(int index) const=0;
    
    virtual QuantumNumber* Clone() const=0;
    static  QuantumNumber* Factory(std::istream&);
};


#endif
