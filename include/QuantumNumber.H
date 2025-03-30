// File: QNs.h  Abstract interface for collection of Quantum Numbers.
#ifndef _QNs_H_
#define _QNs_H_

#include "Imp/Misc/pmstream.h"

class QNs
    : public virtual PMStreamableObject
{
public:
    virtual ~QNs() {};
    virtual size_t SequenceIndex() const=0; //Used for op<
    virtual bool Match(const QNs&) const=0;
    virtual bool MatchType(const QNs& b) const=0;
    friend  bool operator==(const QNs& a, const QNs& b)
    {
        return a.Match(b);
    }
    //! Does not include spin degeneracy which is handles separately
    virtual int GetDegeneracy() const=0;
    //! Combine with a principle QN from eigen levels with this Irrep.
    virtual QNs* AddPrincipleQN(int index) const=0;
    
    virtual QNs* Clone() const=0;
    static  QNs* Factory(std::istream&);
};



#endif //_QNs_H_
