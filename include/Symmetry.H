// File: QNs.h  Abstract interface for collection of Quantum Numbers.
#ifndef _QNs_H_
#define _QNs_H_

#include "Imp/Misc/pmstream.h"

class Symmetry
    : public virtual PMStreamableObject
{
public:
    virtual ~Symmetry() {};
    virtual size_t SequenceIndex() const=0; //Used for op<
    virtual bool Match(const Symmetry&) const=0;
    virtual bool MatchType(const Symmetry& b) const=0;
    friend  bool operator==(const Symmetry& a, const Symmetry& b)
    {
        return a.Match(b);
    }
    //! Does not include spin degeneracy which is handles separately
    virtual int GetDegeneracy() const=0;
    
    virtual Symmetry* Clone() const=0;
    static  Symmetry* Factory(std::istream&);
};



#endif //_QNs_H_
