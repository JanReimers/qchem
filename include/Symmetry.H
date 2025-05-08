// File: Symmetry.H  Abstract interface for symmetries that do not include spin.
#ifndef _Symmetry_H_
#define _Symmetry_H_

#include "Imp/Misc/pmstream.h"

class Symmetry
    : public virtual PMStreamableObject
{
public:
    virtual ~Symmetry() {};
    virtual size_t SequenceIndex() const=0; //Used for op<
    virtual bool   Match     (const Symmetry&) const=0;
    virtual bool   MatchType (const Symmetry&) const=0;
    friend  bool   operator==(const Symmetry& a, const Symmetry& b)
    {
        return a.Match(b);
    }
    //! Does not include spin degeneracy which is handled separately
    virtual int GetDegeneracy     () const=0;
    virtual int GetPrincipleOffset() const=0; //Add to principle QN.  For atoms this is just l.
    std::string GetLabel          () const;
    virtual Symmetry* Clone() const=0;
    static  Symmetry* Factory(std::istream&);
};



#endif //_Symmetry_H_
