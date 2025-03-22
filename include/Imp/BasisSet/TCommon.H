// File: BasisSet/TCommon.H
#ifndef _BasisSet_TCommon_H_
#define _BasisSet_TCommon_H_

#include <BasisSet.H>
#include "Imp/Misc/UniqueID/UniqueIDImp.H"
class SlaterRadialIntegralTests;
class DiracIntegralTests;

class BasisSetImp
    : public virtual BasisSet
    , private UniqueIDImp
{
public:
    typedef BasisSet Base;
    typedef Base::bs_t bs_t;
    typedef Base::const_iterator const_iterator;

    BasisSetImp() {};
    virtual ~BasisSetImp() {};
    virtual void Set(const LAParams&);

    virtual size_t GetNumFunctions() const;
    virtual void Insert(bs_t*);

//
//  StreamableObject stuff.
//
    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
protected:
    virtual const_iterator begin() const {return itsBasisSets.begin();}
    virtual const_iterator end  () const {return itsBasisSets.end  ();}
    friend class SlaterRadialIntegralTests;
    friend class DiracIntegralTests;

    optr_vector1<bs_t*> itsBasisSets;
};

#endif // _BasisSet_TCommon_H_
