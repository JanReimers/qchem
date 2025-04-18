// File: BasisSet/TCommon.H
#ifndef _BasisSet_TCommon_H_
#define _BasisSet_TCommon_H_

#include <BasisSet.H>
#include <Irrep_BS.H>
#include "Imp/Misc/UniqueID/UniqueIDImp.H"
class SlaterRadialIntegralTests;
class DiracIntegralTests;

class BS_Common
    : public virtual BasisSet
    , private UniqueIDImp
{
public:
    BS_Common() {};
    virtual ~BS_Common() {};
    virtual void Set(const LAParams&);

    virtual size_t GetNumFunctions() const;
    virtual symv_t GetSymmetries  () const;

//
//  StreamableObject stuff.
//
    virtual std::ostream&  Write(std::ostream&    ) const;
protected:
    virtual void Insert(bs_t*);

    virtual const_iterator begin() const {return itsBasisSets.begin();}
    virtual const_iterator end  () const {return itsBasisSets.end  ();}
    friend class SlaterRadialIntegralTests;
    friend class DiracIntegralTests;

    optr_vector1<bs_t*> itsBasisSets;
};

#endif // _BasisSet_TCommon_H_
