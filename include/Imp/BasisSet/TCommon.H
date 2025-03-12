// File: BasisSet/TCommon.H
#ifndef _BasisSet_TCommon_H_
#define _BasisSet_TCommon_H_

#include <BasisSet.H>
class SlaterRadialIntegralTests;
class DiracIntegralTests;

class BasisSetImp
    : public virtual BasisSet
   // , public virtual Integrals_BS_2E<double>
{
public:
    typedef BasisSet Base;
    typedef Base::bs_t bs_t;
    typedef Base::const_iterator const_iterator;

    BasisSetImp();
    BasisSetImp(AnalyticIE<double>*);
    
    virtual ~BasisSetImp();
    virtual size_t GetNumFunctions() const;
    virtual void Insert(bs_t*);

    virtual const_iterator begin() const {return itsBasisSets.begin();}
    virtual const_iterator end  () const {return itsBasisSets.end  ();}
//
//  StreamableObject stuff.
//
    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
protected:
    friend class SlaterRadialIntegralTests;
    friend class DiracIntegralTests;

    AnalyticIE<double>* itsIE;
    optr_vector1<bs_t*> itsBasisSets;
};

#endif // _BasisSet_TCommon_H_
