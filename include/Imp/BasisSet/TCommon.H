// File: BasisSet/TCommon.H
#ifndef _BasisSet_TCommon_H_
#define _BasisSet_TCommon_H_

#include <BasisSet.H>
class SlaterRadialIntegralTests;

class BasisSetImp
    : public virtual BasisSet
{
public:
    typedef BasisSet Base;
    typedef Base::const_iterator const_iterator;

    BasisSetImp();
    BasisSetImp(AnalyticIE<double>*);
    
    virtual ~BasisSetImp();
    virtual size_t GetNumFunctions() const;
    virtual size_t GetNumIrreps() const;
    virtual void Insert(IrrepBasisSet*);

    virtual const_iterator begin() const {return itsBasisSets.begin();}
    virtual const_iterator end  () const {return itsBasisSets.end  ();}
    virtual IntegralDataBase<double>* GetDataBase() const {return itsDB;}
//
//  StreamableObject stuff.
//
    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
protected:
    friend class SlaterRadialIntegralTests;
    AnalyticIE<double>* itsIE;
private:
    mutable IntegralDataBase<double>* itsDB;
    optr_vector1<IrrepBasisSet*> itsBasisSets;
};

#endif // _BasisSet_TCommon_H_
