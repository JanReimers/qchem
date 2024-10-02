// File: TBasisSet.H  Abstract interface for a real/complex valued basis set.
#ifndef _TBasisSet_H_
#define _TBasisSet_H_

#include "types.H"
#include "BasisSet.H"
#include "Functions/VectorFunction.H"
#include "Misc/ptr_vector1.h"
//----------------------------------------------------------------------------
//
//  Extend basis to be a set of real or complex valued functions
//
template <class T> class TBasisSet
    : public virtual BasisSet
    , public virtual VectorFunction<T>
{
public:
    typedef dynamic_cast_iterator<BasisFunction*,TBasisFunction<T>* > const_iterator;

    size_t GetVectorSize() const
    {
        return GetNumFunctions();
    }
    virtual SMat GetRepulsion(const FittedFunction*) const=0; //unhide
    virtual SMat GetRepulsion(const SMat& Dcd, const TBasisSet* bs_cd) const=0;
    virtual SMat GetExchange (const SMat& Dcd, const TBasisSet* bs_cd) const=0;

    virtual IntegralDataBase<T>* GetDataBase() const=0;
    //virtual LASolver<T>*         GetLASolver() const=0;
   
    const_iterator beginT() const {return const_iterator(begin());};
};

#endif //_TBasisSet_H_
