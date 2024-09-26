// File: TBasisSet.H  Abstract interface for a real/complex valued basis set.
#ifndef _TBasisSet_H_
#define _TBasisSet_H_



#include "BasisSet/BasisSet.H"
#include "Functions/VectorFunction.H"

template <class T> class IntegralDataBase;
template <class T> class EigenSolver;

//----------------------------------------------------------------------------
//
//  Extend basis to be a set of real or complex valued functions
//
template <class T> class TBasisSet
    : public virtual BasisSet
    , public virtual VectorFunction<T>
{
public:
    index_t GetVectorSize() const
    {
        return GetNumFunctions();
    }
    virtual SMat GetRepulsion(const SMat& Dcd, const TBasisSet* bs_cd) const=0;
    virtual SMat GetExchange (const SMat& Dcd, const TBasisSet* bs_cd) const=0;

    virtual IntegralDataBase<T>* GetDataBase   () const=0;
    virtual EigenSolver<T>*      GetEigenSolver() const=0;

};

#endif //_TBasisSet_H_
