// File: TBasisFunction.H  Real or complex valued basis function interface.
#ifndef _TBasisFunction_H_
#define _TBasisFunction_H_



#include "BasisFunction.H"
#include "Functions/ScalarFunction.H"

//---------------------------------------------------------------
//
//  Here the basis becomes a real or complex valued function
//  depending of T.
//
template <class T> class TBasisFunction
    : public virtual BasisFunction
    , public virtual ScalarFunction<T>
{
};


#endif //_TBasisFunction_H_
