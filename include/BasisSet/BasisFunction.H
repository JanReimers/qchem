// File: BasisFunction.H
#ifndef _BasisFunction_H_
#define _BasisFunction_H_

#include <iosfwd>
import qchem.Streamable;

//--------------------------------------------------------------
//
//  This interface of the basis function provides no notion
//  of weather the function is real or complex (or vector!)
//  valued.  Normalization and total charge should always
//  be real.
//
class BasisFunction
    : public virtual Streamable
{
public:
    virtual ~BasisFunction()  {};

    virtual BasisFunction* Clone  (             ) const=0;
    static  BasisFunction* Factory(std::istream&)        ;
};

import qchem.ScalarFunction;
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

#endif //_BasisFunction_H_
