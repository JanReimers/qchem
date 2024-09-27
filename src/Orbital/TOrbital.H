// File: TOrbital.H  Interface for Orbital functions.
#ifndef _TOrbital_H_
#define _TOrbital_H_



#include "Functions/ScalarFunction.H"
#include "Orbital/Orbital.H"

//---------------------------------------------------------
//
//  Templated depending or whether it is a real or
//  complex valued orbital.
//
template <class T> class TOrbital
    : public virtual Orbital
    , public virtual ScalarFunction<T>
{
    typedef typename ScalarFunction<T>::SMat SMat;
public:
    virtual void AddDensityMatrix(SMat&) const=0;
};

#endif //_TOrbital_H_
