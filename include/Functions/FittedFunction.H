// File: FittedFunction.H  Linear fitted function interface
#ifndef _FittedFunction_H_
#define _FittedFunction_H_

#include "types.H"
#include "Functions/ScalarFunction.H"

//----------------------------------------------------------------------------------
//
//  Interface for a real space function that can be fit with a basis set.
//  The DoFit member function calculate the expansion coefficients.
//  The GetFunctionOverlap is supplied by the client (derived) class, and
//  returns the overlap of the funciont to be fitted with all the fit basis
//  functions.  Members for calculating the overlap of the fit with another
//  basis set are also supplied.
//
class FittedFunction
    : public virtual PlottableScalarFunction
{
public:
    virtual double DoFit(const FittedFunctionClient&)=0;
    virtual void   ReScale         (double factor               )      =0; //Fit *= factor
    virtual void   ShiftOrigin     (const RVec3& newCenter      )      =0;
    virtual void   FitMixIn        (const FittedFunction&,double)      =0; // this = this*(1-c) + that*c.
    virtual double FitGetChangeFrom(const FittedFunction&       ) const=0;
};

#endif //_FittedFunction_H_
