// File: FittedFunctionClient.H  Linear fitted function interface
#ifndef _FittedFunctionClient_H_
#define _FittedFunctionClient_H_

#include <types.H>
#include <Common/ScalarFunction.H>

//-------------------------------------------------------------------
//
//  Abstract interfaces used by the FittedFunction class.
//  This is all the fitting routines needs to know in order
//  to execute the LS fit.
//

//
//  Fit to a simple scalar function.  Requires numerical integration.
// 
class ScalarFFClient
{
public:
    virtual const ScalarFunction<double>* GetScalarFunction() const=0;

};

//
//  Fit to a density matrix function like ro(r)=a(r)*b(r)*Dab 
//
class DensityFFClient
{
public:
    typedef Vector<double> RVec;
    virtual double FitGetConstraint() const=0;
    virtual RVec   GetRepulsion3C(const Fit_IBS*) const=0;
};

#endif //_FittedFunctionClient_H_
