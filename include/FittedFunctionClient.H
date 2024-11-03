// File: FittedFunctionClient.H  Linear fitted function interface
#ifndef _FittedFunctionClient_H_
#define _FittedFunctionClient_H_

#include <types.H>
#include <ScalarFunction.H>

//-------------------------------------------------------------------
//
//  Abstract interfaces used by the FittedFunction class.
//

//
//  Fit so a simple scalar function.  Requires numerical integration.
// 
class ScalarFFClient
{
public:
    virtual const ScalarFunction<double>* GetScalarFunction() const=0;

};

//
//  Fit to a density matrix function ro(r)=a(r)*b(r)*Dab 
//
class DensityFFClient
{
public:
    typedef Vector<double> RVec;
    virtual double FitGetConstraint() const=0;
    virtual RVec   GetRepulsion3C(const IrrepBasisSet* theFitBasisSet) const=0;
};

#endif //_FittedFunctionClient_H_
