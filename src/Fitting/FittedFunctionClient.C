// File: FittedFunctionClient.C  Linear fitted function interface
export module qchem.FittedFunctionClient;
import qchem.Fit_IBS;
import qchem.ScalarFunction;
//-------------------------------------------------------------------
//
//  Abstract interfaces used by the FittedFunction class.
//  This is all the fitting routines needs to know in order
//  to execute the LS fit.
//

//
//  Fit to a simple scalar function.  Requires numerical integration.
// 
export class ScalarFFClient
{
public:
    virtual const ScalarFunction<double>* GetScalarFunction() const=0;

};

//
//  Fit to a density matrix function like ro(r)=a(r)*b(r)*Dab 
//
export class DensityFFClient
{
public:
    typedef Vector<double> RVec;
    virtual double FitGetConstraint() const=0;
    virtual RVec   GetRepulsion3C(const Fit_IBS*) const=0;
};

