// File: FittedFunctionClient.C  Linear fitted function interface
export module qchem.FittedFunctionClient;
import qchem.ScalarFunction;

import qchem.Fit_IBS;
// import qchem.BasisSet1.Fit_IBS;

//-------------------------------------------------------------------
//
//  Abstract interfaces used by the FittedFunction class.
//  This is all the fitting routines needs to know in order
//  to execute the LS fit.
//

//
//  Fit to a simple scalar function.  Requires numerical integration.
// 
export namespace qchem::Fitting
{

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
    typedef Fit_IBS fbs_t;
    virtual double FitGetConstraint() const=0;
    virtual rvec_t GetRepulsion3C(const fbs_t*) const=0;
};

} //namespace