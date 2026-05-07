// File: FittedFunctionClient.C  Linear fitted function interface
export module qchem.FittedFunctionClient;
import qchem.Fitting.Types;
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
    virtual double FitGetConstraint() const=0;
    virtual rvec_t GetRepulsion3C(const fbs_t*) const=0;
};

} //namespace