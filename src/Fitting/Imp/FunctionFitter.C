// File: Fitting/Imp/FunctionFitter.C  Factory for the fitter faces -- the only place a concrete impl is named.
module;
#include <memory>
module qchem.Fitting.FunctionFitter;
import qchem.Fitting.Internal.FunctionFitterImp;   // FunctionFitterImp (Scalar) + IntegralConstrainedFF (Density)

namespace qchem::Fitting
{

std::unique_ptr<FunctionFitter_Scalar<double>>
MakeScalarFitter(std::shared_ptr<const BasisSet::FIT_SF_ABS>& bs)
{
    return std::make_unique<FunctionFitterImp<double>>(bs);
}

std::unique_ptr<FunctionFitter_Density<double>>
MakeDensityFitter(std::shared_ptr<const BasisSet::FIT_CD_ABS>& bs)
{
    return std::make_unique<IntegralConstrainedFF<double>>(bs);
}

} //namespace
