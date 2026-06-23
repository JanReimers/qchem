// File: Fitting/Imp/FunctionFitter.C  Factory for FunctionFitter -- the only place the concrete impl is named.
module;
#include <memory>
module qchem.Fitting.FunctionFitter;
import qchem.Fitting.Internal.FunctionFitterImp;   // FunctionFitterImp + IntegralConstrainedFF (the concrete fitters)

namespace qchem::Fitting
{

std::unique_ptr<FunctionFitter<double>>
MakeFunctionFitter(FitFlavour flavour, std::shared_ptr<const fbs_t>& bs, std::shared_ptr<const Mesh>& m)
{
    switch (flavour)
    {
        case FitFlavour::Unconstrained:     return std::make_unique<FunctionFitterImp<double>>     (bs,m);
        case FitFlavour::ChargeConstrained: return std::make_unique<IntegralConstrainedFF<double>>(bs,m);
    }
    return nullptr;
}

} //namespace
