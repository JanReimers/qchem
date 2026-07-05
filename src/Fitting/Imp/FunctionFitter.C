// File: Fitting/Imp/FunctionFitter.C  Factory for the fitter faces -- the only place a concrete impl is named.
module;
#include <memory>
#include <cassert>
module qchem.Fitting.FunctionFitter;
import qchem.Fitting.Internal.FunctionFitterImp;   // FunctionFitterImp (Scalar) + IntegralConstrainedFF (Density)

namespace qchem::Fitting
{

std::unique_ptr<FunctionFitter_Scalar<double>>
MakeScalarFitter(std::shared_ptr<const BasisSet::FIT_SF_ABS>& bs)
{
    return std::make_unique<FunctionFitterImp<double>>(bs);
}

std::unique_ptr<FunctionFitter_Density_NonOrtho<double>>
MakeDensityFitter(std::shared_ptr<const BasisSet::rFIT_CD_ABS>& bs)
{
    // The non-ortho (Gaussian) density fit needs the Coulomb metric-solve face; recover it from the neutral
    // fit-basis handle (a sanctioned abstract->abstract cross-cast).  An orthonormal (PW) fit basis takes the
    // reciprocal-space fitter instead (a separate dcmplx route), so for the double path this always holds.
    auto nonOrtho = std::dynamic_pointer_cast<const BasisSet::FIT_CD_NonOrtho>(bs);
    assert(nonOrtho && "MakeDensityFitter: a Gaussian density fit requires a FIT_CD_NonOrtho fit basis");
    return std::make_unique<IntegralConstrainedFF<double>>(nonOrtho);
}

} //namespace
