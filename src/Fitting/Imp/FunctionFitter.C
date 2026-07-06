// File: Fitting/Imp/FunctionFitter.C  Factory for the fitter faces -- the only place a concrete impl is named.
module;
#include <memory>
#include <cassert>
module qchem.Fitting.FunctionFitter;
import qchem.Fitting.Internal.FunctionFitterImp;   // FunctionFitterImp (Scalar) + IntegralConstrainedFF (Density)
import qchem.Fitting.Internal.OrthoFunctionFitter; // OrthoFunctionFitter (the orthonormal G-space density fit)

namespace qchem::Fitting
{

std::unique_ptr<FunctionFitter_Scalar<double>>
MakeScalarFitter(std::shared_ptr<const BasisSet::rFIT_SF_ABS>& bs)
{
    // The Gaussian potential fit needs the overlap metric-solve face; recover it from the neutral fit-basis
    // handle (a sanctioned abstract->abstract cross-cast), mirroring MakeDensityFitter.  An orthonormal (PW)
    // fit basis takes the reciprocal-space scalar fitter instead (a separate dcmplx route).
    auto nonOrtho = std::dynamic_pointer_cast<const BasisSet::FIT_SF_NonOrtho>(bs);
    assert(nonOrtho && "MakeScalarFitter: a Gaussian potential fit requires a FIT_SF_NonOrtho fit basis");
    return std::make_unique<FunctionFitterImp<double>>(nonOrtho);
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

std::unique_ptr<FunctionFitter_Density<dcmplx>>
MakeDensityFitter(std::shared_ptr<const BasisSet::cFIT_CD_ABS>& bs)
{
    // An orthonormal (plane-wave) {G} fit basis: the projection IS the fit, so the minimal core fitter
    // (no metric-solve refinement).  It holds the fit basis but delegates Repulsion to the orbital basis.
    return std::make_unique<OrthoFunctionFitter>(bs);
}

} //namespace
