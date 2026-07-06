// File: Fitting/Imp/FunctionFitter.C  Factory for the fitter faces + the ProjectedDensity_AO metric defaults.
module;
#include <memory>
#include <cassert>
module qchem.Fitting.FunctionFitter;
import qchem.Fitting.Internal.FunctionFitterImp;   // FunctionFitterImp (Scalar) + IntegralConstrainedFF (Density)
import qchem.Fitting.Internal.OrthoFunctionFitter; // OrthoFunctionFitter (the orthonormal G-space density fit)
import qchem.BasisSet.Fit_IBS;                     // FIT_CD_NonOrtho (the Coulomb metric-solve face)
import qchem.Blaze;                                // rsmat_t * rvec_t (the J^-1 solve)

namespace qchem::Fitting
{

//---------------------------------------------------------- ProjectedDensity_AO metric defaults
//
//  The default unconstrained-fit STRATEGY (a density with a matrix): the Coulomb-metric solve
//  c0 = J^-1 <rho|c>.  The seed overrides GetUnconstrainedFit (overlap metric) and never reaches here.
//
rvec_t ProjectedDensity_AO::GetUnconstrainedFit(const BasisSet::rFIT_CD_ABS* fbs) const
{
    // "I want more": broaden the neutral CD-fit face to its Coulomb metric-solve capability (the sanctioned
    // request pattern -- a real density matrix genuinely needs J^-1).
    auto* no = dynamic_cast<const BasisSet::FIT_CD_NonOrtho*>(fbs);
    assert(no && "ProjectedDensity_AO::GetUnconstrainedFit: the Coulomb-metric default needs a FIT_CD_NonOrtho fit basis");
    return no->InvRepulsion() * GetRepulsion3C(fbs);   // c0 = J^-1 <rho|c>
}

rvec_t ProjectedDensity_AO::GetRepulsion3C(const BasisSet::rFIT_CD_ABS*) const
{
    // Reached only if a projection uses the Coulomb-metric default WITHOUT providing a Coulomb RHS -- i.e. a
    // matrix-free projection that forgot to override GetUnconstrainedFit.  A seed overrides that and never lands here.
    assert(false && "ProjectedDensity_AO::GetRepulsion3C: this projection has no Coulomb RHS (a matrix-free "
                    "seed must override GetUnconstrainedFit with its own metric instead)");
    return rvec_t{};
}

std::unique_ptr<FunctionFitter_Scalar<double>>
Factory(std::shared_ptr<const BasisSet::rFIT_SF_ABS>& bs)
{
    // A real (Gaussian/Slater/BSpline) potential fit needs the overlap metric-solve face.  isOrtho()==false is
    // the basis's CONTRACT that it carries that face, so the down-cast is GUARANTEED (not an assumption); the
    // second assert just belt-and-suspenders the contract.  (An orthonormal PW basis is dcmplx -> other overload.)
    assert(!bs->isOrtho() && "Fitting::Factory(rFIT_SF_ABS): a real potential-fit basis must be non-orthonormal");
    auto nonOrtho = std::dynamic_pointer_cast<const BasisSet::FIT_SF_NonOrtho>(bs);
    assert(nonOrtho && "isOrtho()==false contract broken: a real potential-fit basis must IS-A FIT_SF_NonOrtho");
    return std::make_unique<FunctionFitterImp<double>>(nonOrtho);
}

std::unique_ptr<FunctionFitter_Scalar<dcmplx>>
Factory(std::shared_ptr<const BasisSet::cFIT_SF_ABS>& bs)
{
    // An orthonormal (plane-wave) {G} scalar fit basis: the projection IS the fit, so the minimal core
    // scalar fitter (no metric solve).  It holds the fit basis but delegates the assembly to the orbital basis.
    assert(bs->isOrtho() && "Fitting::Factory(cFIT_SF_ABS): a plane-wave potential-fit basis must be orthonormal");
    return std::make_unique<OrthoScalarFitter>(bs);
}

std::unique_ptr<FunctionFitter_Density_NonOrtho<double>>
Factory(std::shared_ptr<const BasisSet::rFIT_CD_ABS>& bs)
{
    // A real (Gaussian/Slater/BSpline) density fit needs the Coulomb metric-solve face.  isOrtho()==false is the
    // basis's CONTRACT that it carries that face, so the down-cast is GUARANTEED.  (An orthonormal PW density
    // fit basis is dcmplx -> the other overload.)
    assert(!bs->isOrtho() && "Fitting::Factory(rFIT_CD_ABS): a real density-fit basis must be non-orthonormal");
    auto nonOrtho = std::dynamic_pointer_cast<const BasisSet::FIT_CD_NonOrtho>(bs);
    assert(nonOrtho && "isOrtho()==false contract broken: a real density-fit basis must IS-A FIT_CD_NonOrtho");
    return std::make_unique<IntegralConstrainedFF<double>>(nonOrtho);
}

std::unique_ptr<FunctionFitter_Density<dcmplx>>
Factory(std::shared_ptr<const BasisSet::cFIT_CD_ABS>& bs)
{
    // An orthonormal (plane-wave) {G} fit basis: the projection IS the fit, so the minimal core fitter
    // (no metric-solve refinement).  It holds the fit basis but delegates Repulsion to the orbital basis.
    assert(bs->isOrtho() && "Fitting::Factory(cFIT_CD_ABS): a plane-wave density-fit basis must be orthonormal");
    return std::make_unique<OrthoFunctionFitter>(bs);
}

} //namespace
