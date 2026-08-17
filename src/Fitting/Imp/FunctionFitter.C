// File: Fitting/Imp/FunctionFitter.C  Factory for the fitter faces + the ProjectedDensity_AO metric defaults.
module;
#include <memory>
#include <cassert>
module qchem.Fitting.FunctionFitter;
import qchem.Fitting.Internal.FunctionFitterImp;   // FunctionFitterImp (Scalar) + IntegralConstrainedFF (Density)
import qchem.Fitting.Internal.OrthoFunctionFitter; // OrthoFunctionFitter (the orthonormal G-space density fit)
import qchem.BasisSet.Fit_IBS;                     // FIT_CD_NonOrtho (the Coulomb metric-solve face)
import qchem.BasisSet.G_FieldEvaluator;            // G_FieldEvaluator (the FFT grid half of the ortho contract)
import qchem.Blaze;                                // rsmat_t * rvec_t (the J^-1 solve)

namespace qchem::Fitting
{

//---------------------------------------------------------- the Coulomb metric solve (V1.16)
//
//  The ONLY shared body left.  There is no longer a metric DEFAULT on the base and no poisoned sibling:
//  a projection has the Coulomb face or the overlap face, and each face's own method is pure virtual.
//
rvec_t CoulombMetric_ProjectedDensity::GetUnconstrainedFit(const BasisSet::rFIT_CD_ABS* fbs) const
{
    // "I want more": broaden the neutral CD-fit face to its Coulomb metric-solve capability (the sanctioned
    // request pattern -- a real density matrix genuinely needs J^-1).
    auto* no = dynamic_cast<const BasisSet::FIT_CD_NonOrtho*>(fbs);
    assert(no && "CoulombMetric_ProjectedDensity: the Coulomb-metric solve needs a FIT_CD_NonOrtho fit basis");
    return no->InvRepulsion() * GetRepulsion3C(fbs);   // c0 = J^-1 <rho|c>
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

std::unique_ptr<GriddedScalarFitter>
Factory(std::shared_ptr<const BasisSet::cFIT_SF_ABS>& bs)
{
    // An orthonormal (plane-wave) {G} scalar fit basis: the projection IS the fit, so no metric solve.  Returns
    // the GriddedScalarFitter face -- it owns the FFT quadrature grid (from bs) and exposes it, so the XC term
    // borrows ONE grid instead of cross-casting bs a second time.
    //
    // BOTH halves of the contract are checked HERE, at the one seam that builds the object.  They are
    // INDEPENDENT capabilities and neither implies the other:
    //   isOrtho()         -- the METRIC axis: the projection IS the fit.  Genuinely general; an orthonormal
    //                        WAVELET basis would satisfy it too.
    //   G_Quadrature      -- the QUADRATURE axis: the FFT raster the fit is sampled on and E_xc integrated
    //                        on.  NOT general -- that face is reciprocal-space BY INTERFACE (ΔG_Map,
    //                        ForwardFFT, GridCoeff keyed by an integer ivec3_t reciprocal-index difference),
    //                        so only a G-space basis can implement it.  This is the BINDING requirement.
    //                        (The fitter's op(r) additionally consumes the G_FieldEvaluator evaluate face --
    //                        checked with it, same seam, since the concrete engine carries both.)
    // The grid check used to live in OrthoScalarFitter::FitGrid, i.e. at FIRST GRID USE: an ortho-but-not-
    // G-space fit basis constructed happily and tripped later, somewhere else.  Two-phase contract, same
    // smell as R2.10's SetMesh -- so it is established once, where the object is made.
    assert(bs->isOrtho() && "Fitting::Factory(cFIT_SF_ABS): a plane-wave potential-fit basis must be orthonormal");
    assert(dynamic_cast<const BasisSet::G_Quadrature*>(bs.get())
           && "Fitting::Factory(cFIT_SF_ABS): the fit basis must ALSO provide the G_Quadrature FFT grid "
              "engine -- orthonormality alone is not enough, the grid IS the XC quadrature");
    assert(dynamic_cast<const BasisSet::G_FieldEvaluator*>(bs.get())
           && "Fitting::Factory(cFIT_SF_ABS): the fit basis must provide G_FieldEvaluator so the fitted "
              "v_xc,fit(r) stays evaluatable (GUI / fit-residual)");
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
