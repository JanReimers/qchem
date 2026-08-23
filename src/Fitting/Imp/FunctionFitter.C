// File: Fitting/Imp/FunctionFitter.C  Factory for the fitter faces + the ProjectedDensity_AO metric defaults.
module;
#include <memory>
#include <cassert>
module qchem.Fitting.FunctionFitter;
import qchem.Fitting.Internal.FunctionFitterImp;   // FunctionFitterImp (Scalar) + IntegralConstrainedFF (Density)
import qchem.Fitting.Internal.OrthoFunctionFitter;       // OrthoFunctionFitter (the orthonormal G-space density fit)
import qchem.Fitting.Internal.OrthoNormalFunctionFitter; // OrthoNormalScalarFitter (the orthoNORMAL potential fit)
import qchem.BasisSet.Fit_IBS;                     // FIT_CD_NonOrtho (the Coulomb metric-solve face)
import qchem.BasisSet.G_FieldEvaluator;            // G_FieldEvaluator / G_RasterTransform (the raster half of the contract)
import qchem.BasisSet.Quadrature;                  // BasisSet::Quadrature (the points+weights half)
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

std::unique_ptr<FunctionFitter_Scalar>
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

std::unique_ptr<FunctionFitter_Scalar>
Factory(std::shared_ptr<const BasisSet::cFIT_SF_ABS>& bs)
{
    // An orthonormal (plane-wave) {G} scalar fit basis: the projection IS the fit, so no metric solve.
    //
    // BOTH halves of the contract are checked HERE, at the one seam that builds the object.  They are
    // INDEPENDENT capabilities and neither implies the other:
    //   isOrtho()         -- the METRIC axis: DIAGONAL, so no linear solve.  Genuinely general; an
    //                        orthogonal wavelet basis, or a delta basis (metric diag(w_g)), satisfies it
    //                        too -- the projection IS the fit up to that known diagonal.
    //   Quadrature        -- the POINTS+WEIGHTS axis: where the fit is sampled and E_xc integrated.  This
    //                        one IS general -- any point set answers it (that is the whole point of the
    //                        2026-08-22 split), and a fit basis without it is not defined at all.
    //   G_RasterTransform -- the TRANSFORM axis: the FFT pair the projection is batched through.  NOT
    //                        general -- that face is reciprocal-space BY INTERFACE (ΔG_Map, ForwardFFT,
    //                        GridCoeff keyed by an integer ivec3_t reciprocal-index difference), so only a
    //                        raster-backed G-space basis can implement it.  THIS fitter needs it (its
    //                        DoFit batch-projects by FFT), hence the check; a δ-basis fitter would not.
    //                        (The fitter's op(r) additionally consumes the G_FieldEvaluator evaluate face --
    //                        checked with it, same seam, since the concrete engine carries both.)
    // The grid check used to live in OrthoScalarFitter::FitGrid, i.e. at FIRST GRID USE: an ortho-but-not-
    // G-space fit basis constructed happily and tripped later, somewhere else.  Two-phase contract, same
    // smell as R2.10's SetMesh -- so it is established once, where the object is made.
    assert(bs->isOrtho() && "Fitting::Factory(cFIT_SF_ABS): a plane-wave potential-fit basis must be orthogonal");
    assert(dynamic_cast<const BasisSet::Quadrature*>(bs.get())
           && "Fitting::Factory(cFIT_SF_ABS): the fit basis must CARRY ITS QUADRATURE (BasisSet::Quadrature) "
              "-- a fit basis is a family of weight vectors over shared points, so it is not defined without them");
    assert(dynamic_cast<const BasisSet::G_RasterTransform*>(bs.get())
           && "Fitting::Factory(cFIT_SF_ABS): the fit basis must ALSO provide the G_RasterTransform FFT "
              "pair -- this fitter batch-projects by forward FFT, which only a raster can do");
    assert(dynamic_cast<const BasisSet::G_FieldEvaluator*>(bs.get())
           && "Fitting::Factory(cFIT_SF_ABS): the fit basis must provide G_FieldEvaluator so the fitted "
              "v_xc,fit(r) stays evaluatable (GUI / fit-residual)");
    return std::make_unique<OrthoNormalScalarFitter>(bs);
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
    assert(bs->isOrtho() && "Fitting::Factory(cFIT_CD_ABS): a plane-wave density-fit basis must be orthogonal");
    return std::make_unique<OrthoFunctionFitter>(bs);
}

} //namespace
