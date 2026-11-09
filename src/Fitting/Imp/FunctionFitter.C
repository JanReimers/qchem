// File: Fitting/Imp/FunctionFitter.C  Factory for the fitter faces.
module;
#include <memory>
#include <cassert>
#include <stdexcept>
#include <type_traits>
module qchem.Fitting.FunctionFitter;
import qchem.Fitting.Internal.FunctionFitterImp;   // FunctionFitterImp (the Gaussian scalar fit)
import qchem.Fitting.Internal.OrthoNormalFunctionFitter; // OrthoNormalScalarFitter (the orthoNORMAL potential fit)
import qchem.Fitting.Internal.DeltaFunctionFitter;       // DeltaScalarFitter (the delta/identity potential fit)
import qchem.BasisSet.Orbital_DFT_IBS;                     // FIT_CD_NonOrtho (the Coulomb metric-solve face)
import qchem.BasisSet.G_FieldEvaluator;            // G_FieldEvaluator / G_RasterTransform (the raster half of the contract)

namespace qchem::Fitting
{

//  (2026-09-19) The Coulomb metric solve c0 = J^-1 <rho|c> no longer lives on the projection: it is the
//  FITTER's -- ConstrainedFF::DoFitUnconstrained applies its own metric to the RHS the density hands it
//  through the DensityProjector face (R1.0q: the fitter owns the metric, the density owns D, one integrator
//  object owns the tensor).  This TU keeps only the factories.

// THE SCALAR FITTER, ONE FACTORY (2026-10-31).  There were two overloads, on rFIT_SF_ABS and cFIT_SF_ABS, and the real one
// ASSERTED !isOrtho and down-cast to FIT_SF_NonOrtho: the scalar T stood in for the metric axis, so a real δ basis (diagonal
// metric) could not be fitted at all.  The complex overload's own rule -- CAPABILITY decides -- now serves both axes.  The
// three capabilities are INDEPENDENT axes of a fit basis, asked in order of specificity:
//   G_RasterTransform -- the TRANSFORM axis: the FFT pair the projection is batched through.  Reciprocal-space BY INTERFACE
//                        (ΔG_Map, ForwardFFT, GridCoeff keyed by an integer reciprocal-index difference), so only a raster-backed
//                        G-space basis implements it, and its fitter is complex.  (The fitter's op(r) additionally consumes the
//                        G_FieldEvaluator evaluate face -- checked with it, same seam, since the concrete engine carries both.)
//   isOrtho()         -- the METRIC axis: DIAGONAL, so no linear solve -- the projection IS the fit up to that known diagonal.
//                        A δ basis (metric diag(w_g)) on either axis; any orthogonal basis without a raster.
//   FIT_SF_NonOrtho   -- the overlap-metric SOLVE face (S, S^-1, its norms): the Gaussian/Slater/BSpline auxiliary bases.
// The checks used to live at FIRST USE (OrthoScalarFitter::FitGrid): a basis constructed happily and tripped later, somewhere
// else -- so every contract is established here, once, where the object is made.
template <class T> std::unique_ptr<ScalarFit<T>>
Factory(std::shared_ptr<const BasisSet::FIT_SF_ABS<T>>& bs)
{
    if (dynamic_cast<const BasisSet::G_RasterTransform*>(bs.get()))
    {
        assert(bs->isOrtho() && "Fitting::Factory: a raster fit basis is orthonormal by construction");
        assert(dynamic_cast<const BasisSet::G_FieldEvaluator*>(bs.get())
               && "Fitting::Factory: a raster fit basis must provide G_FieldEvaluator so the "
                  "fitted v_xc,fit(r) stays evaluatable (GUI / fit-residual)");
        // The raster fitter is reciprocal-space -- complex -- BY ITS FACE, not by a choice made here: a REAL basis claiming a
        // raster transform has no fitter to go to, and that is a design defect to name, not an input to refuse quietly.
        if constexpr (std::is_same_v<T,dcmplx>) return std::make_unique<OrthoNormalScalarFitter>(bs);
        else throw std::logic_error("Fitting::Factory: a REAL fit basis offers G_RasterTransform, whose fitter is reciprocal-space "
                                    "(complex) by its face -- there is no real raster fitter (an OOD gap, pin 30)");
    }
    if (bs->isOrtho()) return std::make_unique<DeltaScalarFitter<T>>(bs);   // no raster: projection / metric diagonal
    if constexpr (std::is_same_v<T,double>)
        if (auto nonOrtho=std::dynamic_pointer_cast<const BasisSet::FIT_SF_NonOrtho>(bs))
            return std::make_unique<FunctionFitterImp<double>>(nonOrtho);
    throw std::logic_error("Fitting::Factory: the scalar fit basis '"+bs->Name()+"' is neither orthogonal (isOrtho) nor offers the "
                           "overlap-metric solve face (FIT_SF_NonOrtho) -- no fitter can use it (an OOD gap, pin 30)");
}
template std::unique_ptr<ScalarFit<double>> Factory<double>(std::shared_ptr<const BasisSet::FIT_SF_ABS<double>>&);
template std::unique_ptr<ScalarFit<dcmplx>> Factory<dcmplx>(std::shared_ptr<const BasisSet::FIT_SF_ABS<dcmplx>>&);

template <class T> std::unique_ptr<DeltaScalarFit<T>>
MakeDeltaScalarFitter(std::shared_ptr<const BasisSet::FIT_SF_ABS<T>>& bs)
{
    if (!bs->isOrtho() || dynamic_cast<const BasisSet::G_RasterTransform*>(bs.get()))
        throw std::logic_error("Fitting::MakeDeltaScalarFitter: the fit basis '"+bs->Name()+"' is not a delta basis (orthogonal, no raster)");
    return std::make_unique<DeltaScalarFitter<T>>(bs);
}
template std::unique_ptr<DeltaScalarFit<double>> MakeDeltaScalarFitter<double>(std::shared_ptr<const BasisSet::FIT_SF_ABS<double>>&);
template std::unique_ptr<DeltaScalarFit<dcmplx>> MakeDeltaScalarFitter<dcmplx>(std::shared_ptr<const BasisSet::FIT_SF_ABS<dcmplx>>&);

} //namespace
