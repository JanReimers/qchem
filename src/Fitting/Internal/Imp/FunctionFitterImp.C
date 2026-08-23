    // File: FunctionFitterImp.C  Shared coefficient machinery + the Scalar (overlap-metric) fitter.
module;
#include <iostream>
#include <cassert>
#include <vector>
module qchem.Fitting.Internal.FunctionFitterImp;
import qchem.Fitting.Types;

import qchem.Streamable;
import qchem.Blaze;

namespace qchem::Fitting
{

//=========================================================================== FitImpBase (shared)
//
//  Coefficient + real-space machinery common to both faces.  The metric-specific DoFit + contraction
//  live in the leaf impls; everything here touches only itsFitCoeff (+ the basis for real-space eval).
//
template <class T, class Face, class FBS> void FitImpBase<T,Face,FBS>::ReScale(double factor)
{
    itsFitCoeff*=factor;
}

// The fitted FIELD: ASK THE BASIS to evaluate its own expansion (doc/CleanupCandidates.md R1.0 step 1).
// This used to reach in for the basis's op(r)/Gradient and do the sum here, which is the pattern that let
// a basis promise pointwise values it could not honestly give -- and it was the LAST consumer of op(r) on
// a fit basis anywhere, so retiring it is what let VectorFunction<T> come off IrrepBasisSet<T>.
template <class T, class Face, class FBS> double FitImpBase<T,Face,FBS>::operator()(const rvec3_t& r) const
{
    return itsBasisSet->EvalField(itsFitCoeff, r);
}

template <class T, class Face, class FBS> rvec3_t FitImpBase<T,Face,FBS>::Gradient(const rvec3_t& r) const
{
    return itsBasisSet->EvalFieldGradient(itsFitCoeff, r);
}

template <class T, class Face, class FBS> std::ostream& FitImpBase<T,Face,FBS>::Write(std::ostream& os) const
{
    os << "Fit Function: " << std::endl;
    os << *itsBasisSet;
    os << std::endl;
    os << "  Coeff=" << itsFitCoeff << std::endl;
    return os;
}

//=========================================================================== FunctionFitterImp (Scalar)
//
//  Overlap-metric (potential) fit:  c = S^-1 <f_a|f>, contracted as Sum_a c_a <Oi|f_a|Oj>.
//
template <class T> void FunctionFitterImp<T>::DoFit(const ProjectedScalar_R& r)
{
    auto Sinv=this->itsBasisSet->InvOverlap();
    this->itsFitCoeff = Sinv * this->itsBasisSet->Overlap(*r.GetScalarFunction());
}

// ASK THE FIT BASIS for <i|f_a|j>, then contract MY coefficients into it (2026-08-23).  Same object and
// same cache as the old orb->Overlap3C(fitBasis): FIT_SF_ABS::Overlap3C's DEFAULT is exactly that call.
// Spelling it this way is what makes all three representations one shape -- a delta basis, which builds
// the tensor itself from op(r) rather than delegating, is reached by the identical line.
template <class T> hmat_t<T> FunctionFitterImp<T>::Overlap(const robs_t<T>* bs) const
{
    assert(bs);
    const Projector3<T>& O3=this->itsBasisSet->Overlap3C(*bs);
    hmat_t<T> J=blazem::zeroH<T>(bs->GetNumFunctions());
    size_t i=0;
    for (auto c:this->itsFitCoeff) J+=c*O3.dense[i++];
    assert(!blazem::isnan(J));
    return J;
}

// Both FitImpBase faces are emitted HERE, where the shared member definitions above are visible (the
// Density-face members would otherwise be undefined: ConstrainedFF.C can't emit what it can't see).
template class FitImpBase<double, FunctionFitter_Scalar, BasisSet::FIT_SF_NonOrtho>;
template class FitImpBase<double, FunctionFitter_Density_NonOrtho<double>, BasisSet::FIT_CD_NonOrtho>;
template class FunctionFitterImp<double>;

} //namespace
