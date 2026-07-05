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

template <class T, class Face, class FBS> void FitImpBase<T,Face,FBS>::FitMixIn(const Face& ff,double c)
{
    const FitImpBase* ffi = dynamic_cast<const FitImpBase*>(&ff);
    assert(ffi);
    assert(itsBasisSet->GetID() == ffi->itsBasisSet->GetID());
    itsFitCoeff = itsFitCoeff*(1-c) + ffi->itsFitCoeff*c;
}

template <class T, class Face, class FBS> double FitImpBase<T,Face,FBS>::FitGetChangeFrom(const Face& ff) const
{
    const FitImpBase* ffi = dynamic_cast<const FitImpBase*>(&ff);
    assert(ffi);
    assert(itsBasisSet->GetID() == ffi->itsBasisSet->GetID());
    return blazem::max(blazem::abs(itsFitCoeff - ffi->itsFitCoeff));
}

template <class T, class Face, class FBS> double FitImpBase<T,Face,FBS>::operator()(const rvec3_t& r) const
{
    return blazem::trans(itsFitCoeff) * (*itsBasisSet)(r);
}

template <class T, class Face, class FBS> rvec3_t FitImpBase<T,Face,FBS>::Gradient(const rvec3_t& r) const
{
    vec_t<rvec3_t> br = itsBasisSet->Gradient(r);
    rvec3_t ret(0,0,0);
    auto c(itsFitCoeff.begin());
    auto b(br.begin());
    for (; b!=br.end()&&c!=itsFitCoeff.end(); b++,c++) ret+=(*c) * (*b);
    return ret;
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
template <class T> void FunctionFitterImp<T>::DoFit(const ScalarFFClient& ffc)
{
    auto Sinv=this->itsBasisSet->InvOverlap();
    this->itsFitCoeff = Sinv * this->itsBasisSet->Overlap(*ffc.GetScalarFunction());
}

template <class T> hmat_t<T> FunctionFitterImp<T>::Overlap(const robs_t<T>* bs) const
{
    auto dftbs=dynamic_cast<const BasisSet::Orbital_DFT_IBS<T>*>(bs); // robs_t is the 1E base; need the 3-centre one
    assert(dftbs && "FunctionFitterImp::Overlap: Gaussian fitting needs an Orbital_DFT_IBS (3-centre) basis");
    const ERI3<T>& O3=dftbs->Overlap3C(*this->itsBasisSet);
    hmat_t<T> J=blazem::zeroH<T>(bs->GetNumFunctions());
    size_t i=0;
    for (auto c:this->itsFitCoeff) J+=c*O3[i++];
    assert(!blazem::isnan(J));
    return J;
}

// Both FitImpBase faces are emitted HERE, where the shared member definitions above are visible (the
// Density-face members would otherwise be undefined: ConstrainedFF.C can't emit what it can't see).
template class FitImpBase<double, FunctionFitter_Scalar <double>, BasisSet::FIT_SF_ABS>;
template class FitImpBase<double, FunctionFitter_Density<double>, BasisSet::FIT_CD_NonOrtho>;
template class FunctionFitterImp<double>;

} //namespace
