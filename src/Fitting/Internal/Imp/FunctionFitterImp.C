    // File: FunctionFitterImp.C  Common imp for Fitted Functions.
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

//---------------------------------------------------------------------
//
//  Construction zone.  The CDFit flag indicates to fit the electric
//  instead of the charge density.  The only difference in practice
//  is that all overlap integrals are replaced with repulsion integrals.
//
template <class T> FunctionFitterImp<T>::
FunctionFitterImp(bs_t& fbs)
    : itsBasisSet(fbs)
    , itsFitCoeff(fbs->GetNumFunctions(),0.0)
{
    itsFitCoeff[0]=1.0/itsBasisSet->Charge()[0]; //Wild guess with the correct total charge.
};

template <class T> FunctionFitterImp<T>::FunctionFitterImp()
    : itsBasisSet (    )
    , itsFitCoeff (    )
{

};

template <class T> FunctionFitterImp<T>::~FunctionFitterImp()
{
}


//--------------------------------------------------------------------------
//
//  Implement all DoFit functions.  The overlaps will be accumulated in
//  itsFitCoeff by the call to GetRepulsions or GetOverlap.
//
template <class T> void FunctionFitterImp<T>::DoFit(const ScalarFFClient& ffc)
{
    DoFitInternal(ffc,0); //No contraint.
}
template <class T> void FunctionFitterImp<T>::DoFit(const DensityFFClient& ffc)
{
    DoFitInternal(ffc,0); //No contraint.
}
template <class T> void FunctionFitterImp<T>::DoFit(const FourierMap&)
{
    assert(false && "FunctionFitterImp::DoFit(FourierMap): the Gaussian fitter fits via a client callback, "
                    "not pre-computed Fourier coefficients (that is the plane-wave path).");
}

template <class T> void FunctionFitterImp<T>::DoFitInternal(const ScalarFFClient& ffc,double constraint)
{
    auto Sinv=itsBasisSet->InvOverlap();
    itsFitCoeff= Sinv * itsBasisSet->Overlap(*ffc.GetScalarFunction());
}

template <class T> void FunctionFitterImp<T>::DoFitInternal(const DensityFFClient& ffc,double constraint)
{   
    auto Sinv=itsBasisSet->InvRepulsion();
    itsFitCoeff=Sinv * ffc.GetRepulsion3C(itsBasisSet.get());
}

//---------------------------------------------------------------------------
//
//  Fit-derived quantities the clients query (the "what's your overlap/repulsion with this basis?" side).
//
template <class T> hmat_t<T> FunctionFitterImp<T>::
Overlap(const obs_t<T>* bs) const
{
    auto dftbs=dynamic_cast<const BasisSet::Orbital_DFT_IBS<T>*>(bs); // obs_t is the 1E base; need the 3-centre one
    assert(dftbs && "FunctionFitterImp::Overlap: Gaussian fitting needs an Orbital_DFT_IBS (3-centre) basis");
    const ERI3<T>& O3=dftbs->Overlap3C(*itsBasisSet);
    hmat_t<T> J=blazem::zeroH<T>(bs->GetNumFunctions());
    size_t i=0;
    for (auto c:itsFitCoeff) J+=c*O3[i++];
    assert(!blazem::isnan(J));
    return J;
}

template <class T> hmat_t<T> FunctionFitterImp<T>::
Repulsion(const obs_t<T>* bs) const
{
    auto dftbs=dynamic_cast<const BasisSet::Orbital_DFT_IBS<T>*>(bs); // obs_t is the 1E base; need the 3-centre one
    assert(dftbs && "FunctionFitterImp::Repulsion: Gaussian fitting needs an Orbital_DFT_IBS (3-centre) basis");
    const ERI3<T>& R3=dftbs->Repulsion3C(*itsBasisSet);
    hmat_t<T> J=blazem::zeroH<T>(bs->GetNumFunctions());
    size_t i=0;
    for (auto c:itsFitCoeff) J+=c*R3[i++];
    assert(!blazem::isnan(J));
    return J;
}

template <class T> double FunctionFitterImp<T>::
FitGetRepulsion(const FunctionFitterImp<T>* ffi) const
{
    return
        blazem::trans(itsFitCoeff) * itsBasisSet->Repulsion(*ffi->itsBasisSet.get()) *
        ffi->itsFitCoeff;
}

template <class T> double FunctionFitterImp<T>::FitGetSelfRepulsion() const
{
    return FitGetRepulsion(this);   // <fit|1/r12|fit>
}

template <class T> double FunctionFitterImp<T>::Integral() const
{
    return blazem::trans(itsFitCoeff) * itsBasisSet->Charge();
}

//------------------------------------------------------------------------
//
//  Handy utilities for fitted functions.
//
template <class T> void FunctionFitterImp<T>::FitMixIn(const FunctionFitter<T>& ff,double c)
{
    const FunctionFitterImp<T>* ffi = dynamic_cast<const FunctionFitterImp<T>*>(&ff);
    assert(ffi);
    assert(itsBasisSet->GetID() == ffi->itsBasisSet->GetID());
    itsFitCoeff = itsFitCoeff*(1-c) + ffi->itsFitCoeff*c;
}

template <class T> double FunctionFitterImp<T>::FitGetChangeFrom(const FunctionFitter<T>& ff) const
{
    const FunctionFitterImp<T>* ffi = dynamic_cast<const FunctionFitterImp<T>*>(&ff);
    assert(ffi);
    assert(itsBasisSet->GetID() == ffi->itsBasisSet->GetID());
    return blazem::max(blazem::abs(itsFitCoeff - ffi->itsFitCoeff));
}

template <class T> void FunctionFitterImp<T>::ReScale(double factor)
{
    itsFitCoeff*=factor;
}

//-------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> double  FunctionFitterImp<T>::operator()(const rvec3_t& r) const
{
    return blazem::trans(itsFitCoeff) * (*itsBasisSet)(r);
}

template <class T> rvec3_t  FunctionFitterImp<T>::Gradient(const rvec3_t& r) const
{
    vec_t<rvec3_t> br = itsBasisSet->Gradient(r);
    rvec3_t ret(0,0,0);
    auto c(itsFitCoeff.begin());
    auto b(br.begin());
    for (; b!=br.end()&&c!=itsFitCoeff.end(); b++,c++) ret+=(*c) * (*b);
    return ret;
}

//-----------------------------------------------------------------------
//
//  Streamable stuff.
//
template <class T> std::ostream& FunctionFitterImp<T>::Write(std::ostream& os) const
{
    os << "Fit Function: " << std::endl;

    os << *itsBasisSet;
    os << std::endl;
    os << "  Coeff=" << itsFitCoeff << std::endl;

    return os;
}


template class FunctionFitterImp<double>;

} //namespace