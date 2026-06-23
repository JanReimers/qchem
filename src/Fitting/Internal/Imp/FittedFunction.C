    // File: FittedFunctionImp.C  Common imp for Fitted Functions.
module;
#include <iostream>
#include <cassert>
#include <vector>
module qchem.Fitting.Internal.FittedFunctionImp;
import qchem.Fitting.Types;

import qchem.Mesh;
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
template <class T> FittedFunctionImp<T>::
FittedFunctionImp(bs_t& fbs,mesh_t& m)
    : itsBasisSet(fbs)
    , itsFitCoeff(fbs->GetNumFunctions(),0.0)
    , itsMesh    (m)
{
    assert(itsMesh);
    itsFitCoeff[0]=1.0/itsBasisSet->Charge()[0]; //Wild guess with the correct total charge.
};

template <class T> FittedFunctionImp<T>::FittedFunctionImp()
    : itsBasisSet (    )
    , itsFitCoeff (    )
    , itsMesh     (0   )
{
    
};

template <class T> FittedFunctionImp<T>::~FittedFunctionImp()
{
}


//--------------------------------------------------------------------------
//
//  Implement all DoFit functions.  The overlaps will be accumulated in
//  itsFitCoeff by the call to GetRepulsions or GetOverlap.
//
template <class T> void FittedFunctionImp<T>::DoFit(const ScalarFFClient& ffc)
{
    DoFitInternal(ffc,0); //No contraint.
}
template <class T> void FittedFunctionImp<T>::DoFit(const DensityFFClient& ffc)
{
    DoFitInternal(ffc,0); //No contraint.
}

template <class T> void FittedFunctionImp<T>::DoFitInternal(const ScalarFFClient& ffc,double constraint)
{
    auto Sinv=itsBasisSet->InvOverlap();
    itsFitCoeff= Sinv * itsBasisSet->Overlap(itsMesh.get(),*ffc.GetScalarFunction());
}

template <class T> void FittedFunctionImp<T>::DoFitInternal(const DensityFFClient& ffc,double constraint)
{   
    auto Sinv=itsBasisSet->InvRepulsion();
    itsFitCoeff=Sinv * ffc.GetRepulsion3C(itsBasisSet.get());
}

//---------------------------------------------------------------------------
//
//  Fit-derived quantities the clients query (the "what's your overlap/repulsion with this basis?" side).
//
template <class T> smat_t<T> FittedFunctionImp<T>::
FitGet3CenterOverlap(const obs_t<T>* bs) const
{
    const ERI3<T>& O3=bs->Overlap3C(*itsBasisSet);
    smat_t<T> J=blazem::zero<T>(bs->GetNumFunctions());
    size_t i=0;
    for (auto c:itsFitCoeff) J+=c*O3[i++];
    assert(!blazem::isnan(J));
    return J;
}

template <class T> smat_t<T> FittedFunctionImp<T>::
FitGet3CenterRepulsion(const obs_t<T>* bs) const
{
    const ERI3<T>& R3=bs->Repulsion3C(*itsBasisSet);
    smat_t<T> J=blazem::zero<T>(bs->GetNumFunctions());
    size_t i=0;
    for (auto c:itsFitCoeff) J+=c*R3[i++];
    assert(!blazem::isnan(J));
    return J;
}

template <class T> double FittedFunctionImp<T>::
FitGetRepulsion(const FittedFunctionImp<T>* ffi) const
{
    return
        blazem::trans(itsFitCoeff) * itsBasisSet->Repulsion(*ffi->itsBasisSet.get()) *
        ffi->itsFitCoeff;
}

template <class T> double FittedFunctionImp<T>::FitGetSelfRepulsion() const
{
    return FitGetRepulsion(this);   // <fit|1/r12|fit>
}

template <class T> double FittedFunctionImp<T>::FitGetCharge() const
{
    return blazem::trans(itsFitCoeff) * itsBasisSet->Charge();
}

//------------------------------------------------------------------------
//
//  Handy utilities for fitted functions.
//
template <class T> void FittedFunctionImp<T>::FitMixIn(const FunctionFitter<T>& ff,double c)
{
    const FittedFunctionImp<T>* ffi = dynamic_cast<const FittedFunctionImp<T>*>(&ff);
    assert(ffi);
    assert(itsBasisSet->GetID() == ffi->itsBasisSet->GetID());
    itsFitCoeff = itsFitCoeff*(1-c) + ffi->itsFitCoeff*c;
}

template <class T> double FittedFunctionImp<T>::FitGetChangeFrom(const FunctionFitter<T>& ff) const
{
    const FittedFunctionImp<T>* ffi = dynamic_cast<const FittedFunctionImp<T>*>(&ff);
    assert(ffi);
    assert(itsBasisSet->GetID() == ffi->itsBasisSet->GetID());
    return blazem::max(blazem::abs(itsFitCoeff - ffi->itsFitCoeff));
}

template <class T> void FittedFunctionImp<T>::ReScale(double factor)
{
    itsFitCoeff*=factor;
}

//-------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> double  FittedFunctionImp<T>::operator()(const rvec3_t& r) const
{
    return blazem::trans(itsFitCoeff) * (*itsBasisSet)(r);
}

template <class T> rvec3_t  FittedFunctionImp<T>::Gradient(const rvec3_t& r) const
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
template <class T> std::ostream& FittedFunctionImp<T>::Write(std::ostream& os) const
{
    os << "Fit Function: " << std::endl;

    os << *itsBasisSet;
    os << std::endl;
    os << "  Coeff=" << itsFitCoeff << std::endl;

    return os;
}


template class FittedFunctionImp<double>;

} //namespace