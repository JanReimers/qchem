// File: FittedFunctionImp.C  Common imp for Fitted Functions.
module;
#include <iostream>
#include <cassert>
#include <vector>
#include "blaze/Math.h"
module qchem.FittedFunctionImp;
import qchem.FittedFunction;
import qchem.Fit_IBS;
import qchem.Orbital_DFT_IBS;
import qchem.Mesh;
import qchem.Streamable;
import qchem.Conversions;
import oml;

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
    itsFitCoeff=itsBasisSet->InvOverlap() * itsBasisSet->Overlap(itsMesh.get(),*ffc.GetScalarFunction());
}

template <class T> void FittedFunctionImp<T>::DoFitInternal(const DensityFFClient& ffc,double constraint)
{
    itsFitCoeff=itsBasisSet->InvRepulsion() * convert(ffc.GetRepulsion3C(itsBasisSet.get()));
}

//---------------------------------------------------------------------------
//
//  Provide Overlap and Repulsion matricies for derived classes.
//
template <class T> vec_t<T> FittedFunctionImp<T>::
FitGet2CenterOverlap(const Fit_IBS* bs) const
{
    return trans(itsBasisSet->Overlap(itsMesh.get(),*bs))*itsFitCoeff;
}

template <class T> vec_t<T> FittedFunctionImp<T>::
FitGet2CenterRepulsion(const Fit_IBS* bs) const
{
    return trans(itsBasisSet->Repulsion(*bs))*itsFitCoeff;
}

template <class T> smat_t<T> FittedFunctionImp<T>::
FitGet3CenterOverlap(const Orbital_DFT_IBS<double>* bs) const
{
    const ERI3<T>& O3=bs->Overlap3C(*itsBasisSet);
    smat_t<T> J=zero<T>(bs->GetNumFunctions());
    size_t i=0;
    for (auto c:itsFitCoeff) J+=c*O3[i++];
    assert(!isnan(J));
    return J;
}

//-------------------------------------------------------------------------------------
//
//  Get overlap and repulsion with a charge density.  And total charge.
//  Again only for derived classes.
//
template <class T> double FittedFunctionImp<T>::
FitGetOverlap(const FittedFunctionImp<T>* ffi) const
{
    return
        trans(itsFitCoeff) *
        itsBasisSet->Overlap(itsMesh.get(),*ffi->itsBasisSet) *
        ffi->itsFitCoeff;
}

template <class T> double FittedFunctionImp<T>::
FitGetRepulsion(const FittedFunctionImp<T>* ffi) const
{
    return
        trans(itsFitCoeff) * itsBasisSet->Repulsion(*ffi->itsBasisSet.get()) *
        ffi->itsFitCoeff;
}

template <class T> double FittedFunctionImp<T>::FitGetCharge() const
{
    return trans(itsFitCoeff) * itsBasisSet->Charge();
}

//------------------------------------------------------------------------
//
//  Handy utilities for fitted functions.
//
template <class T> void FittedFunctionImp<T>::FitMixIn(const FittedFunction& ff,double c)
{
    const FittedFunctionImp<T>* ffi = dynamic_cast<const FittedFunctionImp<T>*>(&ff);
    assert(ffi);
    assert(itsBasisSet->GetID() == ffi->itsBasisSet->GetID());
    itsFitCoeff = itsFitCoeff*(1-c) + ffi->itsFitCoeff*c;
}

template <class T> double FittedFunctionImp<T>::FitGetChangeFrom(const FittedFunction& ff) const
{
    const FittedFunctionImp<T>* ffi = dynamic_cast<const FittedFunctionImp<T>*>(&ff);
    assert(ffi);
    assert(itsBasisSet->GetID() == ffi->itsBasisSet->GetID());
    return max(abs(itsFitCoeff - ffi->itsFitCoeff));
}

template <class T> void FittedFunctionImp<T>::ReScale(double factor)
{
    itsFitCoeff*=factor;
}

//-------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> double  FittedFunctionImp<T>::operator()(const RVec3& r) const
{
    return trans(itsFitCoeff) * convert((*itsBasisSet)(r));
}

template <class T> RVec3  FittedFunctionImp<T>::Gradient(const RVec3& r) const
{
    Vector<RVec3> br = itsBasisSet->Gradient(r);
    RVec3 ret(0,0,0);
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
