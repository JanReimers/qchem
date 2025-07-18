// File: FittedFunctionImp.C  Common imp for Fitted Functions.

#include <iostream>
#include <cassert>
#include <vector>
#include "FittedFunction.H"
#include <ChargeDensity/ChargeDensity.H>
#include <BasisSet/Fit_IBS.H>
#include <BasisSet/DFT_IBS.H>

import Mesh;
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
    , itsFitCoeff(fbs->GetNumFunctions())
    , itsMesh    (m)
    , itsLAParams({qchem::Lapack,qchem::SVD,1e-10,1e-12})
{
    assert(itsMesh);
    Fill(itsFitCoeff,0.0);
    itsFitCoeff(1)=1.0/itsBasisSet->Charge()(1); //Wild guess with the correct total charge.
};

template <class T> FittedFunctionImp<T>::FittedFunctionImp()
    : itsBasisSet (    )
    , itsFitCoeff (    )
    , itsMesh     (0   )
{
    Fill(itsFitCoeff,0.0);
};

template <class T> FittedFunctionImp<T>::~FittedFunctionImp()
{
}


//--------------------------------------------------------------------------
//
//  Implement all DoFit functions.  The overlaps will be accumulated in
//  itsFitCoeff by the call to GetRepulsions or GetOverlap.
//
template <class T> double FittedFunctionImp<T>::DoFit(const ScalarFFClient& ffc)
{
    return DoFitInternal(ffc,0); //No contraint.
}
template <class T> double FittedFunctionImp<T>::DoFit(const DensityFFClient& ffc)
{
    return DoFitInternal(ffc,0); //No contraint.
}

template <class T> double FittedFunctionImp<T>::DoFitInternal(const ScalarFFClient& ffc,double constraint)
{
    SMat Sinv=itsBasisSet->InvOverlap(itsLAParams);
    itsFitCoeff=Sinv*itsBasisSet->Overlap(itsMesh.get(),*ffc.GetScalarFunction());;
    return 0;
}

template <class T> double FittedFunctionImp<T>::DoFitInternal(const DensityFFClient& ffc,double constraint)
{
    SMat Sinv=itsBasisSet->InvRepulsion(itsLAParams);
    itsFitCoeff=Sinv*ffc.GetRepulsion3C(itsBasisSet.get());
    return 0;
}

//---------------------------------------------------------------------------
//
//  Provide Overlap and Repulsion matricies for derived classes.
//
template <class T> typename FittedFunctionImp<T>::Vec FittedFunctionImp<T>::
FitGet2CenterOverlap(const Fit_IBS* bs) const
{
    return itsFitCoeff * (itsBasisSet->Overlap(itsMesh.get(),*bs));
}

template <class T> typename FittedFunctionImp<T>::Vec FittedFunctionImp<T>::
FitGet2CenterRepulsion(const Fit_IBS* bs) const
{
    return itsFitCoeff * (itsBasisSet->Repulsion(*bs));
}

template <class T> typename FittedFunctionImp<T>::SMat FittedFunctionImp<T>::
FitGet3CenterOverlap(const TOrbital_DFT_IBS<double>* bs) const
{
    const std::vector<SMat>& O3=bs->Overlap3C(*itsBasisSet);
    int n=bs->GetNumFunctions();
    SMat J(n,n);
    Fill(J,0.0);
    size_t i=0;
    for (auto c:itsFitCoeff) J+=SMat(c*O3[i++]);
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
        itsFitCoeff *
        itsBasisSet->Overlap(itsMesh.get(),*ffi->itsBasisSet) *
        ffi->itsFitCoeff;
}

template <class T> double FittedFunctionImp<T>::
FitGetRepulsion(const FittedFunctionImp<T>* ffi) const
{
    return
        itsFitCoeff * itsBasisSet->Repulsion(*ffi->itsBasisSet.get()) *
        ffi->itsFitCoeff;
}

template <class T> double FittedFunctionImp<T>::FitGetCharge() const
{
    return itsFitCoeff * itsBasisSet->Charge();
}

//------------------------------------------------------------------------
//
//  Handy utilities for fitted functions.
//
template <class T> void FittedFunctionImp<T>::ShiftOrigin(const RVec3& newCenter)
{
    itsBasisSet.reset(itsBasisSet->Clone(newCenter)); //TOT Clone
}

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
    return Max(fabs(itsFitCoeff - ffi->itsFitCoeff));
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
    return itsFitCoeff * (*itsBasisSet)(r);
}

template <class T> void  FittedFunctionImp<T>::Eval(const Mesh& m, Vec& v) const
{
    v += Vec(itsFitCoeff * (*itsBasisSet)(m));
}

template <class T> typename FittedFunctionImp<T>::RVec3  FittedFunctionImp<T>::Gradient(const RVec3& r) const
{
    Vec3Vec br = itsBasisSet->Gradient(r);
    RVec3 ret(0,0,0);
    typename Vec    ::const_iterator c(itsFitCoeff.begin());
    typename Vec3Vec::const_iterator b(br.begin());
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
