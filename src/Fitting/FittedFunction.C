// File: FittedFunctionImp.C  Common imp for Fitted Functions.

#include "Imp/Fitting/FittedFunction.H"
#include <IntegralDataBase.H>
#include <ChargeDensity.H>
#include <BasisSet.H>
#include <Mesh.H>
#include "oml/vector.h"
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "oml/imp/binio.h"
#include <iostream>
#include <cassert>

//---------------------------------------------------------------------
//
//  Construction zone.  The CDFit flag indicates to fit the electric
//  instead of the charge density.  The only difference in practice
//  is that all overlap integrals are replaced with repulsion integrals.
//
template <class T> FittedFunctionImp<T>::
FittedFunctionImp(bs_t& theFitBasisSet,mesh_t& m)
    : itsBasisSet(theFitBasisSet)
    , itsFitCoeff(theFitBasisSet->GetNumFunctions())
    , itsMesh    (m)
    , itsLAParams({qchem::Lapack,qchem::SVD,1e-10,1e-12})
{
    assert(itsMesh);
    Fill(itsFitCoeff,0.0);
    //itsFitCoeff(1)=1.0/CastBasisSet()->GetCharge()(1);
    itsFitCoeff(1)=1.0/itsBasisSet->Charge()(1);
    itsInvOvlp=itsBasisSet->InvOverlap(itsLAParams);
    itsInvRepl=itsBasisSet->InvRepulsion(itsLAParams);
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
    itsFitCoeff=itsInvOvlp*itsBasisSet->GetOverlap(&*itsMesh,ffc.GetScalarFunction());;
    return 0;
}

template <class T> double FittedFunctionImp<T>::DoFitInternal(const DensityFFClient& ffc,double constraint)
{
    itsFitCoeff=itsInvRepl*ffc.GetRepulsion3C(&*itsBasisSet);
    return 0;
}

//---------------------------------------------------------------------------
//
//  Provide Overlap and Repulsion matricies for derived classes.
//
template <class T> typename FittedFunctionImp<T>::Vec FittedFunctionImp<T>::
FitGet2CenterOverlap(const IrrepBasisSet* bs) const
{
    return itsFitCoeff * (itsBasisSet->GetOverlap(&*itsMesh,bs));
}

template <class T> typename FittedFunctionImp<T>::Vec FittedFunctionImp<T>::
FitGet2CenterRepulsion(const IrrepBasisSet* bs) const
{
    return itsFitCoeff * (itsBasisSet->Repulsion(*bs));
}

template <class T> typename FittedFunctionImp<T>::SMat FittedFunctionImp<T>::
FitGet3CenterOverlap(const IrrepBasisSet* bs) const
{
    const std::vector<SMat>& O3=bs->GetOverlap3C(itsBasisSet.get());
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
        itsBasisSet->GetOverlap(&*itsMesh,ffi->itsBasisSet.get()) *
        ffi->itsFitCoeff;
}

template <class T> double FittedFunctionImp<T>::
FitGetRepulsion(const FittedFunctionImp<T>* ffi) const
{
    const IrrepBasisSet* bs=ffi->itsBasisSet.get();
    return
        itsFitCoeff * itsBasisSet->Repulsion(*bs) *
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
    return itsFitCoeff * (*CastBasisSet())(r);
}

template <class T> void  FittedFunctionImp<T>::Eval(const Mesh& m, Vec& v) const
{
    v += Vec(itsFitCoeff * (*CastBasisSet())(m));
}

template <class T> typename FittedFunctionImp<T>::RVec3  FittedFunctionImp<T>::Gradient(const RVec3& r) const
{
    Vec3Vec br = CastBasisSet()->Gradient(r);
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
//    if (StreamableObject::Binary())
//       
//    else if (StreamableObject::Ascii())
//        
//    else
        os << "Fit Function: " << std::endl;

    os << *itsBasisSet;
    if (!StreamableObject::Pretty())
    os << itsFitCoeff;
    else
    {
        os << std::endl;
        os << "  Coeff=" << itsFitCoeff << std::endl;
    }

    return os;
}

template <class T> std::istream& FittedFunctionImp<T>::Read (std::istream& is)
{
//    if (StreamableObject::Binary())
//        BinaryRead(itsCDFitFlag,is);
//    else
//        is >> itsCDFitFlag;

    //itsBasisSet.reset(IrrepBasisSet::Factory(is));
    //is >> *itsBasisSet >> itsFitCoeff;
    return is;
}

template class FittedFunctionImp<double>;
