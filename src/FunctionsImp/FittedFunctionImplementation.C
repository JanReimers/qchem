// File: FittedFunctionImplementation.C  General Fitted Function.


#include "FunctionsImp/FittedFunctionImplementation.H"
#include "BasisSet/IntegralDataBase.H"
#include "ChargeDensity/ChargeDensity.H"
#include "BasisSet/TBasisSet.H"
//#include "oml/isnan.h"
#include "oml/vector.h"
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "oml/imp/binio.h"
//#include "oml/vector_io.h"
#include <iostream>
#include <cassert>

//---------------------------------------------------------------------
//
//  Construction zone.  The CDFit flag indicates to fit the electric
//  instead of the charge density.  The only difference in practice
//  is that all overlap integrals are replaced with repulsion integrals.
//
template <class T> FittedFunctionImplementation<T>::
FittedFunctionImplementation(const rc_ptr<BasisSet>& theFitBasisSet, bool CDfit)
    : itsBasisSet(theFitBasisSet)
    , itsFitCoeff(theFitBasisSet->GetNumFunctions())
    , itsCDFitFlag(CDfit)
{
    Fill(itsFitCoeff,0.0);
    itsFitCoeff(1)=1.0/CastBasisSet()->GetDataBase()->GetCharge()(1);
};

template <class T> FittedFunctionImplementation<T>::FittedFunctionImplementation()
    : itsBasisSet (    )
    , itsFitCoeff (    )
    , itsCDFitFlag(true)
{
    Fill(itsFitCoeff,0.0);
};

//---------------------------------------------------------------------
//
//  Implementation stuff for derived constrained fit classes.
//
template <class T> typename FittedFunctionImplementation<T>::SMat FittedFunctionImplementation<T>::
GetInverseOverlap() const
{
//    assert(!isnan(CastBasisSet()->GetDataBase()->GetInverseRepulsion()));
//    assert(!isnan(CastBasisSet()->GetDataBase()->GetInverseOverlap()));
    return itsCDFitFlag ? CastBasisSet()->GetDataBase()->GetInverseRepulsion()
           : CastBasisSet()->GetDataBase()->GetInverseOverlap();
}

template <class T> void FittedFunctionImplementation<T>::SetFitCoeff(const Vec& fc)
{
    assert(!isnan(fc));
    itsFitCoeff=fc;
}

//--------------------------------------------------------------------------
//
//  Implement all DoFit functions.  The overlaps will be accumlated in
//  itsFitCoeff by the call to InjectRepulsions or InjectOverlaps.
//
#include <typeinfo>
template <class T> double FittedFunctionImplementation<T>::DoFit(const FittedFunctionClient& ffc)
{
    Fill(itsFitCoeff,0.0);
    if (itsCDFitFlag)
        ffc.InjectRepulsions(this,&*itsBasisSet); //itsFitCoeff gets repulsion vector.
    else
        ffc.InjectOverlaps  (this,&*itsBasisSet); //itsFitCoeff gets overlap vector.
//	cout << "DoFit " <<  typeid(*this).name() << " " << typeid(ffc).name() << " " << itsFitCoeff << std::endl;
    return DoFit(0,itsFitCoeff);
}

template <class T> double FittedFunctionImplementation<T>::DoFit(double constraint, const FittedFunctionClient& ffc)
{
    Fill(itsFitCoeff,0.0);
    if (itsCDFitFlag)
        ffc.InjectRepulsions(this,&*itsBasisSet); //itsFitCoeff gets repulsion vector.
    else
        ffc.InjectOverlaps  (this,&*itsBasisSet); //itsFitCoeff gets overlap vector.
//	cout << "DoFit constrained" << constraint << " " << itsFitCoeff << std::endl;
    return DoFit(constraint,itsFitCoeff);
}

template <class T> double FittedFunctionImplementation<T>::DoFit(double,  const Vec& overlap)
{
//	cout << "fit overlap" << overlap << std::endl;
    assert(!isnan(overlap));
    SetFitCoeff(GetInverseOverlap()*overlap);
//	cout << "Fit " << itsFitCoeff << std::endl;
    return 0;
}

//---------------------------------------------------------------------------
//
//  Provide Overlap and Repulsion matricies for derived classes.
//
template <class T> typename FittedFunctionImplementation<T>::Vec FittedFunctionImplementation<T>::
FitGet2CenterOverlap(const BasisSet* bs) const
{
    const TBasisSet<T>* tbs=dynamic_cast<const TBasisSet<T>*>(bs);
    assert(tbs);
    return itsFitCoeff * (CastBasisSet()->GetDataBase()->GetOverlap(*tbs));
}

template <class T> typename FittedFunctionImplementation<T>::Vec FittedFunctionImplementation<T>::
FitGet2CenterRepulsion(const BasisSet* bs) const
{
    const TBasisSet<T>* tbs=dynamic_cast<const TBasisSet<T>*>(bs);
    assert(tbs);
    return itsFitCoeff * (CastBasisSet()->GetDataBase()->GetRepulsion(*tbs));
}

//-------------------------------------------------------------------------------------
//
//  Get overlap and repulsion with a charge density.  And total charge.
//  Again only for derived classes.
//
template <class T> double FittedFunctionImplementation<T>::
FitGetOverlap(const FittedFunctionImplementation<T>* ffi) const
{
    return
        itsFitCoeff *
        CastBasisSet()->GetDataBase()->GetOverlap(*ffi->CastBasisSet()) *
        ffi->itsFitCoeff;
}

template <class T> double FittedFunctionImplementation<T>::
FitGetRepulsion(const FittedFunctionImplementation<T>* ffi) const
{
    return
        itsFitCoeff *
        CastBasisSet()->GetDataBase()->GetRepulsion(*ffi->CastBasisSet()) *
        ffi->itsFitCoeff;
}

template <class T> double FittedFunctionImplementation<T>::FitGetCharge() const
{
    return itsFitCoeff * CastBasisSet()->GetDataBase()->GetCharge();
}

//------------------------------------------------------------------------
//
//  Handy utilities for fitted functions.
//
template <class T> void FittedFunctionImplementation<T>::ShiftOrigin(const RVec3& newCenter)
{
    itsBasisSet.reset(itsBasisSet->Clone(newCenter)); //TOT Clone
}

template <class T> void FittedFunctionImplementation<T>::FitMixIn(const FittedFunction& ff,double c)
{
    const FittedFunctionImplementation<T>* ffi = dynamic_cast<const FittedFunctionImplementation<T>*>(&ff);
    assert(ffi);
    assert(itsBasisSet->GetID() == ffi->itsBasisSet->GetID());
    itsFitCoeff = itsFitCoeff*(1-c) + ffi->itsFitCoeff*c;
}

template <class T> double FittedFunctionImplementation<T>::FitGetChangeFrom(const FittedFunction& ff) const
{
    const FittedFunctionImplementation<T>* ffi = dynamic_cast<const FittedFunctionImplementation<T>*>(&ff);
    assert(ffi);
    assert(itsBasisSet->GetID() == ffi->itsBasisSet->GetID());
    return Max(fabs(itsFitCoeff - ffi->itsFitCoeff));
}

template <class T> void FittedFunctionImplementation<T>::ReScale(double factor)
{
    itsFitCoeff*=factor;
}

//-------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> double  FittedFunctionImplementation<T>::operator()(const RVec3& r) const
{
    return itsFitCoeff * (*CastBasisSet())(r);
}

template <class T> void  FittedFunctionImplementation<T>::Eval(const Mesh& m, Vec& v) const
{
    v += Vec(itsFitCoeff * (*CastBasisSet())(m));
}

template <class T> typename FittedFunctionImplementation<T>::RVec3  FittedFunctionImplementation<T>::Gradient(const RVec3& r) const
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
template <class T> std::ostream& FittedFunctionImplementation<T>::Write(std::ostream& os) const
{
    if (Binary())
        BinaryWrite(itsCDFitFlag,os);
    else if (Ascii())
        os << itsCDFitFlag << " ";
    else
        os << "Fit Function: " << std::endl << "  Fit flag = " << itsCDFitFlag << std::endl;

    os << *itsBasisSet;
    if (!Pretty())
    os << itsFitCoeff;
    else
    {
        os << std::endl;
        os << "  Coeff=" << itsFitCoeff << std::endl;
    }

    return os;
}

template <class T> std::istream& FittedFunctionImplementation<T>::Read (std::istream& is)
{
    if (Binary())
        BinaryRead(itsCDFitFlag,is);
    else
        is >> itsCDFitFlag;

    itsBasisSet.reset(BasisSet::Factory(is));
    is >> *itsBasisSet >> itsFitCoeff;
    return is;
}

template class FittedFunctionImplementation<double>;
